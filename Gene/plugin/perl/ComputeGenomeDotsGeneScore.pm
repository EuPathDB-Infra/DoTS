package DoTS::Gene::Plugin::ComputeGenomeDotsGeneScore;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

use DoTS::Gene::Util;
use DoTS::Gene::ConfidenceScore::Signals;
use DoTS::Gene::ConfidenceScore::Composition;
use DoTS::Gene::ConfidenceScore::Coding;
use DoTS::Gene::ConfidenceScore::ESTPlotScore;
use DoTS::Gene::ConfidenceScore::Overall;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'Calculates quality score for genome-based DoTS Genes';
  
    my $purpose = $purposeBrief;

    my $tablesDependedOn = [];
    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
give comma separated list of completed chroms (from log file) 
RESTART

    my $failureCases = <<FAILURE_CASES;
Not known yet.
FAILURE_CASES

    my $notes = <<NOTES;
=pod

=head1 DESCRIPTION

Calculates confidence scores of I<genome-based DoTS Genes> as follows:

=over 4

=item I)
Look at these signals for aligned genes:

=over 4

=item 1.
splice signal: GT-AG (in constituent DoTS assemblies)

=item 2.
polyA signal type (AATAAA or ATTAAA and the reverse complements)

=item 3.
polyA track

=back

(assume blat_alignment_signals.pl has been run to populate
the BlatAlignmentSignals table)

=item II)
Look at EST composition of aligned genes:

=over 4

=item 1.
number of RNAs

=item 2.
whether contains_mRNA

=item 3.
number of distinct EST clones

=item 4.
number of EST libraries

=item 5.
number of 5'-3' pairs of ESTs from the same clones

=back

=item III)
Look at protein coding potential (FrameFinder/DIANA_ATG predictions)

=over 4

=item 1.
max (FrameFinder) ORF length (among all DTs)

=item 2.
min (FrameFinder) ORF p_value

=item 3.
max DIANA_ATG score

=back

=item IV)
Look at EST stacking consistency

=over 4

=item 1.
EST plot score (5' ends of 5' ESTs and 3' ends of 3' ESTs)

=back

=back

=head1 AUTHOR

Y. Thomas Gan <ygan\@pcbi.upenn.edu>, October 20, 2003

=head1 COPYRIGHT

Copyright, Trustees of University of Pennsylvania 2004. 

=cut

NOTES
    
    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
    my $argsDeclaration  =
    [
     stringArg({name => 'skip_chrs',
		descr => 'comma separated chromosomes for which to skip process (for restart)',
		constraintFunc=> undef,
		reqd  => 0,
		isList => 1, 
	    }),
    
     stringArg({name => 'chr',
		descr => 'chromosome on which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0, 
	    }),

     integerArg({name => 'start',
		descr => 'start of chromosome region in which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0 
		}),

     integerArg({name => 'end',
		descr => 'end of chromosome region in which to do analysis',
		constraintFunc=> undef,
		reqd => 0,
		isList => 0 
		}),

     stringArg({name => 'temp_login',
		descr => 'login for temp table usage',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0, 
	    }),

    stringArg({name => 'genome_dots_gene_cache',
		descr => 'table that caches info about genome dots genes',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'GenomeDotsGene',
		isList => 0 
		}),

    stringArg({name => 'genome_dots_transcript_cache',
		descr => 'table that caches info about genome dots transcripts',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'GenomeDotsTranscript',
		isList => 0 
		}),

    stringArg({name => 'blat_signals_cache',
		descr => 'table that caches info about BLAT alignment signals',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'BlatAlignmentSignals',
		isList => 0 
		}),

    integerArg({name => 'taxon_id',
		descr => 'taxon id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),

    integerArg({name => 'genome_db_rls_id',
		descr => 'genome external database release id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		})
    ];

    $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$',
		     cvsTag => '$Name$',
		     name => ref($self),
		     revisionNotes => '',
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation
		    });
    return $self;
}


sub run {
    my $self = shift;
    
    $self->logAlgInvocationId();
    $self->logArgs();

    my $args = $self->getArgs();
    my $dbh = $self->getQueryHandle();
    my $genomeId = $self->getArg('genome_db_rls_id');

    my ($coords, $skip_chrs) = &DoTS::Gene::Util::getCoordSelectAndSkip($dbh, $genomeId, $args);
    my @done_chrs = @$skip_chrs;

    my $c = 0;
    foreach my $coord (@$coords) {
        $self->log("processing chr" . $coord->{chr} . ":"
		   . $coord->{start} . "-" . $coord->{end}); 
	$c += &processRegion($dbh, $coord, $args);
	$dbh->commit;
	push @done_chrs, $coord->{chr};
	$self->log("completed/skipped chromosomoes: " . join(', ', @done_chrs));
    }
    return "finished gDG score calculation for $c genome-based dots genes";
}

#----------------
#
# sub-routines
#
#----------------
sub processRegion {
    my ($dbh, $coord, $args) = @_;

    my $gdgids = &selectGenomeDotsGenes(@_);
    my $gdg_tab = $args->{'temp_login'} . '.' . $args->{'genome_dots_gene_cache'};
    my $gdt_tab = $args->{'temp_login'} . '.' . $args->{'genome_dots_transcript_cache'};
    my $bs_tab = $args->{'temp_login'} . '.' . $args->{'blat_signals_cache'};

    my $stats = {};
    my $tot = scalar(@$gdgids);
    my $c = 0;
    for (my $i=0; $i<$tot; $i++) {
	my $gid = $gdgids->[$i];
	
	# confidence score info
	my $cs_info = {};
	# 1. splice signal, polya signal, polya track
	DoTS::Gene::ConfidenceScore::Signals::setSignals($dbh, $gdt_tab, $bs_tab,
							 $gid, $cs_info);
	# 2. input seq composition (mRNA, total RNA, EST clones, EST libs, EST 5-3 pairs)
	DoTS::Gene::ConfidenceScore::Composition::setComposition($dbh, $gdg_tab, $gdt_tab,
								 $gid, $cs_info);
	# 3. protein coding potential
	DoTS::Gene::ConfidenceScore::Coding::setCoding($dbh, $gdg_tab, $gdt_tab,
						       $gid, $cs_info);
	# 4. EST 5'-3' plot score
	DoTS::Gene::ConfidenceScore::ESTPlotScore::setESTPlotScore($dbh, $gdg_tab, $gdt_tab,
						       $gid, $cs_info);
	# overall score
	DoTS::Gene::ConfidenceScore::Overall::setScore($dbh, $gdg_tab,
						       $gid, $cs_info, $stats);

	&dbUpdate($dbh, $gdg_tab, $gid, $cs_info);

	unless (++$c % 1000) {
	    print "# processed $c entries\n";
	    $dbh->commit;
	}
    }
    $c;
}

sub selectGenomeDotsGenes {
    my ($dbh, $coord, $args) = @_;

    my $genomeId = $args->{'genome_db_rls_id'};
    my $tempLogin = $args->{'temp_login'};
    my $taxonId = $args->{'taxon_id'};
    my $gdgCache = $args->{'genome_dots_gene_cache'};
    my $isRerun = ($args->{'skip_chrs'} ? 1 : 0);

    my ($start, $end) = ($coord->{start}, $coord->{end});

    my $sql = "select genome_dots_gene_id from $gdgCache "
            . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	    . "and chromosome = '" . $coord->{chr} . "' "
	    . ($start ? "and chromosome_end >= $start " : "")
	    . ($end ? "and chromosome_start <= $end " : "")
	    . ($isRerun ? "and confidence_score is null" : "");
    my $sth = $dbh->prepareAndExecute($sql);
    my @gdgids;
    while (my ($gdgid) = $sth->fetchrow_array) { push @gdgids, $gdgid; }
    return \@gdgids;
}

=item dbUpdate:

update database with scoring attributes

=cut

sub dbUpdate {
    my ($dbh, $tab, $gid, $cs_info) = @_;

    my @set_clauses;
    foreach (keys %$cs_info) { push @set_clauses, $_ . ' = ' . $cs_info->{$_}; }

    my $sql = "update $tab set " . join(', ', @set_clauses) . ' '
	    . "where genome_dots_gene_id = $gid";

    $dbh->sqlexec($sql);
}

1;
