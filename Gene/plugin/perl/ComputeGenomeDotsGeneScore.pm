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
;
    my $failureCases = <<FAILURE_CASES;
Not known yet.
FAILURE_CASES
;
    my $notes = &_getNotes;

    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
;
    my $argsDeclaration  =
	[
	 stringArg({name => 'skip_chrs',
		    descr => 'comma separated chromosomes for which to skip process (for restart)',
		    constraintFunc=> undef,
		    reqd  => 0,
		    isList => 1, 
		}),
    
	 integerArg({name => 'subset_selector',
		     descr => '1 for spliced or contains_mrna, -1 for others, 0 for do not care',
		     constraintFunc=> sub { die "must be -1, 0, or 1 if defined, but is $_[0]"
						if defined($_[0]) && $_[0] !~ /^-1|0|1$/; },
		     reqd  => 0,
		     default => 0,
		     isList => 0 
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
    my @done_chrs; push @done_chrs, @$skip_chrs;

    my $c = 0;
    foreach my $coord (@$coords) {
        print "processing chr" . $coord->{chr} . ":" . $coord->{start} . "-" . $coord->{end} . "\n"; 
	$c += $self->processRegion($dbh, $coord, $args);
	$dbh->commit;
	$self->getDb->undefPointerCache();
	push @done_chrs, $coord->{chr};
	print "completed/skipped chromosomoes: " . join(',', @done_chrs) . "\n";
    }
    return "finished gDG score calculation for $c genome-based dots genes";
}

#----------------
#
# sub-routines
#
#----------------
sub processRegion {
    my ($self, $dbh, $coord, $args) = @_;

    my $gdgids = &selectGenomeDotsGenes(@_);
    my $gdg_tab = $args->{'temp_login'} . '.' . $args->{'genome_dots_gene_cache'};
    my $gdt_tab = $args->{'temp_login'} . '.' . $args->{'genome_dots_transcript_cache'};
    my $chr_id = $coord->{'chr_id'};
    my $bs_tab = $args->{'temp_login'} . '.' . $args->{'blat_signals_cache'};

    my $stats = {};
    my $tot = scalar(@$gdgids);
    my $c = 0;
    for (my $i=0; $i<$tot; $i++) {
	my $gid = $gdgids->[$i];
	$self->log("processing GDG.$gid ...") if $self->getArg('debug');
	
	# confidence score info
	my $cs_info = {};
	# 1. splice signal, polya signal, polya track
	$self->log("\tsignals ...") if $self->getArg('debug');
	DoTS::Gene::ConfidenceScore::Signals::setSignals($dbh, $gdt_tab, $bs_tab,
							 $gid, $cs_info);
	# 2. input seq composition (mRNA, total RNA, EST clones, EST libs, EST 5-3 pairs)
	$self->log("\tcomposition ...") if $self->getArg('debug');
	DoTS::Gene::ConfidenceScore::Composition::setComposition($dbh, $gdg_tab, $gdt_tab,
								 $gid, $cs_info);
	# 3. protein coding potential
	$self->log("\tcoding potential ...") if $self->getArg('debug');
	DoTS::Gene::ConfidenceScore::Coding::setCoding($dbh, $gdg_tab, $gdt_tab,
						       $gid, $cs_info);
	# 4. EST 5'-3' plot score
        # $self->log("\tEST plot score ...") if $self->getArg('debug');
	# DoTS::Gene::ConfidenceScore::ESTPlotScore::setESTPlotScore($dbh, $gdg_tab, $gdt_tab,
	#					       $gid, $chr_id, $cs_info);
	# overall score
	$self->log("\toverall score ...") if $self->getArg('debug');
	DoTS::Gene::ConfidenceScore::Overall::setScore($dbh, $gdg_tab,
						       $gid, $cs_info, $stats);
	$self->log("\tdb update ...") if $self->getArg('debug');
	&dbUpdate($dbh, $gdg_tab, $gid, $cs_info);

	unless (++$c % 100) {
	    $self->log("processed $c");
	    $dbh->commit;
	}
    }
    $c;
}

sub selectGenomeDotsGenes {
    my ($self, $dbh, $coord, $args) = @_;

    my $genomeId = $args->{'genome_db_rls_id'};
    my $tempLogin = $args->{'temp_login'};
    my $taxonId = $args->{'taxon_id'};
    my $ss_sel = $args->{'subset_selector'};
    my $gdgCache = $args->{'genome_dots_gene_cache'};
    my $isRerun = ($args->{'skip_chrs'} ? 1 : 0);

    my ($start, $end) = ($coord->{start}, $coord->{end});

    my $sql = "select genome_dots_gene_id from $gdgCache "
            . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	    . "and chromosome = '" . $coord->{chr} . "' "
	    . ($start ? "and chromosome_end >= $start " : "")
	    . ($end ? "and chromosome_start <= $end " : "")
	    . ($isRerun ? "and confidence_score is null" : "");
    if ($ss_sel) {
	$sql .= " and " . ($ss_sel == 1 ? "(max_intron >= 47 or contains_mrna = 1)"
                                        : "max_intron < 47 and (contains_mrna = 0 or contains_mrna is null)"); 
    }
    print "running sql: $sql ...\n";
    my $sth = $dbh->prepareAndExecute($sql);
    my @gdgids;
    while (my ($gdgid) = $sth->fetchrow_array) { push @gdgids, $gdgid; }
    print "found " . scalar(@gdgids) . " genome dots genes to process\n";
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

sub _getNotes {
    my $notes =<<NOTES;
Calculates confidence scores of genome-based DoTS Genes as follows:

Look at these signals for aligned genes:
-- splice signal: GT-AG (in constituent DoTS assemblies)
-- polyA signal type (AATAAA or ATTAAA and the reverse complements)
-- polyA track
(assume FindGenomeSignals plugin has been run)

Look at EST composition of aligned genes:
-- number of RNAs
-- whether contains_mRNA
-- number of distinct EST clones
-- number of EST libraries
-- number of 5\'-3\' pairs of ESTs from the same clones

Look at protein coding potential (FrameFinder/DIANA_ATG predictions)
-- max (FrameFinder) ORF length (among all DTs)
-- min (FrameFinder) ORF p_value
-- max DIANA_ATG score

Look at EST stacking consistency
-- EST plot score (5\' ends of 5\' ESTs and 3\' ends of 3\' ESTs)

NOTES
    $notes;
}

1;
