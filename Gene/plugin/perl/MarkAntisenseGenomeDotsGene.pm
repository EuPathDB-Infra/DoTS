package DoTS::Gene::Plugin::MarkAntisenseGenomeDotsGene;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

use DoTS::Gene::Util;
use DoTS::Gene::GenomeFeature;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'Find candidate antisense genes and mark in database';
  
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
		}),

    stringArg({name => 'exon_overlap_cutoff',
		descr => 'minimum exonic overlap required to make a call. eg 15bp, 0.25 (25\% of the smaller)',
		constraintFunc=> undef,
		reqd  => 0,
		default => '0.10',
		isList => 0 
		}),
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
    my $tab = $self->getArg('temp_login') . '.' . $self->getArg('genome_dots_gene_cache');

    my ($coords, $skip_chrs) = &DoTS::Gene::Util::getCoordSelectAndSkip($dbh, $genomeId, $args);
    my @done_chrs = @$skip_chrs;

    foreach my $coord (@$coords) {
        print "processing chr" . $coord->{chr} . ":" . $coord->{start} . "-" . $coord->{end} . "\n"; 

	# the plus/minus strand genes
	my $genesp = $self->getGenomeDotsGenes($dbh, $tab, $coord,'+');
	my $genesm = $self->getGenomeDotsGenes($dbh, $tab, $coord,'-');

	# the overlapping pairs
	my $rev_pairs = $self->getReverseOverlapPairs($genesp, $genesm);

	# pick one from each pair to mark, and mark them
	my $anti = &pickAntisense($rev_pairs);
	&markAntisense($dbh, $tab, $rev_pairs, $anti);
	$dbh->commit();
	$self->getDb->undefPointerCache();
	push @done_chrs, $coord->{chr};
	print "completed/skipped chromosomoes: " . join(',', @done_chrs) . "\n";
    }
    return "finished";
}

############### subroutines ###########################

sub getGenomeDotsGenes {
    my ($self, $dbh, $tab, $coord, $strand) = @_;

    my $sql = $self->getTheQuery($tab, $coord, $strand);

    my $sth = $dbh->prepareAndExecute($sql);

    my @genes;
    while (my ($gid, $gs, $starts, $ends, $numass, $polya_type) = $sth->fetchrow_array) {
	my ($s, $e, $exons) = &_makeGeneModel($starts, $ends);
	$polya_type = ($polya_type ? $polya_type : 0);
	my $gene = [$gid, $numass, $gs, $s, $e, $exons, $polya_type];
	push @genes, $gene;
    }
    $sth->finish;
    \@genes;
}

sub getTheQuery {
    my ($self, $tab, $coord, $strand) = @_;

    my $chr = $coord->{chr};
    my $start = $coord->{start};
    my $end = $coord->{end};

    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');

    my $sql = "select genome_dots_gene_id, gene_size, exonstarts, exonends, number_of_rnas, polya_signal_type "
	    . "from $tab where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	    . "and chromosome = '$chr' and strand = '$strand' "
	    . ($start ? "and chromosome_end >= $start " : "")
            . ($end ? "and chromosome_start <= $end " : "")
	    . "order by chromosome_start asc, chromosome_end asc";
    return $sql;
}

sub _makeGeneModel {
    my ($starts, $ends) = @_;
    my @ss = split(/,/, $starts);
    my @es = split(/,/, $ends);
    my $num_exons = scalar(@ss);
    die "ERROR: number of start and end coordinates do not match"
	unless scalar(@es) == $num_exons;
    
    my @gene;
    for (my $i = 0; $i < $num_exons; $i++) {
	push @gene, [$ss[$i], $es[$i]];
    }

    return ($ss[0], $es[$num_exons-1], \@gene);
}


sub markAntisense {
    my ($dbh, $tab, $rev_pairs, $antis) = @_;

    my $tot = scalar(@$antis);
    for (my $i=0; $i<$tot; $i++) {
	my $idx = $antis->[$i];
	my $dep_gid = $rev_pairs->[$i]->[$idx]->[0];

	my $sql = "update $tab set deprecated = 1 "
	        . "where genome_dots_gene_id = $dep_gid";
	$dbh->sqlexec($sql);
    }
}

sub pickAntisense {
    my ($rev_pairs) = @_;

    my @antis;

    foreach (@$rev_pairs) {
	my ($gp, $gm) = @$_;

	my ($idp, $numassp, $gsp, $sp, $ep, $exonsp, $polya_typep, $pop) = @$gp;
	my ($idm, $numassm, $gsm, $sm, $em, $exonsm, $polya_typem, $pom) = @$gm;

	my $numexonsp = scalar(@$exonsp);
	my $numexonsm = scalar(@$exonsm);

	# todo: consider polya_type and pct overlap too
	if ($numassp != $numassm) {
	    push @antis, ($numassp > $numassm ? 1 : 0);
	} elsif ($pop != $pom) {
	    push @antis, ($pop < $pom ? 1 : 0);
	} elsif ($numexonsp != $numexonsm) {
	    push @antis, ($numexonsp > $numexonsm ? 1 : 0);
	} elsif ($gsp != $gsm) {
	    push @antis, ($gsp > $gsm ? 1 : 0);
	} else {
	    my $i = int rand(1);
	    push @antis, $i;
	}
    }

    \@antis;
}

sub getReverseOverlapPairs {
    my ($self, $genesp, $genesm) = @_;

    my $overlap_cutoff = $self->getArg('exon_overlap_cutoff'); 
    my ($cutoff, $by_bases);
    if ($overlap_cutoff =~ /(\d+)bp/) {
	$cutoff = $1;
	$by_bases = 1;
    } else {
	$cutoff = $overlap_cutoff * 100;
    }

    my $tot_gsp = scalar(@$genesp);
    my $tot_gsm = scalar(@$genesm);
    my $lastI = 0;

    my @rev_pairs;
    foreach my $gp (@$genesp) {
	my ($idp, $numassp, $gsp, $sp, $ep, $exonsp) = @$gp;
	for (my $i = $lastI; $i < $tot_gsm; $i++) {
	    my $gm = $genesm->[$i];
	    my ($idm, $numassm, $gsm, $sm, $em, $exonsm) = @$gm;
	    $lastI = $i+1 if $em <= $sp;
	    last if $sm >= $ep;
	    my $overlap = &DoTS::Gene::GenomeFeature::getSpanOverlap($exonsp, $exonsm);
	    my $pct_p = sprintf("%3.1f", 100*$overlap/$gsp);
	    my $pct_m = sprintf("%3.1f", 100*$overlap/$gsm);

	    my $isRev = 0;
	    if ($by_bases) {
		$isRev = ($overlap >= $cutoff);
	    } else {
		$isRev = ($pct_p > $cutoff || $pct_m > $cutoff); 
	    }

	    if ($isRev) {
		if ($self->getArg('debug')) {
		    $self->log("gDG $idp ($gsp bp) & $idm ($gsm bp), $overlap bp exon overlap,"
			       . "this does" . ($isRev ? '' : ' not') . " meet $cutoff cutoff");
		}
		my @ngp = @$gp; push @ngp, $pct_p;
		my @ngm = @$gm; push @ngm, $pct_m;
		push @rev_pairs, [\@ngp, \@ngm];
	    }
	}
    }
    \@rev_pairs;
}

sub _getNotes {
    my $notes =<<NOTES;
To fill in later.
NOTES
;
    $notes;
}

1;
