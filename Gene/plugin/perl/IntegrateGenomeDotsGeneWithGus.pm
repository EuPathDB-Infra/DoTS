package DoTS::Gene::Plugin::IntegrateGenomeDotsGeneWithGus;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'Integrate genome-based DoTS Genes from temp space to GUS';
  
    my $purpose = $purposeBrief;

    my $tablesDependedOn = [];
    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
Can  not restart yet.
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
    my $tempLogin = $self->{'temp_login'};
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');
    my $gdgTab = $tempLogin . '.' . $self->getArg('genome_dots_gene_cache');
    my $gdtTab = $tempLogin . '.' . $self->getArg('genome_dots_transcript_cache');

    # HAC for now, just copy to allgenes schema
    $self->moveToAllgenes($dbh, $tempLogin, $taxonId, $genomeId, $gdgTab, $gdtTab);
    return "finished moving $gdgTab and $gdtTab into allgenes";

    # TODO: really integrate with GUS central dogma tables
    # $self->moveToGus($dbh, $taxonId, $genomeId, $gdgTab, $gdtTab);
}

#----------------
#
# sub-routines
#
#----------------

sub moveToAllgenes {
    my ($self, $dbh, $tempLogin, $taxonId, $genomeId, $gdgTab, $gdtTab) = @_;

    my $aid = $self->getAnalysisId($dbh, $taxonId, $genomeId);
    die;
    $self->archiveOldResults($dbh, $tempLogin, $aid);

    my ($maxg, $maxt) = $self->getMaxUsedIds($dbh);

    $self->moveNewResults($dbh, $aid, $taxonId, $genomeId, $gdgTab, $gdtTab, $maxg, $maxt);
}

sub moveNewResults {
    my ($self, $dbh, $aid, $taxonId, $genomeId, $gdgTab, $gdtTab, $maxg, $maxt) = @_;
    
    my $sql = "insert into Allgenes.AlignedGene "
	. "select $aid as aligned_gene_analysis_id, "
	. "(genome_dots_gene_id + $maxg) as aligned_gene_id, "
	. "CHROMOSOME, CHROMOSOME_START, CHROMOSOME_END, STRAND, GENE_SIZE, "
	. "NUMBER_OF_EXONS, MIN_EXON, MAX_EXON, MIN_INTRON, MAX_INTRON, "
	. "EXONSTARTS, EXONENDS, DEPRECATED, EST_PLOT_SCORE, NUMBER_OF_SPLICE_SIGNALS "
	. "HAS_HUMAN_MOUSE_ORTHOLOGY, POLYA_SIGNAL_TYPE, HAS_POLYA_TRACK, "
	. "NUMBER_OF_EST_LIBRARIES, NUMBER_OF_EST_CLONES, NUMBER_OF_EST_P53PAIRS "
	. "CONFIDENCE_SCORE, NUMBER_OF_RNAS, CONTAINS_MRNA, MAX_ORF_LENGTH, "
	. "MAX_ORF_SCORE, MIN_ORF_PVAL, GENE_ID "
	. "from $gdgTab where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $self->log("moving new gene results: sql=$sql");
    $dbh->sqlexec();

    $sql = "insert into Allgenes.AlignedGeneAssembly "
	. "select (genome_dots_transcript_id + $maxt) as aligned_gene_assembly_id, "
	. "(genome_dots_gene_id + $maxg) as aligned_gene_id, na_sequence_id, "
	. "similarity_dots_gene_id as cluster_id, blat_alignment_id "
	. "from $gdtTab where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $self->log("moving new transcript results: sql=$sql");
    $dbh->sqlexec();
    $self->log("finished moving");
}

sub getMaxUsedIds {
    my ($self, $dbh) = @_;

    my ($sql, $sth, $maxg, $maxt);

    $sql = "select max(aligned_gene_id) from Allgenes.AlignedGene";
    $sth = $dbh->prepareAndExecute($sql);
    $maxg = $sth->fetchrow_array;

    $sql = "select max(aligned_gene_assembly_id) from Allgenes.AlignedGeneAssembly";
    $sth = $dbh->prepareAndExecute($sql);
    $maxt = $sth->fetchrow_array;

    ($maxg, $maxt);
}

sub archiveOldResults {
    my ($self, $dbh, $tempLogin, $aid) = @_;

    my $sql = "create table ${tempLogin}.AlignedGene_$aid as "
	. "(select * from Allgenes.AlignedGene where aligned_gene_analysis_id = $aid)";
    $self->log("saving old genes: sql=$sql");
    $dbh->sqlexec($sql);

    $sql = "create table ${tempLogin}.AlignedGeneAssembly_$aid as "
	. "(select * from Allgenes.AlignedGeneAssembly where aligned_gene_id in "
	. " (select aligned_gene_id from Allgenes.AlignedGene where aligned_gene_analysis_id = $aid))";
    $self->log("saving old transcripts: sql=$sql");
    $dbh->sqlexec($sql);

    $sql = "delete Allgenes.AlignedGeneAssembly where aligned_gene_id in "
	. "(select aligned_gene_id from Allgenes.AlignedGene where aligned_gene_analysis_id = $aid)";
    $self->log("deleting old transcripts: sql=$sql");
    $dbh->sqlexec($sql);

    $sql = "delete Allgenes.AlignedGene where aligned_gene_analysis_id = $aid";
    $self->log("deleting old transcripts: sql=$sql");
    $dbh->sqlexec($sql);
}

sub getAnalysisId {
    my ($self, $dbh, $taxonId, $genomeId) = @_;

    my $sql = "select aligned_gene_analysis_id from Allgenes.AlignedGene "
	. "where parameters like '%--t $taxonId --dbr $genomeId%'";
    my $sth = $dbh->prepareAndExecute($sql);
    my @aids = ();
    while (my $aid = $sth->fetchrow_array) { push @aids, $aid; }
    my $c = scalar(@aids);
    die "expecting one analysis id but found $c" unless $c == 1;
    $self->log("aligned gene analysis id is $aids[0]");

    return $aids[0];
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

sub _getNotes {
    my $notes =<<NOTES;
To be filled in.
NOTES
    $notes;
}

1;
