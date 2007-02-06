package DoTS::Gene::Plugin::MoveGenomeDotsGeneToAllgenes;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;
use DoTS::Gene::Util;

use GUS::PluginMgr::Plugin;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'Move genome-based DoTS Genes from temp space to Allgenes';
  
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
		     }),

     stringArg({name => 'copy_table_suffix',
		descr => 'suffix of names for tables to keep a local copy of result, relieve next build from archiving',
		constraintFunc=> undef,
		isList => 0, 
		reqd => 0,
		})
	 ];

    $self->initialize({requiredDbVersion => 3.5,
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
    my $tempLogin = $self->getArg('temp_login');
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');
    my $gdgTab = $tempLogin . '.' . $self->getArg('genome_dots_gene_cache');
    my $gdtTab = $tempLogin . '.' . $self->getArg('genome_dots_transcript_cache');
    my $copyTabSuffix = $self->getArg('copy_table_suffix');

    if ($copyTabSuffix) {
	$self->log("making a copy of results to relieve next build from archiving");
	$self->makeCopy($dbh, $gdgTab, 'gDG', $gdtTab, 'gDT', $copyTabSuffix);
    }

    $self->moveToAllgenes($dbh, $tempLogin, $taxonId, $genomeId, $gdgTab, $gdtTab);

    return "finished moving $gdgTab and $gdtTab into allgenes";
}

#----------------
#
# sub-routines
#
#----------------

sub makeCopy {
    my ($self, $dbh, $gDGtab, $gDGpref, $gDTtab, $gDTpref, $copyTabSuf) = @_;

    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');

    my $sql = "create table $gDGpref$copyTabSuf as "
    	. "(select * from $gDGtab where taxon_id = $taxonId"
	. " and genome_external_db_release_id = $genomeId)";
    $self->log("making a copy of gene result: sql=$sql");
    $dbh->sqlexec($sql);

    my $sql = "create table $gDTpref$copyTabSuf as "
    	. "(select * from $gDTtab where taxon_id = $taxonId"
	. " and genome_external_db_release_id = $genomeId)";
    $self->log("making a copy of transcript result: sql=$sql");
    $dbh->sqlexec($sql);
}

sub moveToAllgenes {
    my ($self, $dbh, $tempLogin, $taxonId, $genomeId, $gdgTab, $gdtTab) = @_;

    my $aid = &DoTS::Gene::Util::getAnalysisId($dbh, $taxonId, $genomeId);

    $self->deleteOldResults($dbh, $tempLogin, $aid);

    my ($maxg, $maxt) = $self->getMaxUsedIds($dbh);

    $self->moveNewResults($dbh, $aid, $taxonId, $genomeId, $gdgTab, $gdtTab, $maxg, $maxt);
}

sub moveNewResults {
    my ($self, $dbh, $aid, $taxonId, $genomeId, $gdgTab, $gdtTab, $maxg, $maxt) = @_;

    $maxg = $maxg ? $maxg : 0;

    $maxt = $maxt ? $maxt : 0;

    my $sql = "insert into Allgenes.AlignedGene "
	. "(select (genome_dots_gene_id + $maxg) as aligned_gene_id, "
	. " $aid as aligned_gene_analysis_id, '' as aligned_gene_accession, "
	. "CHROMOSOME, CHROMOSOME_START, CHROMOSOME_END, STRAND, GENE_SIZE, "
	. "NUMBER_OF_EXONS, MIN_EXON, MAX_EXON, MIN_INTRON, MAX_INTRON, "
	. "EXONSTARTS, EXONENDS, DEPRECATED, EST_PLOT_SCORE, NUMBER_OF_SPLICE_SIGNALS, "
	. "HAS_HUMAN_MOUSE_ORTHOLOGY, POLYA_SIGNAL_TYPE, HAS_POLYA_TRACK, "
	. "NUMBER_OF_EST_LIBRARIES, NUMBER_OF_EST_CLONES, NUMBER_OF_EST_P53PAIRS, "
	. "CONFIDENCE_SCORE, NUMBER_OF_RNAS, CONTAINS_MRNA, MAX_ORF_LENGTH, "
	. "MAX_ORF_SCORE, MIN_ORF_PVAL, GENE_ID "
	. "from $gdgTab where taxon_id = $taxonId "
	. "and genome_external_db_release_id = $genomeId)";
    $self->log("moving new gene results to allgenes tables: sql=$sql");
    $dbh->sqlexec($sql);

    $sql = "insert into Allgenes.AlignedGeneAssembly "
	. "(select (genome_dots_transcript_id + $maxt) as aligned_gene_assembly_id, "
	. "(genome_dots_gene_id + $maxg) as aligned_gene_id, na_sequence_id, "
	. "'' as cluster_id, blat_alignment_id "
	. "from $gdtTab where taxon_id = $taxonId "
	. "and genome_external_db_release_id = $genomeId)";
    $self->log("moving new transcript results: sql=$sql");
    $dbh->sqlexec($sql);
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

sub deleteOldResults {
    my ($self, $dbh, $tempLogin, $aid) = @_;

    my $sql = "delete Allgenes.AlignedGeneAssembly where aligned_gene_id in "
	. "(select aligned_gene_id from Allgenes.AlignedGene where aligned_gene_analysis_id = $aid)";
    $self->log("deleting old transcripts: sql=$sql");
    $dbh->sqlexec($sql);

    $sql = "delete Allgenes.AlignedGene where aligned_gene_analysis_id = $aid";
    $self->log("deleting old genes: sql=$sql");
    $dbh->sqlexec($sql);
}

sub _getNotes {
    my $notes =<<NOTES;
To be filled in.
NOTES
    $notes;
}

1;
