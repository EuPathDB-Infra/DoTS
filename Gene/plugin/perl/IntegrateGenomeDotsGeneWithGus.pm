package DoTS::Gene::Plugin::IntegrateGenomeDotsGeneWithGus;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;
use DoTS::Gene::Util;

use GUS::PluginMgr::Plugin;

use GUS::Model::DoTS::Gene;
use GUS::Model::DoTS::GeneInstance;
use GUS::Model::DoTS::GeneInstanceCategory;
use GUS::Model::DoTS::GeneFeature;
use GUS::Model::DoTS::ExonFeature;

use GUS::Model::DoTS::RNA;
use GUS::Model::DoTS::RNAInstance;
use GUS::Model::DoTS::RNAInstanceCategory;
use GUS::Model::DoTS::RNAFeature;

use GUS::Model::DoTS::RNAFeatureExon;
use GUS::Model::DoTS::NALocation;

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
		     }),

	 integerArg({name => 'test_number',
		     descr => 'number of entries to do test on',
		     constraintFunc=> undef,
		     reqd  => 0,
		     isList => 0
		     }),

	 booleanArg({name => 'is_restart',
		     descr => 'whether this is a restart (no effect if only_delete is set)',
		     constraintFunc=> undef,
		     reqd  => 0,
		     isList => 0
		     }),

	 booleanArg({name => 'only_delete',
		     descr => 'flag, only delete exist gDG result',
		     constraintFunc=> undef,
		     reqd  => 0,
		     isList => 0
		     }),

	 booleanArg({name => 'only_insert',
		     descr => 'flag, only insert new gDG result',
		     constraintFunc=> undef,
		     reqd  => 0,
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
    my $tempLogin = $self->getArg('temp_login');
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');
    my $gdgTab = $tempLogin . '.' . $self->getArg('genome_dots_gene_cache');
    my $gdtTab = $tempLogin . '.' . $self->getArg('genome_dots_transcript_cache');
    my $testNum = $self->getArg('test_number');
    my $isRestart = $self->getArg('is_restart');
    my $onlyDelete = $self->getArg('only_delete');
    my $onlyInsert = $self->getArg('only_insert');

    my ($gis, $gfs, $gdgs) = $self->getExistGiGfGdg($dbh, $genomeId);
    unless ($isRestart || $onlyInsert) {
	my $efs = $self->getExistEfOrRf($dbh, 'DoTS.ExonFeature', $genomeId);
	my $rfs = $self->getExistEfOrRf($dbh, 'DoTS.RnaFeature', $genomeId);
	$self->cleanOldResults($dbh, $gis, $gfs, $efs, $rfs);
    }

    if ($onlyDelete) {
	return "finished deleting old gDG results"; 
    }
    
    $self->moveToGus($dbh, $taxonId, $genomeId, $gdgTab, $gdtTab, $gdgs, $testNum);
    return "finished moving $gdgTab and $gdtTab into GUS central dogma tables";
}

#----------------
#
# sub-routines
#
#----------------

sub getExistGiGfGdg {
    my ($self, $dbh, $genomeId) = @_;

    my $sql = "select gi.gene_instance_id, gf.na_feature_id, gf.name"
	. " from DoTS.GeneFeature gf, DoTS.GeneInstance gi"
	. " where gf.na_feature_id = gi.na_feature_id and gi.gene_instance_category_id = 1"
	. " and gf.external_database_release_id = $genomeId";
    $self->log("finding existing GeneInstance & GeneFeature for gDGs: $sql");
    my $sth = $dbh->prepareAndExecute($sql);
    my ($gis, $gfs, $gdgs) = ([], [], {});
    while (my ($gi, $gf, $nam) = $sth->fetchrow_array) { 
	push @$gis, $gi;
	push @$gfs, $gf;
	my $id = $1 if $nam =~ /gdg\.(\d+)/i;
	$gdgs->{$id} = 1 if defined $id;
    }
    $self->log("found " . scalar(@$gis) . " geneinstances and genefeatures");
    return ($gis, $gfs, $gdgs);
}

sub getExistEfOrRf {
    my ($self, $dbh, $tab, $genomeId) = @_;

    my $sql = "select f.na_feature_id"
	. " from $tab f, DoTS.GeneFeature gf, DoTS.GeneInstance gi"
	. " where f.parent_id = gf.na_feature_id"
	. " and gf.na_feature_id = gi.na_feature_id"
	. " and gi.gene_instance_category_id = 1"
	. " and gf.external_database_release_id = $genomeId";
    $self->log("finding existing $tab for gDGs: $sql");
    my $sth = $dbh->prepareAndExecute($sql);
    my $res = [];
    while (my ($f) = $sth->fetchrow_array) { 
	push @$res, $f;
    }
    $self->log("found " . scalar(@$res) . " $tab");
    return $res;
}

sub cleanOldResults {
    my ($self, $dbh, $gis, $gfs, $efs, $rfs) = @_;

    my $c = 0;
    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.GeneInstance', 'gene_instance_id', $gis);
    $self->log("deleted $c GeneInstance entries by gene_instance_id");

    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.NALocation', 'na_feature_id', $gfs);
    $self->log("deleted $c NALocation entries by na_feature_id of GeneFeature");

    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.NALocation', 'na_feature_id', $efs);
    $self->log("deleted $c NALocation entries by na_feature_id of ExonFeature");

    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.ExonFeature', 'na_feature_id', $efs);
    $self->log("deleted $c ExonFeature entries by na_feature_id");

    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.RnaFeature', 'na_feature_id', $rfs);
    $self->log("deleted $c RnaFeature entries by na_feature_id");

    $c = &DoTS::Gene::Util::delete($dbh, 'DoTS.GeneFeature', 'na_feature_id', $gfs);
    $self->log("deleted $c GeneFeature entries by na_feature_id");
}

sub moveToGus {
    my ($self, $dbh, $taxonId, $genomeId, $gdgTab, $gdtTab, $skip, $testNum) = @_;

    my $sql = "select genome_dots_gene_id from $gdgTab"
	. " where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $sql .= " and rownum <= $testNum" if $testNum;

    $self->log("finding gDGs to move: $sql");
    my $sth = $dbh->prepareAndExecute($sql);
    my @gdg_ids;
    while (my ($id) = $sth->fetchrow_array()) { push @gdg_ids, $id; }
    $sth->finish();

    my $gi_cat = $self->getGeneInstanceCategory;
    my $fake_gene = $self->getFakeGene; # place holder for gDGs w/o sDG ids
    my $chrIds = DoTS::Gene::Util::getChromToId($dbh, $genomeId);

    my $tally = 0;
    my $total = scalar(@gdg_ids);
    foreach my $gdg_id (@gdg_ids) {
	$tally++;
	next if $skip && $skip->{$gdg_id};
	$sql = "select gdg.gene_id, gdg.genome_dots_gene_id, gdg.number_of_exons, "
	    . " gdg.chromosome, gdg.chromosome_start, gdg.chromosome_end, gdg.strand, "
	    . " gdg.exonstarts, gdg.exonends, "
	    . " gdt.na_sequence_id, gdt.blat_alignment_id, gdt.genome_dots_transcript_id"
	    . " from $gdgTab gdg, $gdtTab gdt "
	    . " where gdg.genome_dots_gene_id = gdt.genome_dots_gene_id"
	    . " and gdg.genome_dots_gene_id = $gdg_id";
	$sth = $dbh->prepareAndExecute($sql);
	my @rows;
	while (my $row = $sth->fetchrow_hashref('NAME_lc')) {  push @rows, $row; }
	die "unexpected error: no info about gDG.$gdg_id" unless @rows;
	my $r1 = $rows[0];

	my $gene = $fake_gene;
	my $geneId = $r1->{gene_id};
	if ($geneId) {
	    $gene = GUS::Model::DoTS::Gene->new({'gene_id'=>$geneId });
	    unless ($gene->retrieveFromDB()) {
		$self->log("WARNING: geneId $geneId no longer in DoTS.Gene");
		# for now, create gene instance pointing to place holder
		$gene = $fake_gene;
	    }
	}
	my $gf = GUS::Model::DoTS::GeneFeature->new();
	$gf->setName('gDG.' . $r1->{genome_dots_gene_id});
	$gf->setNumberOfExons($r1->{number_of_exons});
	$gf->setScore($r1->{confidence_score});
	$gf->setExternalDatabaseReleaseId($genomeId);
	my $cid = $chrIds->{$r1->{chromosome}};
	die "no chrom id found for chr" . $r1->{chromosome} unless $cid;
	$gf->setNaSequenceId($cid);
	my $success = $gf->submit(0);
	$self->error("could not submit new GeneFeature " . $gf->getName) if !$success;

	my $gi = GUS::Model::DoTS::GeneInstance->new();
	$gi->setGeneId($gene->getGeneId);
	$gi->setGeneInstanceCategoryId($gi_cat->getGeneInstanceCategoryId);
	$gi->setNaFeatureId($gf->getNaFeatureId);
	$gi->setIsReference(0);
	$gi->setReviewStatusId(0);
	$success = $gi->submit(0);
	$self->error("could not submit new GeneInstance for " . $gf->getName) if !$success;

	my $gfLoc = GUS::Model::DoTS::NALocation->new();
	$gfLoc->setNaFeatureId($gf->getNaFeatureId);
	$gfLoc->setStartMin($r1->{chromosome_start});
	$gfLoc->setEndMax($r1->{chromosome_end});
	$gfLoc->setIsReversed($r1->{strand} eq '-' ? 1 : 0);
	$success = $gfLoc->submit(0);
	$self->error("could not submit new NALocation for " . $gf->getName) if !$success;
	$gfLoc->undefPointerCache();

	my $exoncount = $gf->getNumberOfExons;
	my @exonstarts = split(/,/, $r1->{exonstarts});
	my @exonends = split(/,/, $r1->{exonends});
	for (my $i=1; $i<=$exoncount; $i++) {
	    my $ef = GUS::Model::DoTS::ExonFeature->new();
	    $ef->setParent($gf);
	    $ef->setName("E$i");
	    $ef->setNaSequenceId($gf->getNaSequenceId);
	    $ef->setOrderNumber($i);
	    $ef->setExternalDatabaseReleaseId($genomeId);
	    $success = $ef->submit(0);
	    $self->error("could not submit new ExonFeature for " . $gf->getName . " exon $i") if !$success;
	    

	    my $efLoc = GUS::Model::DoTS::NALocation->new();
	    $efLoc->setNaFeatureId($ef->getNaFeatureId);
	    my ($s, $e) = ($exonstarts[$i-1], $exonends[$i-1]);
	    $self->error("unexpected exon start and end ($s, $e)") unless defined $s && $e >= $s;
	    $efLoc->setStartMin($s);
	    $efLoc->setEndMax($e);
	    $efLoc->setIsReversed($gfLoc->getIsReversed);
	    $success = $efLoc->submit(0);
	    $self->error("could not submit new NALocation for " . $gf->getName . " exon $i") if !$success;
	    $efLoc->undefPointerCache();
	    $ef->undefPointerCache();
	}

	foreach (@rows) {
	    my $dt = $_->{na_sequence_id}; 
	    my $blat = $_->{blat_alignment_id};
	    my $rf = GUS::Model::DoTS::RNAFeature->new();
	    $rf->setParent($gf);
	    $rf->setName("BLAT.$blat");
	    $rf->setNaSequenceId($dt);
	    $rf->setExternalDatabaseReleaseId($genomeId);
	    # TODO: make NALocation for coords within dt
	    # TODO: make RNAInstance, associate with RNAInstanceCategory & RNA
	    # TODO: make RNAFeatureExon that associate exons with RNAs
	    # (not yet supported because the current gDG creation process looses this info.)
	    $success = $rf->submit(0);
	    $self->error("could not submit new RNAFeature for " . $gf->getName . " DT.$dt") if !$success;
	    $rf->undefPointerCache();
	}

	$gf->undefPointerCache();
	$gi->undefPointerCache();

	$self->log("integrated $tally of $total gDGs") unless $tally % 200;
    }
    $self->log("integrated $tally of $total gDGs");
}

sub getGeneInstanceCategory {
    my ($self) = @_;

    my $gic = GUS::Model::DoTS::GeneInstanceCategory->new({ 'gene_instance_category_id' => 1 });
    return $gic if $gic->retrieveFromDB();

    $gic->setName('gDG');
    $gic->setDescription('DoTS gene created from genome alignments of DTs');
    if ($gic->submit()) {
	$self->log("made/found GeneInstanceCategory " . $gic->getName . '(id=1)');
	return $gic;
    } else { 
	die "not able to make/find GeneInstanceCategory " . $gic->getName . '(id=1)';
    }
}

sub getFakeGene {
    my ($self) = @_;
    my $gene = GUS::Model::DoTS::Gene->new({'gene_id'=>-1});
    return $gene if $gene->retrieveFromDB();

    $gene->setReviewStatusId(0);
    my $success = $gene->submit;
    die "not able to submit fake gene (id: -1)" unless $success;

    return $gene;
}

sub addGene {
  my ($self,$rnaArray) = @_;
  my $count;
  foreach my $rna_id (@$rnaArray) {
    my $rna = GUS::Model::DoTS::RNA->new({'rna_id' => $rna_id});
    $rna->retrieveFromDB();
    $count += $self->makeGene($rna);
    $self->log ("$count : rna_id : $rna_id\n");
    $rna->undefPointerCache();
  }
  $self->log ("$count total gene rows inserted and rna rows updated\n");  
}

sub makeGene {
  my ($self) = shift;
  my ($rna) = @_;
  my $review_status_id = 0;  
  my $gene = GUS::Model::DoTS::Gene->new({'review_status_id'=>$review_status_id});
  $gene->addChild($rna);
  my $num = $gene->submit();
  return $num;
}

sub _getNotes {
    my $notes =<<NOTES;
To be filled in.
NOTES
    $notes;
}

1;
