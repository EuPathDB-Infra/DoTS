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

    $self->moveToGus($dbh, $taxonId, $genomeId, $gdgTab, $gdtTab, $testNum);
    return "finished moving $gdgTab and $gdtTab into GUS central dogma tables";
}

#----------------
#
# sub-routines
#
#----------------

sub moveToGus {
    my ($self, $dbh, $taxonId, $genomeId, $gdgTab, $gdtTab, $testNum) = @_;

    my $sql = "select genome_dots_gene_id from $gdgTab"
	. " where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $sql .= " and rownum <= $testNum" if $testNum;

    my $sth = $dbh->prepareAndExecute($sql);
    my @gdg_ids;
    while (my ($id) = $sth->fetchrow_array()) { push @gdg_ids, $id; }
    $sth->finish();

    my $gi_cat = $self->getGeneInstanceCategory;
    my $fake_gene = $self->getFakeGene; # place holder for gDGs w/o sDG ids
    my @chrs = DoTS::Gene::Util::getChromsInfo($dbh, $genomeId);
    my %chrIds; foreach (@chrs) { $chrIds{$_->{chr}} => $_->{chr_id}; }

    foreach my $gdg_id (@gdg_ids) {
	$sql = "select gdg.*, gdt.na_sequence_id, gdt.blat_alignment_id,"
	    . " gdt.genome_dots_transcript_id"
	    . " from $gdgTab gdg, $gdtTab gdt "
	    . " where gdg.genome_dots_gene_id = gdt.genome_dots_gene_id"
	    . " and gdg.genome_dots_gene_id = $gdg_id";
	$sth = $dbh->prepareAndExecute($sql);
	my @rows;
	while (my $row = $sth->fetchrow_hashref('NAMElc')) {  push @rows, $row; }
	die "unexpected error: not info about gDG.$gdg_id" unless @rows;
	my $r1 = $rows[0];

	my $gene = $fake_gene;
	if ($r1->{gene_id}) {
	    $gene = GUS::Model::DoTS::Gene->new({'gene_id'=>$r1->{gene_id}}); 
	}
	my $gf = GUS::Model::DoTS::GeneFeature->new();
	$gf->setName('gDG.' . $r1->{genome_dots_gene_id});
	$gf->setNumberOfExons($r1->{number_of_exons});
	$gf->setScore($r1->{confidence_score});
	$gf->setNaSequenceId($chrIds{$r1->{chromosome}});

	my $gi = GUS::Model::DoTS::GeneInstance->new();
	$gi->addChild($gene);
	$gi->addChild($gi_cat);
	$gi->addChild($gf);
	my $success = $gi->sumbit(1);
	die "could not submit new GeneInstance for " . $gf->getName if !$success;

	my $gfLoc = GUS::Model::DoTS::NALocation->new();
	$gfLoc->addChild($gf);
	$gfLoc->setStartMin($r1->{chromosome_start});
	$gfLoc->setEndMax($r1->{chromosome_end});
	$gfLoc->setIsReversed($r1->{strand} eq '-' ? 1 : 0);
	$success = $gfLoc->sumbit(0);
	die "could not submit new NALocation for " . $gf->getName if !$success;

	my $exoncount = $gf->getNumberOfExons;
	my @exonstarts = split(/,/, @{ $r1->{exonstarts} });
	my @exonends = split(/,/, @{ $r1->{exonends} });
	for (my $i=1; $i<=$exoncount; $i++) {
	    my $ef = GUS::Model::DoTS::ExonFeature->new();
	    $ef->setParentId($gf->getGeneFeatureId);
	    $ef->setName('');
	    $ef->setNaSequenceId($gf->getNaSequenceId);
	    $ef->setOrderNumber($i);
	    my $efLoc = GUS::Model::DoTS::NALocation->new();
	    $efLoc->addChild($ef);
	    $efLoc->setStartMin($exonstarts[$i-1]);
	    $efLoc->setEndMax($exonends[$i-1]);
	    $efLoc->setIsReversed($gf->getIsReversed);
	    $success = $efLoc->sumbit(1);
	    die "could not submit new NALocation and ExonFeature for "
		. $gf->getName . " exon $i" if !$success;
	}

	foreach (@rows) {
	    my $dt = $_->{na_sequence_id}; 
	    # my $blat = $_->{blat_alignment_id};
	    my $rf = GUS::Model::DoTS::RnaFeature->new();
	    $rf->setParentId($gf->getGeneFeatureId);
	    $rf->setName('');
	    $rf->setNaSequenceId($dt);
	    # TODO: make NALocation for coords within dt
	    # TODO: make RNAInstance, associate with RNAInstanceCategory & RNA
	    # TODO: make RNAFeatureExon that associate exons with RNAs
	    # (not yet supported because the current gDG creation process looses this info.)
	    $success = $rf->sumbit(0);
	    die "could not submit new RNAFeature for " . $gf->getName . " DT.$dt" if !$success;
	}

    }
}

sub getGeneInstanceCategory {
    my ($self) = @_;

    my $gic = GUS::Model::DoTS::GeneInstanceCategory->new({ 'gene_instance_category_id' => 1 });
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
    $gene->setReviewStatusId(0);
    $gene->submit;
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
