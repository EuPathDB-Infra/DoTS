package DoTS::DotsBuild::Plugin::AssemblyProteinInstance;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use GUS::Model::DoTS::RNA;
use GUS::Model::DoTS::Assembly;
use GUS::Model::DoTS::TranslatedAAFeature;
use GUS::Model::DoTS::ProteinInstance;
use GUS::Model::DoTS::RNAFeature;
use GUS::Model::DoTS::Protein;
use GUS::Model::DoTS::RNAInstance;
use GUS::Model::DoTS::TranslatedAASequence;


$| = 1;

sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self,$class);
    
    my $usage = 'Plug_in to populate the ProteinInstance table relating Protein to TranslatedAAFeature for the assemblies';
    
    my $easycsp =
	[{o => 'testnumber',
	  t => 'int',
	  h => 'number of iterations for testing',
         },
	 {o => 'taxon_id',
	  t => 'int',
	  h => 'the taxon_id of the assemblies to be used',
         }];
    
    $self->initialize({requiredDbVersion => {},
		       cvsRevision => '$Revision$',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		       name => ref($self),
		       revisionNotes => 'make consistent with GUS 3.0',
		       easyCspOptions => $easycsp,
		       usage => $usage
		       });
    
    return $self;
}


sub run {
    my $self = shift;
    
    $self->log ($self->getCla{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");
    $self->log ("Testing on $self->getCla{'testnumber'}\n" if $self->getCla{'testnumber'});
    
    unless ($self->getCla{'taxon_id'}) {
	die "you must provide a taxon_id\n";
    }
    
    $self->log ("Starting entries\n");
    
    my $ids = $self->getIds();
    
    my $count = $self->processIds($ids);
    
    $self->log ("$count entries have been made to the ProteinInstance table: $time\n");
    return "$count entries to the ProteinInstance table in GUS\n";
}

sub getIds {
    
    my $self = shift;
    
    my @ids;
    my $i=0;
    
    my $dbh = $self->getQueryHandle();
    
    my $st1 = $dbh->prepareAndExecute("select  /** RULE */ a.na_sequence_id from dots.assembly a where a.taxon_id = $ctx->{'cla'}->{'taxon_id'} and a.description != 'DELETED'");
    
    my $st2 = $dbh->prepare("select /** RULE */ p.protein_sequence_id from dots.proteininstance p, dots.translatedaafeature f, dots.rnafeature r where r.na_sequence_id = ? and r.na_feature_id = f.na_feature_id and f.aa_feature_id = p.aa_feature_id");  
    
    while (my $na_sequence_id = $st1->fetchrow_array) {
	if ($ctx->{'cla'}->{'testnumber'} && $i >= $ctx->{'cla'}->{'testnumber'}) {
	    last;
	}
	$st2->execute($na_sequence_id);
	if ($st2->fetchrow_array) {
	    $st2->finish();
	    next;
	}
	else {
	    $ids[$i]=$na_sequence_id;
	    $i++;
	}
    }
    return \@ids;
} 

sub processIds {
    
    my $self = shift;
    
    my $ids = @_;
    
    my $count = 0;
    
    foreach my $id (@$ids){ 
	my $ass = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $id});
	if( $ass->retrieveFromDB()){
	    my $rnafeat = $ass->getChild('RNAFeature',1) ? $ass->getChild('RNAFeature') : $self->makeRNAFeature($ass);
	    my $rnainst = $rnafeat->getChild('RNAInstance',1) ? $rnafeat->getChild('RNAInstance') : $self->makeRNAInstance($rnafeat);
	    my $rna = $rnainst->getParent('RNA',1) ? $rnainst->getParent('RNA') : $self->makeRNA($rnainst);
	    $ass->addToSubmitList($rna);
	    my $prot = $rna->getChild('Protein',1) ? $rna->getChild('Protein') : $self->makeProtein($rna);
	    my $protInst = $prot->getChild('ProteinInstance',1) ? $prot->getChild('ProteinInstance') : $self->makeProteinInstance($prot);
	    
	    my $trAF = $rnafeat->getChild('TranslatedAAFeature',1) ? $rnafeat->getChild('TranslatedAAFeature') : $self->makeTranAAFeat($rnafeat);
	    $trAF->addChild($protInst);
	    my $trAS = $trAF->getParent('TranslatedAASequence',1) ?  $trAF->getParent('TranslatedAASequence',1) : $self->makeTranAASeq($trAF);
	    $ass->addToSubmitList($trAS);
	    $count += $ass->submit();
	    $ass->undefPointerCache();
	    $self->log ("Processed $count ending with na_sequence_id: $id\n") if $count % 1000 == 0;
	}
    }
    
    return $count;
}


sub makeRNAFeature {
    my $self = shift;
    my ($ass) = @_;
    my $name = 'assembly';
    my $rnafeat = GUS::Model::DoTS::RNAFeature->new({'name'=>$name});
    $rnafeat->setParent($ass);
    return $rnafeat;
}

sub makeRNAInstance {
    my $self = shift;
    my ($rnaf) = @_;
    my $is_reference = 0;
    my $review_status_id = 0;  
    my $rna_instance_category_id) = 0; 
my %attHash = ('review_status_id'=>$review_status_id, 'is_reference'=>$is_reference, 'rna_instance_category_id'=>$rna_instance_category_id);
my $rnainst = GUS::Model::DoTS::RNAInstance->new(\%attHash);
$rnainst->setParent($rnaf);
return $rnainst;
}

sub makeRNA {
    my $self = shift;
    my ($rnai) = @_;
    my $review_status_id = 0;  
    my $rna = GUS::Model::DoTS::RNA->new({'review_status_id'=>$review_status_id});
    $rna->addChild($rnai);
    return $rna;
}

sub makeProtein {
    my $self = shift;
    my ($rna) = @_;
    my $review_status_id = 0; 
    my $protein = GUS::Model::DoTS::Protein->new({'review_status_id'=>$review_status_id});
    $protein->setParent($rna);
    return $protein;
}

sub makeProteinInstance {
    my $self = shift;
    my ($prot) = @_;
    my $review_status_id = 0; 
    my $is_reference = 0;
    my %attHash = ('review_status_id'=>$review_status_id, 'is_reference'=>$is_reference);
    my $proteininst = GUS::Model::DoTS::ProteinInstance->new(\%attHash);
    $proteininst->setParent($prot);
    return $proteininst;
}

sub makeTranAAFeat {
    my $self = shift;
    my($rnaf) = @_;
    my $is_predicted = 0;
    my $manually_reviewed = 0;
    my $TranAAFeat = GUS::Model::DoTS::TranslatedAAFeature->new({'manually_reviewed'=>$manually_reviewed, 'is_predicted'=>$is_predicted});
    $TranAAFeat->setParent($rnaf);
    return $TranAAFeat;
}

sub makeTranAASeq {
    my $self = shift;
    my $transaafeat = @_;
    my $sub_view = 'TranslatedAASequence';
    my $seq_ver = 0;
    my %attHash = ('subclass_view'=> $sub_view, 'sequence_version'=>$seq_ver);
    my $tranaaseq = GUS::Model::DoTS::TranslatedAASequence->new(\%attHash);
    $tranaaseq->addChild($transaafeat);
    return $tranaaseq;
}

1;

__END__
    
    =pod
    =head1 Description
B<AssemblyProteinInstance> - a plug-in that creates a ProteinInstance entry that corresponds to each assembly.
    
    =head1 Purpose
    B<AssemblyProteinInstance> is a 'plug-in' GUS application that makes a ProteinInstance entry for each assembly.
    
    
