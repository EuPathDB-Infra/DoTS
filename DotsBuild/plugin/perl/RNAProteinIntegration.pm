package DoTS::DotsBuild::Plugin::RNAProteinIntegration;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use GUS::Model::DoTS::ExternalNASequence;
use GUS::Model::DoTS::RNAFeature;
use GUS::Model::DoTS::RNAInstance;
use GUS::Model::DoTS::TranslatedAAFeature;
use GUS::Model::DoTS::ProteinInstance;
use GUS::Model::DoTS::RNA;
use GUS::Model::DoTS::Protein;



my $ctx;
my $count=0;
my $debug = 0;
my $algoInvo;
my $dbh;

$| = 1;

sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self,$class);
    
    my $usage = 'Plug_in to populate the RNAFeature, RNAInstance,TranslatedAAFeatureProteinInstance table relating Protein to TranslatedAAFeature for the assemblies';
    
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


sub Run {
    my $M   = shift;
    $ctx = shift;
    
    $algoInvo = $ctx->{self_inv};
    
    print STDERR $ctx->{cla}->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
    print STDERR "Testing on $ctx->{'cla'}->{'testnumber'}\n" if $ctx->{'cla'}->{'testnumber'};
    
    unless ($ctx->{'cla'}->{'taxon_id'}) {
	die "you must provide a taxon_id\n";
    }
    
    $dbh = $ctx->{self_inv}->getQueryHandle();
    
    my $time = `date`;
    chomp($time);
    print STDERR ("Starting entries : $time\n");
    
#call the subroutine to make an array of na_sequence_ids for the mRNA  
    my @ids = &getmRNA();

#loop through each mRNA na_sequence_id in @ids  
    foreach my $id (@ids) {
	my $extNAseq = GUS::Model::DoTS::ExternalNASequence->new({'na_sequence_id' => $id});
	$extNAseq->retrieveFromDB();
	my $rnafeat = $extNAseq->getChild('GUS::Model::DoTS::RNAFeature',1) ? $extNAseq->getChild('GUS::Model::DoTS::RNAFeature') : &makeRNAFeature ($extNAseq);
	my $rnainst = $rnafeat->getChild('GUS::Model::DoTS::RNAInstance',1) ? $rnafeat->getChild('GUS::Model::DoTS::RNAInstance') : &makeRNAInstance ($rnafeat);
	my $rna; 
	next unless ($rna = $rnainst->getParent('GUS::Model::DoTS::RNA',1) ? $rnainst->getParent('GUS::Model::DoTS::RNA') : &getRNA($id, $dbh, $rnainst));
	$extNAseq->addToSubmitList($rna);
	my $prot;
	next unless ($prot = $rna->getChild('GUS::Model::DoTS::Protein',1));
	my $transaafeat;
	next unless $transaafeat = &makeTransAAFeat($id, $dbh, $rnafeat);

	my $protinst = $transaafeat->getChild('GUS::Model::DoTS::ProteinInstance', 1) ? $transaafeat->getChild('GUS::Model::DoTS::ProteinInstance') : &makeProteinInstance($transaafeat);

	$protinst->setParent($prot);
	$extNAseq->addToSubmitList($prot);
	$count += $extNAseq->submit();
	$extNAseq->undefPointerCache();
	if ($count%1000==0) {
	    $time = `date`;
	    chomp($time); 
	    print STDERR ("$count    $time\n");
	}   
    }
    $time = `date`;
    chomp($time);
    print STDERR ("Finishing entries: $count completed:  $time\n");
    
    return ("$count mRNA entries have complete RNA/protein integration.\n");
}

#subroutine that gets all the mRNA na_sequence_ids and puts them in @ids
sub getmRNA {
  my @ids;
  my $st = $dbh->prepareAndExecute("select na_sequence_id from dots.externalnasequence where taxon_id in ($ctx->{'cla'}->{'taxon_id'}) and sequence_type_id in (2,7) and external_database_release_id in (select external_database_release_id from sres.externaldatabaserelease where external_database_id in (992,22,78))"); 
  
  while (my ($na_sequence_id) = $st->fetchrow_array) {
    if ( $ctx->{'cla'}->{'testnumber'} && @ids >= $ctx->{'cla'}->{'testnumber'}) {
      last;
    }
    push(@ids,$na_sequence_id); 
  } 
  my $length = @ids;
  print STDERR ("There are $length mRNA ids to process\n");
  return @ids;
}

#subroutine that puts an entry into RNAfeature that represents the mRNA and returns the na_feature_id   
sub makeRNAFeature {
  my ($extNAseq) = @_;
  my $external_database_release_id = $extNAseq->get('external_database_release_id');
  my $source_id = $extNAseq->get('source_id');
  my $name;
  my $newExternalDatabaseRelease = GUS::Model::SRes::ExternalDatabaseRelease->new({'external_database_release_id'=>$external_database_release_id});
  my $external_database_id = $newExternalDatabaseRelease->get('external_database_id');
  if ($external_database_id == 992){
    $name = "REFSEQ";
  }
  else {
    $name = "mRNA";
  }
  my $is_predicted = 0;
  my $review_status_id = 0;
  
  my %attHash =('name'=>$name, 'is_predicted'=>$is_predicted, 'review_status_id'=>$review_status_id, 'source_id'=>$source_id, 'external_database_release_id'=>$external_database_release_id ); 
  my $newRNAFeature = GUS::Model::DoTS::RNAFeature->new(\%attHash);

  $newRNAFeature->setParent($extNAseq);

  return $newRNAFeature;
}

#subroutine that makes an entry into RNAInstance for the mRNA  
sub makeRNAInstance {
  my ($rnafeat) = @_;
  my $is_reference = 0;
  my $review_status_id = 0;
  my $rna_instance_category_id = 1;
  my %attHash = ('is_reference'=>$is_reference, 'review_status_id'=>$review_status_id, 'rna_instance_category_id'=>$rna_instance_category_id);
  my $newRNAInstance = GUS::Model::DoTS::RNAInstance->new(\%attHash);
  
  $newRNAInstance->setParent($rnafeat);
  
  return $newRNAInstance;
}

#identify the rna_id that corresponds to the assembly containing the mRNA whose na_sequence_id = $id
sub getRNA {
  my ($id, $dbh, $rnainst) = @_;
  my $st = $dbh->prepare("select s.rna_id from dots.rnafeature f, dots.rnainstance s, dots.assemblysequence a where a.na_sequence_id = ? and a.assembly_na_sequence_id = f.na_sequence_id and f.na_feature_id = s.na_feature_id");

  $st->execute($id);
  if(my ($rna_id) = $st->fetchrow_array) {
    my $rna = GUS::Model::DoTS::RNA->new({'rna_id'=>$rna_id});
    $rna->retrieveFromDB();
    $rna->addChild($rnainst);
    $st->finish();
    return $rna;
  }
  else {
    return;
  }
}

#subroutine to make a TranslatedAAFeature that will correspond to the RNAFeature entry for this mRNA
#note that protein_id is the source_id from a GenBank entry and not the GUS Protein table id
#updates the TranslatedAAFeature table if aa_sequence_id has changed
sub makeTransAAFeat {
  my ($id, $dbh, $rnafeat) = @_;
  my $st1 = $dbh->prepare("select protein_id from dots.transcript where name = 'CDS' and na_sequence_id = ?");
  my $st2 = $dbh->prepare("select aa_sequence_id from dots.nrdbentry where source_id = ?");
  my $is_predicted = 0;
  my $review_status_id = 0;
  $st1->execute($id);
  my $source_id; 
  unless (($source_id) = $st1->fetchrow_array()) {
    return;
  }
  $st1->finish();
  if ($source_id =~ (/(\S+)\.\d+/)) {
    $source_id = $1;
  }
  $st2->execute($source_id);
  my $aa_sequence_id;
  unless (($aa_sequence_id) = $st2->fetchrow_array()) {
    return;
  }
  $st2->finish(); 
  my $newTranslatedAAFeature; 
  if ($newTranslatedAAFeature = $rnafeat->getChild('GUS::Model::DoTS::TranslatedAAFeature', 1)) {
      if ($newTranslatedAAFeature->get('aa_sequence_id') != $aa_sequence_id) {
	  $newTranslatedAAFeature->set('aa_sequence_id', $aa_sequence_id);
      }
  }
  else {
      my %attHash = ('is_predicted'=>$is_predicted, 'review_status_id'=>$review_status_id, 'aa_sequence_id'=>$aa_sequence_id);
      $newTranslatedAAFeature = GUS::Model::DoTS::TranslatedAAFeature->new(\%attHash);
  }
  
  $newTranslatedAAFeature->setParent($rnafeat);

  return $newTranslatedAAFeature;
}

#subroutine to make an entry into ProteinSequence that links the TranslatedAAFeature corresponding to the GenBank 
#translation of the mRNA and the Protein entry that corresponds to the RNA for the assembly containing the mRNA
sub makeProteinInstance {
  my ($transaafeat) = @_;
  my $is_reference = 0;
  my $protein_instance_category_id = 2;
  my $review_status_id = 0;
  my %attHash =('is_reference'=>$is_reference, 'review_status_id'=>$review_status_id, 'protein_instance_category_id'=>$protein_instance_category_id);
  my $newProteinInstance = GUS::Model::DoTS::ProteinInstance->new(\%attHash);
  $newProteinInstance->setParent($transaafeat);
  return $newProteinInstance;
}

1;

__END__

=pod
=head1 Description
B<mRNAandProtein> - a plug-in that creates and integrates RNAFeature, RNAFeature, ProteinFeature, and TranslatedAAFeature entries that corresponds to each mRNA.

=head1 Purpose
B<mRNAandProtein> is a 'plug-in' GUS application that makes and integrates RNAFeature, RNAInstance, ProteinInstance, and TranslatedAAFeature entries for each mRNA.
