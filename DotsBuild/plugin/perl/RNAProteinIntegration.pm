#################################################################################################
##RNAProteinIntegration.pm 
##
##This is a ga plug_in to populate the RNAFeature,RNAInstance,ProteinFeature,and ProteinInstance 
##tables with entries that correspond  to mRNA in assemblies. Presently this is done for the taxon_ids 
##specified on the command line. It is written specifically for entries in ExternalNASequence from 
##EMBL-22, GenBank-78 (including refseq-992) 
##
##Created Aug 8, 2001
##
##
##Deborah Pinney 
##
##algorithm_id=5590       
##algorithm_imp_id=7125
##################################################################################################
package RNAProteinIntegration;

use strict;

use Objects::GUSdev::ExternalNASequence;
use Objects::GUSdev::RNAFeature;
use Objects::GUSdev::RNASequence;   ##comment this out or delete it 
#use Objects::GUSdev::RNAInstance; ##uncomment this
use Objects::GUSdev::TranslatedAAFeature;
use Objects::GUSdev::ProteinSequence; ##comment this out or delete it
#use Objects::GUSdev::ProteinInstance; ##uncomment this
use Objects::GUSdev::RNA;
use Objects::GUSdev::Protein;



my $Cfg;
my $ctx;
my $count=0;
my $debug = 0;
my $algoInvo;
my $dbh;
$| = 1;

sub new {
	my $Class = shift;
	$Cfg = shift;
	return bless {}, $Class;
}

sub Usage {
	my $M   = shift;
	return 'Plug_in to populate the RNAFeature, RNAInstance,TranslatedAAFeatureProteinInstance table relating Protein to TranslatedAAFeature for the assemblies';
}


sub EasyCspOptions {
    my $M   = shift;
	{

  testnumber  => {
                  o => 'testnumber=i',
									h => 'number of iterations for testing',
								 },
  taxon_id   => {
                  o => 'taxon_id=i',
                  h => 'the taxon_id of externalnasequence',
                 }
	}
}

sub Run {
    my $M   = shift;
    $ctx = shift;
    
    $algoInvo = $ctx->{self_inv};
    
    print STDERR $ctx->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
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
	my $extNAseq = ExternalNASequence->new({'na_sequence_id' => $id});
	$extNAseq->retrieveFromDB();
	my $rnafeat = $extNAseq->getChild('RNAFeature',1) ? $extNAseq->getChild('RNAFeature') : &makeRNAFeature ($extNAseq);
	my $rnaseq = $rnafeat->getChild('RNASequence',1) ? $rnafeat->getChild('RNASequence') : &makeRNASequence ($rnafeat);
	my $rna; 
	next unless ($rna = $rnaseq->getParent('RNA',1) ? $rnaseq->getParent('RNA') : &getRNA($id, $dbh, $rnaseq));
	$extNAseq->addToSubmitList($rna);
	my $prot;
	next unless ($prot = $rna->getChild('Protein',1));
	my $transaafeat;
	next unless $transaafeat = &makeTransAAFeat($id, $dbh, $rnafeat);

	my $protseq = $transaafeat->getChild('ProteinSequence', 1) ? $transaafeat->getChild('ProteinSequence') : &makeProteinSequence($transaafeat);

	$protseq->setParent($prot);
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
  my $st = $dbh->prepareAndExecute("select na_sequence_id from externalnasequence where taxon_id in ($ctx->{'cla'}->{'taxon_id'}) and sequence_type_id = 2 and external_db_id in (992,22,78)"); 
  
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
  my $external_db_id = $extNAseq->get('external_db_id');
  my $source_id = $extNAseq->get('source_id');
  my $name;
  if ($external_db_id == 992){
    $name = "REFSEQ";
  }
  else {
    $name = "mRNA";
  }
  my $is_predicted = 0;
  my $manually_reviewed = 0;
  
  my %attHash =('name'=>$name, 'is_predicted'=>$is_predicted, 'manually_reviewed'=>$manually_reviewed, 'source_id'=>$source_id, 'external_db_id'=>$external_db_id ); 
  my $newRNAFeature = RNAFeature->new(\%attHash);

  $newRNAFeature->setParent($extNAseq);

  return $newRNAFeature;
}

#subroutine that makes an entry into RNASequence for the mRNA  
sub makeRNASequence {
  my ($rnafeat) = @_;
  my $is_reference = 0;
  my $manually_reviewed =0; #delete when review_status_id is available
  #my $review_status_id = ?;  #######get this - not created yet
  my $rna_sequence_type_id = 1; 
  #my $rna_instance_category_id = ?;   #######get this - not created yet
  my %attHash = ('is_reference'=>$is_reference, 'manually_reviewed'=>$manually_reviewed, 'rna_sequence_type_id'=>$rna_sequence_type_id);
  my $newRNASequence = RNASequence->new(\%attHash);
  
  $newRNASequence->setParent($rnafeat);
  
  return $newRNASequence;
}

#identify the rna_id that corresponds to the assembly containing the mRNA whose na_sequence_id = $id
sub getRNA {
  my ($id, $dbh, $rnaseq) = @_;
  my $st = $dbh->prepare("select s.rna_id from rnafeature f, rnasequence s, assemblysequence a where a.na_sequence_id = ? and a.assembly_na_sequence_id = f.na_sequence_id and f.na_feature_id = s.na_feature_id");

  $st->execute($id);
  if(my ($rna_id) = $st->fetchrow_array) {
    my $rna = RNA->new({'rna_id'=>$rna_id});
    $rna->retrieveFromDB();
    $rna->addChild($rnaseq);
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
  my $st1 = $dbh->prepare("select protein_id from transcript where name = 'CDS' and na_sequence_id = ?");
  my $st2 = $dbh->prepare("select aa_sequence_id from nrdbentry where source_id = ?");
  my $is_predicted = 0;
  my $manually_reviewed = 0;
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
  if ($newTranslatedAAFeature = $rnafeat->getChild('TranslatedAAFeature', 1)) {
      if ($newTranslatedAAFeature->get('aa_sequence_id') != $aa_sequence_id) {
	  $newTranslatedAAFeature->set('aa_sequence_id', $aa_sequence_id);
      }
  }
  else {
      my %attHash = ('is_predicted'=>$is_predicted, 'manually_reviewed'=>$manually_reviewed, 'aa_sequence_id'=>$aa_sequence_id);
      $newTranslatedAAFeature = TranslatedAAFeature->new(\%attHash);
  }
  
  $newTranslatedAAFeature->setParent($rnafeat);

  return $newTranslatedAAFeature;
}

#subroutine to make an entry into ProteinSequence that links the TranslatedAAFeature corresponding to the GenBank 
#translation of the mRNA and the Protein entry that corresponds to the RNA for the assembly containing the mRNA
sub makeProteinSequence {
  my ($transaafeat) = @_;
  my $is_reference = 0;
  #my $review_status_id = ?; #######get this
  my $manually_reviewed = 0;
  my %attHash =('is_reference'=>$is_reference, 'manually_reviewed'=>$manually_reviewed);
  my $newProteinSequence = ProteinSequence->new(\%attHash);
  $newProteinSequence->setParent($transaafeat);
  return $newProteinSequence;
}

1;

__END__

=pod
=head1 Description
B<mRNAandProtein> - a plug-in that creates and integrates RNAFeature, RNAFeature, ProteinFeature, and TranslatedAAFeature entries that corresponds to each mRNA.

=head1 Purpose
B<mRNAandProtein> is a 'plug-in' GUS application that makes and integrates RNAFeature, RNAInstance, ProteinInstance, and TranslatedAAFeature entries for each mRNA.
