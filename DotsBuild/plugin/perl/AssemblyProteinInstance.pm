#################################################################################################
##AssemblyProteinInstance.pm 
##
##This is a ga plug_in to populate the ProteinInstance table with entries that correspond  
##to the assemblies which generally have entries in RNAFeature,RNAInstance,RNA,Protein,and 
##TranslatedAAFeature linking them to entries in the RNA table.  
##
##Created Aug 7, 2001
##
##
##Deborah Pinney 
##
##algorithm_id=4689       
##algorithm_imp_id=5524
##################################################################################################
package AssemblyProteinInstance;

use strict;

use Objects::GUSdev::RNA;
use Objects::GUSdev::Assembly;
use Objects::GUSdev::TranslatedAAFeature;
use Objects::GUSdev::ProteinSequence;##comment this out or delete 
#use Objects::GUSdev::ProteinInstance;##uncomment this
use Objects::GUSdev::RNAFeature;
use Objects::GUSdev::Protein;
#use Objects::GUSdev::RNAInstance;##uncomment
use Objects::GUSdev::RNASequence; ##comment this outor delete
use Objects::GUSdev::TranslatedAASequence;


my $Cfg;
my $ctx;
my $count=0;
my $debug = 0;
my %idHash;
$| = 1;

sub new {
	my $Class = shift;
	$Cfg = shift;
	return bless {}, $Class;
}

sub Usage {
	my $M   = shift;
	return 'Plug_in to populate the ProteinInstance table relating Protein to TranslatedAAFeature for the assemblies';
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
                  h => 'the taxon_id of the assemblies to be used',
                 }
	}
}

my $algoInvo;

sub Run {
  my $M   = shift;
  $ctx = shift;
  
  $algoInvo = $ctx->{self_inv};
  
  print STDERR $ctx->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print STDERR "Testing on $ctx->{'cla'}->{'testnumber'}\n" if $ctx->{'cla'}->{'testnumber'};
  
  unless ($ctx->{'cla'}->{'taxon_id'}) {
    die "you must provide a taxon_id\n";
  }
  
  my $dbh = $ctx->{'self_inv'}->getQueryHandle();
  
  my $time = `date`;
  chomp($time);
  print STDERR ("Starting entries: $time\n");
  
 #get the na_sequence_ids of all the assemblies with the appropriate taxon_id 
  my @ids = &getIds($dbh);

  foreach my $id (@ids){ 
    my $ass = Assembly->new({'na_sequence_id' => $id});
    if( $ass->retrieveFromDB()){
      my $rnafeat = $ass->getChild('RNAFeature',1) ? $ass->getChild('RNAFeature') : &makeRNAFeature($ass);
      #my $rnainst = $rnafeat->getChild('RNAInstance',1) ? $rnafeat->getChild('RNAInstance') : &makeRNAInstance($rnafeat);**uncomment after schema change
      my $rnainst = $rnafeat->getChild('RNASequence',1) ? $rnafeat->getChild('RNASequence') : &makeRNAInstance($rnafeat);#delete after schema change
      my $rna = $rnainst->getParent('RNA',1) ? $rnainst->getParent('RNA') : &makeRNA($rnainst);
      $ass->addToSubmitList($rna);
      my $prot = $rna->getChild('Protein',1) ? $rna->getChild('Protein') : &makeProtein($rna);
      #my $protInst = $prot->getChild('ProteinInstance',1) ? $prot->getChild('ProteinInstance') : &makeProteinInstance($prot);**uncomment after schema c
      my $protInst = $prot->getChild('ProteinSequence',1) ? $prot->getChild('ProteinSequence') : &makeProteinInstance($prot);#delete after schema change
      my $trAF = $rnafeat->getChild('TranslatedAAFeature',1) ? $rnafeat->getChild('TranslatedAAFeature') : &makeTranAAFeat($rnafeat);
      $trAF->addChild($protInst);
      my $trAS = $trAF->getParent('TranslatedAASequence',1) ?  $trAF->getParent('TranslatedAASequence',1) : &makeTranAASeq($trAF);
      $ass->addToSubmitList($trAS);
      $count += $ass->submit();
      $ass->undefPointerCache();
      print STDERR ("Processed $count ending with na_sequence_id: $id\n") if $count % 1000 == 0;
    }
  }

  $time = `date`;
  chomp($time);
  print STDERR ("$count entries have been made to the ProteinInstance table: $time\n");
  return "$count entries have been made to the ProteinInstance table in GUS\n";
} 

#subroutine that makes an array with the assembly na_sequence_ids-avoids assemblies with existing proteinsequence entries because-important for restart 

sub getIds {
    my ($dbh) = @_;
    my @ids;
    my $i=0;
    my $st1 = $dbh->prepareAndExecute("select  /** RULE */ a.na_sequence_id from assembly a where a.taxon_id = $ctx->{'cla'}->{'taxon_id'} and a.description != 'DELETED'");
    
    my $st2 = $dbh->prepare("select /** RULE */ p.protein_sequence_id from proteinsequence p, translatedaafeature f, rnafeature r where r.na_sequence_id = ? and r.na_feature_id = f.na_feature_id and f.aa_feature_id = p.aa_feature_id");  
    
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
    return @ids;
}

#subroutine that makes an RNAFeature entry corresponding to an assembly-this does not need an alternate
sub makeRNAFeature {
    my ($ass) = @_;
    my $name = 'assembly';
    my $rnafeat = RNAFeature->new({'name'=>$name});
    $rnafeat->setParent($ass);
    return $rnafeat;
}

##subroutine that makes an RNAInstance entry and sets its RNAFeature child*****uncomment this after the schema changes
#sub makeRNAInstance {
#  my ($rnaf) = @_;
#  my $is_reference = 0;
#  my $revigetDatabase()->ew_status_id = ?;  ####not created yet-get this
#  my $rna_instance_category_id) = ?; ####not created yet-get this
#  my %attHash = ('review_status_id'=>$review_status_id, 'is_reference'=>$is_reference, 'rna_instance_category_id'=>$rna_instance_category_id);
#  my $rnainst = RNAInstance->new(\%attHash);
#  $rnainst->setParent($rnaf);
#  return $rnainst;
#}

#subroutine that makes an RNAInstance entry and sets its RNAFeature child*****this is the alternative-delete or comment out after the schema changes
sub makeRNAInstance {
  my ($rnaf) = @_;
  my $is_reference = 0;
  my $manually_reviewed = 0;
  my $rna_sequence_type_id = 0;
  my %attHash = ('manually_reviewed'=>$manually_reviewed, 'is_reference'=>$is_reference, 'rna_sequence_type_id'=>$rna_sequence_type_id);
  my $rnainst = RNASequence->new(\%attHash);
  $rnainst->setParent($rnaf);
  return $rnainst;
}

##subroutine that makes an RNA entry and adds its RNAInstance child*****uncomment this after the schema changes
#sub makeRNA {
#  my ($rnai) = @_;
#  my $review_status_id = ?;  ####not created yet-get this
#  my $rna = RNA->new({'review_status_id'=>$review_status_id});
#  $rna->addChild($rnai);
#  return $rna;
#}

#subroutine that makes an RNA entry and adds its RNAInstance child*****this is the alternative - delete or comment out after the schema changes
sub makeRNA {
  my ($rnai) = @_;
  my $manually_reviewed = 0;
  my $is_reference = 0;
  my $name = 'assembly';
  my $rna = RNA->new({'manually_reviewed'=>$manually_reviewed, 'is_reference'=>$is_reference, 'name'=>$name});
  $rna->addChild($rnai);
  return $rna;
}

##subroutine that makes a Protein entry and sets its RNA parent******uncomment this after the schema change
#sub makeProtein {
#  my ($rna) = @_;
#  my $review_status_id = ?;  ####not created yet-get this
#  my $protein = Protein->new({'review_status_id'=>$review_status_id});
#  $protein->setParent($rna);
#  return $protein;
#}

#subroutine that makes a Protein entry and sets its RNA parent******this is the alternate-delete or comment out after the schema change
sub makeProtein {
  my ($rna) = @_;
  my $is_reference = 0;
  my $manually_reviewed = 0;
  my $protein = Protein->new({'is_reference'=>$is_reference, 'manually_reviewed' => $manually_reviewed });
  $protein->setParent($rna);
  return $protein;
}

##subroutine that makes a ProteinInstance entry and sets its Protein parent*****uncomment this after the schema change
#sub makeProteinInstance {
#  my ($prot) = @_;
#  my $review_status_id = ?;  ####not created yet-get this
#  my $is_reference = 0;
#  my %attHash = ('review_status_id'=>$review_status_id, 'is_reference'=>$is_reference);
#  my $proteininst = ProteinInstance->new(\%attHash);
#  $proteininst->setParent($prot);
#  return $proteininst;
#}

#subroutine that makes a ProteinInstance entry and sets its Protein parent*****this is the alternate-delete or comment out after schema change
sub makeProteinInstance {
  my ($prot) = @_;
  my $is_reference = 0;
  my $manually_reviewed = 0;
  my %attHash = ('manually_reviewed'=>$manually_reviewed, 'is_reference'=>$is_reference);
  my $proteininst = ProteinSequence->new(\%attHash);
  $proteininst->setParent($prot);
  return $proteininst;
}

#subroutine that makes a TranslatedAAFeature entry and sets its RNAFeature parent*****don't need alternate to test
sub makeTranAAFeat {
  my($rnaf) = @_;
  my $is_predicted = 0;
  my $manually_reviewed = 0;
  my $TranAAFeat = TranslatedAAFeature->new({'manually_reviewed'=>$manually_reviewed, 'is_predicted'=>$is_predicted});
  $TranAAFeat->setParent($rnaf);
  return $TranAAFeat;
}

#subroutine that makes a TranslatedAASequence entry and adds its TranslatedAAFeature child*****don't need alternate to test
sub makeTranAASeq {
  my $transaafeat = @_;
  my $sub_view = 'TranslatedAASequence';
  my $seq_ver = 0;
  my %attHash = ('subclass_view'=> $sub_view, 'sequence_version'=>$seq_ver);
  my $tranaaseq = TranslatedAASequence->new(\%attHash);
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

