package DoTS::DotsBuild::Plugin::DeleteAssembliesWithNoAssemblySequences;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::ObjRelP::DbiDatabase;
use GUS::Model::DoTS::Assembly;

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);
  
  my $usage = 'Deletes assemblies and children that have no AssemblySequences';
  
  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'taxon_id',
      t => 'int',
      h => 'taxon_id of assemblies that will be deleted',
     }
    ];
  
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


$| = 1;

sub run {
  my $self   = shift;
  
  print $self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on ".$self->getArgs()->{'testnumber'}."\n" if $self->getArgs()->{'testnumber'};

  my $dbh = $self->getQueryHandle();

  my $taxon_id = $self->getArgs()->{'taxon_id'};

  my $algoInvo = $self->getAlgInvocation;

  $algoInvo->setGlobalDeleteEvidenceOnDelete(0);
  $algoInvo->setGlobalDeleteSimilarityOnDelete(0);

  my $query = "select na_sequence_id from dots.Assembly where taxon_id = $taxon_id MINUS select assembly_na_sequence_id from dots.AssemblySequence"; 
  my @ids;
  my $stmt = $dbh->prepareAndExecute($query);
  while(my ($id) = $stmt->fetchrow_array()){
    push(@ids,$id);
  }
  print STDERR "Deleting ",scalar(@ids), " Assemblies without AssemblySequence children\n";
  
  my $ct = 0;
  my $ctDel = 0;
  foreach my $na_sequence_id (@ids){
    print STDERR "Deleting Assembly $na_sequence_id\n";
    $ct++;
    last if $self->getArgs()->{'testnumber'} && $ct > $self->getArgs()->{testnumber};
    my $ass = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $na_sequence_id});
    if($ass->retrieveFromDB()){
      &markAssemblyAndChildrenDeleted($ass);
      $ctDel += $ass->submit();
      $ass->undefPointerCache();
    }else{
      print STDERR "Unable to retrieve entry for Assembly $na_sequence_id from DB\n";
    }
  }

  my $ret = "Deleted $ctDel Assemblies\n";

  print "$ret";

  return "$ret";
}


sub markAssemblyAndChildrenDeleted {
  my($delAss) = @_;
  foreach my $ts ( $delAss->getTranslatedAASequences(1,1)){
    $delAss->markFromAASequenceToSelfDeleted($ts);
    $delAss->addToSubmitList($ts);
  }
  my $rna = $delAss->getRNA(1,1);
  if($rna){
    $delAss->markFromRNAToSelfDeleted();
    $delAss->addToSubmitList($rna);
  }
  $delAss->retrieveAllChildrenExceptAssSeqsFromDB(1);
  $delAss->markDeleted(1);
}

1;

