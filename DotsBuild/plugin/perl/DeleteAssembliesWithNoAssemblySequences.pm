package DoTS::DotsBuild::Plugin::DeleteAssembliesWithNoAssemblySequences;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::PluginMgr::Plugin;
use GUS::ObjRelP::DbiDatabase;
use GUS::Model::DoTS::Assembly;

my $argsDeclaration =
[
    integerArg({
        name => 'testnumber',
        descr => 'number of iterations for testing',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    integerArg({
        name => 'taxon_id',
        descr => 'taxon_id of assemblies that will be deleted',
        constraintFunc => undef,
        reqd => 1,
        isList => 0
    }),

];

my $purposeBrief = <<PURPOSEBRIEF;
Deletes assemblies and children that have no AssemblySequences
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
Deletes assemblies and children that have no AssemblySequences
PLUGIN_PURPOSE

#check the documentation for this
my $tablesAffected = [
    ['DoTS::Assembly', 'and children']    
];

my $tablesDependedOn = [
    ['DoTS::Assembly', ''],
];

my $howToRestart = <<PLUGIN_RESTART;
PLUGIN_RESTART

my $failureCases = <<PLUGIN_FAILURE_CASES;
Fails if there is no DoTS.SequenceType entry for 'mRNA'
PLUGIN_FAILURE_CASES

my $notes = <<PLUGIN_NOTES;
PLUGIN_NOTES


my $documentation = {
    purposeBrief => $purposeBrief,
    purpose => $purpose,
    tablesAffected => $tablesAffected,
    tablesDependedOn => $tablesDependedOn,
    howToRestart => $howToRestart,
    failureCases => $failureCases,
    notes => $notes
};

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);
      
  $self->initialize({requiredDbVersion => 3.5,
           cvsRevision => '$Revision$', # cvs fills this in!
           name => ref($self),
           argsDeclaration   => $argsDeclaration,
           documentation     => $documentation
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

