package DoTS::DotsBuild::Plugin::DupAssSeqFromClusterFile;



@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use FileHandle;
use GUS::PluginMgr::Plugin;
use GUS::Model::DoTS::AssemblySequence;

my $argsDeclaration =
[
    stringArg({
        name => 'clusterFile',
        descr => 'name of cluster file for input',
        constraintFunc => undef,
        reqd => 1,
        isList => 0
    }),
    stringArg({
        name => 'newClusterFile',
        descr => 'name of new cluster file with duplicated assembly sequence ids replaced',
        constraintFunc => undef,
        reqd => 1,
        isList => 0
    }),
    booleanArg({
	name => 'restart',
        descr => 'option to restart if plugin stopped during previous run',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    })
];

my $purposeBrief = <<PURPOSEBRIEF;
Duplicate dots.AssemblySequence rows and create a new cluster file
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
Duplicate dots.AssemblySequence rows for assembly sequences used in more than one cluster and create a new cluster file
PLUGIN_PURPOSE

#check the documentation for this
my $tablesAffected = [
    ['DoTS::AssemblySequence', '']
];

my $tablesDependedOn = [
    ['DoTS::AssemblySequence', ''],
];

my $howToRestart = <<PLUGIN_RESTART;
Use restart argument
PLUGIN_RESTART

my $failureCases = <<PLUGIN_FAILURE_CASES;
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
		     cvsRevision => '$Revision$',  # cvs fills this in!
		     name => ref($self),
		     argsDeclaration   => $argsDeclaration,
		     documentation     => $documentation
		    });

  return $self;
}


$| = 1;

sub run {
  my $self = shift;

  my $clusterFile = $self->getArg('clusterFile');

  my $clusterFH = new FileHandle("$clusterFile", "<") || $self->userError ("Can't open $clusterFile for reading");

  my $done = $self->getClustersDone() if $self->getArg('restart');

  my $newClusterFile = $self->getArg('newClusterFile');

  my $newClusterFH = new FileHandle("$newClusterFile", ">>") || $self->userError ("Can't open $newClusterFile for writing");

  my $clustersDone;

  my %assSeqs;

  while (<$clusterFH>) {
    my $clusterLine = $_;
    $clusterLine =~ m/(Cluster_\d+)/;
    $clustersDone++;
    if ($done->{$1}) {
      $self->listAssSeqs($clusterLine,\%assSeqs);
      next();
    }

    my $newClusterLine = $self->processCluster($clusterLine,\%assSeqs);

    $newClusterFH->print ("$newClusterLine");
  }

  return "$clustersDone, Number of clusters processed from $clusterFile and printed to $newClusterFile\n";

}

sub getClustersDone {
  my ($self) = @_;

  my $newClusterFile = $self->getArg('newClusterFile');

  my $newClusterFH = new FileHandle("$newClusterFile", "<") || $self->userError ("Can't open $newClusterFile for reading");

  my %done;

  while(<$newClusterFH>) {
    $_ =~ m/(Cluster_\d+)/;

    $done{$1}=1;
  }

  $newClusterFH->close();

  return \%done;
}

sub listAssSeqs {
  my ($self,$clusterLine,$assSeqs) = @_;

  $clusterLine =~ /Cluster_\d+\s\(\d+\ssequences\):\s\(([\S\s]+)\)/;

  my @seqs = split (/,\s/,$1);

  foreach my $seqId (@seqs){
    $assSeqs->{$seqId}++;
  }
}

sub processCluster {
  my ($self,$clusterLine,$assSeqs) = @_;

  #Example of cluster file line: Cluster_1648 (4 sequences): (3598654, 2030263, 2028495, 1174727)

  $clusterLine =~ /Cluster_\d+\s\(\d+\ssequences\):\s\(([\S\s]+)\)/;

  my @seqs = split (/,\s/,$1);

  foreach my $seqId (@seqs){
    my $newSeqId = $self->makeDupAssSeqRow($seqId) if ($assSeqs->{$seqId});
    $clusterLine =~ s/$seqId/$newSeqId/ if $newSeqId;
    $assSeqs->{$seqId}++;
  }

  return $clusterLine
}

sub makeDupAssSeqRow {
  my ($self, $seqId) = @_;

  my $assSeq = GUS::Model::DoTS::AssemblySequence->new({'assembly_sequence_id' => $seqId});

  $assSeq->retrieveFromDB();

  my %attHsh = ('assembly_na_sequence_id'=>$assSeq->get('assembly_na_sequence_id'),
                'na_sequence_id'=>$assSeq->get('na_sequence_id'),
                'sequence_start'=>$assSeq->get('sequence_start'),
                'sequence_end'=>$assSeq->get('sequence_end'),
                'quality_start'=>$assSeq->get('quality_start'),
                'quality_end'=>$assSeq->get('quality_end'),
                'assembly_offset'=>$assSeq->get('assembly_offset'),
                'assembly_strand'=>$assSeq->get('assembly_strand'),
                'gapped_sequence'=>$assSeq->get('gapped_sequence'),
                'have_processed'=>$assSeq->get('have_processed'),
                'processed_category'=>$assSeq->get('processed_category'));

  my $newAssSeq = GUS::Model::DoTS::AssemblySequence->new(\%attHsh);

  $newAssSeq->submit();

  my $newAssSeqId = $newAssSeq->getId();

  $newAssSeq->undefPointerCache();

  return $newAssSeqId;
}



1;

