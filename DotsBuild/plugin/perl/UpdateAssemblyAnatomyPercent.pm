package DoTS::DotsBuild::Plugin::UpdateAssemblyAnatomyPercent;

@ISA = qw(GUS::PluginMgr::Plugin);

use GUS::PluginMgr::Plugin;
use strict;
use GUS::Model::DoTS::AssemblyAnatomyPercent;
use DoTS::DotsBuild::AssemblyAnatomyNode;


$| = 1;

sub new {
  my ($class) = @_;

  my $self = {};
  bless($self,$class);

  my $usage = 'assign library distribution to assemblies';

 my  $argsDeclaration =
    [ integerArg({name => 'testnumber',
		  descr => 'number of iterations for testing',
		  constraintFunc => undef,
		  reqd => 0,
		  isList => 0
		 }),
      integerArg({name => 'taxon_id',
		  descr => 'taxon_id',
		  constraintFunc => undef,
		  reqd => 1,
		  isList => 0
		 }),
      booleanArg({name => 'restart',
                  descr => 'restarts: ignores those assembies already in the AssemblyAnatomyPercent',
                   constraintFunc => undef,
		  reqd => 0
		 }),
      stringArg({name => 'idSQL',
		  descr => 'SQL query that returns Assembly na_sequence_ids, deletes these from AssemblyAnatomyPercent unless >= restart date and then creates new entries for each na_sequence_id',
		  constraintFunc => undef,
		  reqd => 1,
		  isList => 0
		 }),
      booleanArg({name => 'print_tree',
                  descr => 'prints a tree into STDOUT for each DT',
		  constraintFunc => undef,
		  reqd => 0
		 })];


   my $purposeBrief = <<PURPOSEBRIEF;
Assign library distribution to assemblies
PURPOSEBRIEF

  my $purpose = <<PLUGIN_PURPOSE;
Assign library distribution to assemblies
PLUGIN_PURPOSE

  #check the documentation for this
  my $tablesAffected = [['DoTS::AssemblyAnatomyPercent', 'DoTS::AssemblyAnatomyNode']
		       ];

  my $tablesDependedOn = [
			  ['DoTS::AnatomyLibrary']
			 ];

  my $howToRestart = <<PLUGIN_RESTART;
Use boolean argumnet restart
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

  $self->initialize({requiredDbVersion => 3.5,
		     cvsRevision => '$Revision$',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation,
		    });
  return $self;
}



sub run {
  my ($self) = @_;

  $self->logAlgInvocationId;
  $self->logCommit;

  # get args
  my $testnumber = $self->getArg('testnumber');
  my $sql =  $self->getArg('idSQL');
  my $taxonId = $self->getArg('taxon_id');

  my $dbh = $self->getQueryHandle();

  # make anatomy tree
  my ($root, $nodeHash) = $self->makeTree($dbh);

  # place raw EST counts in each node of the tree
  $self->setESTCounts($nodeHash, $dbh, $taxonId);

  # handle$self-> each DT
  print "Testing on $testnumber\n" if $testnumber;

  my $sql = $self->getArg('idSQL');
  $sql .= " and na_sequence_id not in (select na_sequence_id from dots.assemblyanatomypercent where taxon_id = $taxonId)" if $self->getArg('restart');
  my $stmt = $dbh->prepareAndExecute($self->getArg('idSQL'));
  my $count;
  while (my ($dt) = $stmt->fetchrow_array()) {
    die "Error: na_sequence_id '$dt' returned in the result set is not an integer" unless $dt =~ /\d+/;
    $count++;
    print STDERR "Updated $count rows\n" if ($count % 10000) == 0;
    $self->processDT($dt, $nodeHash, $taxonId, $root, $dbh);
    $self->printTree($root) if $self->getArg('print_tree');
    $self->undefPointerCache();
    last if ($testnumber && $count > $testnumber);
  }

  $self->setResultDescr("Updated $count rows");
}

#########################################################################
##################  subroutines  ########################################
#########################################################################


# static method to make a whole tree
# return (rootNode, hashOfNodesByAnatomyId)
sub makeTree {
  my ($self, $dbh) = @_;
 
  my $sql = "select anatomy_id,parent_id from sres.anatomy order by hier_level";

  my $stmt = $dbh->prepareAndExecute($sql) || die "Can't prepareAndExecute  sql: $sql\n";

  my $root;
  my %nodeHash;

  while (my ($anatomy_id, $parent_id) = $stmt-> fetchrow_array()) {
    my $parent = $nodeHash{$parent_id};
    my $node = DoTS::DotsBuild::AssemblyAnatomyNode->new($anatomy_id, $parent);
    $nodeHash{$anatomy_id} = $node;
    if ($parent) {
      $parent->addKid($node);
    } else {
      $root = $node;
    }
  }
  return ($root, \%nodeHash);
}

# set the raw total est counts for each node in the tree
sub setESTCounts {
  my ($self, $nodeHash, $dbh, $taxonId) = @_;

  # issue a query
  my $sql = "select /*+ RULE */ al.anatomy_id, count(e.na_sequence_id) from dots.anatomylibrary al,dots.library l, dots.assemblysequence a, dots.est e where a.assembly_na_sequence_id is not null and a.na_sequence_id = e.na_sequence_id and l.taxon_id = $taxonId and e.library_id = l.library_id and l.dbest_id = al.dbest_library_id group by al.anatomy_id";

  my $stmt = $dbh->prepareAndExecute($sql) || die "Can't prepareAndExecute  sql: $sql\n";

  while (my ($anatomy_id, $count) = $stmt-> fetchrow_array()) {
    $nodeHash->{$anatomy_id}->setESTCount($count);
  }

}

# process a single dt.  
sub processDT {
  my ($self, $dt, $nodeHash,$taxonId, $root, $dbh) = @_;

  # zero out previous DT's junk
  $root->clearDTValues();

  # load this DT's values into the existing anatomy tree
  # return the sum of the effective counts and the sum of the raw counts
  my ($sum_effective, $sum_raw) = $self->loadDT($nodeHash,$dbh,$dt);

  # percolate from bottom up and write out the rows
  $root->percolateAndWrite($dt, $sum_effective, $sum_raw, $taxonId) if $sum_effective > 0;
}

# For a given DT, foreach anatomyID, place in the anatomyId's node:
#   raw est count
#   effective est count
# Also, accumulate sums for those values
# return Sum_effective, Sum_raw
sub loadDT {
  my ($self, $nodeHash, $dbh, $dtId) = @_;

  # issue a query (anatomy_id, count)
  my $sql = "select al.anatomy_id, count(e.est_id) from dots.anatomylibrary al,dots.library l, dots.assemblysequence a, dots.est e where a.assembly_na_sequence_id =$dtId and a.na_sequence_id = e.na_sequence_id and e.library_id = l.library_id and l.dbest_id = al.dbest_library_id group by al.anatomy_id";
  # 
  my $stmt = $dbh->prepareAndExecute($sql) || die "Can't prepareAndExecute  sql: $sql\n";

  my $sum_effective;
  my $sum_raw;
  while (my ($anatomy_id, $count) = $stmt-> fetchrow_array()) {
    $nodeHash->{$anatomy_id}->setDTRaw($count);
    $sum_raw += $count;
    $sum_effective += $nodeHash->{$anatomy_id}->getDTEffective();
  }

  return ($sum_effective, $sum_raw);
}

sub printTree {
  my ($self,$root) = @_;
  $root->printNode();
}

1;
