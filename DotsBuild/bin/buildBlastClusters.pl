#!/usr/bin/perl 

# take in files from blastMatrixNew.pl (latest on e2k2)

# Brian Brunk 1/2001

use strict;
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;

##wannt to change so that does essentailly a graph analysis (allthought without graphs!)
##first build cliques where all members of a click have the threshhold cutoffs to all the other members.
##second stage, join cliques based on number of connections between them
##this last stage (number of connections required) should be a function of the number of members ofthe cliques
#3 ie, the more members, the more connections required...certainly  more than one for large cliques
## to avoid chimeras or other anomalies

##data structures....
## 1.  hash with values being array of identifiers in the cliques, unique random keys <cliques_id>
##    NOTE: cliquess will not merge during the initial build stage...either  an Id will only be in a single
my %cliques;                   ##the cliques....

##          clique and  should be put into the largest clique with which it is consistent
## 2.  hash with identifier => <clique_id> pairs.
my %idAssign;                   ##assignment of ids to cliques

## 3.  Something identifying linkages between  cliques..unidirectionsal
##     $link{clique_id}->{clique_id}->{input_id}
my %link; 

##cache for clone_ids
my %clones;

my $debug;

my ($verbose,$percentCutoff,$lengthCutoff,$chimeraFile,$ignoreFile,$files,
    $consistentEnds, $minConnections,$logBase,$sort,$logBase2,$logBaseMax,
    $iterateCliqueSize,$iterateLogBase,$iterateLogBase2,$iterateLogBaseMax,
    $help,$iterateDescending,$minIterateLogBase,$useAllLinks,$useCloneIds,
    $reassignLinkNodes,$iterateCliqueSizeArray,$iterateLogBaseArray,
   $gusConfigFile);

&GetOptions("verbose!"=> \$verbose,
            "percentCutoff=i"=>\$percentCutoff,
            "lengthCutoff=i" => \$lengthCutoff,
            "minConnections=i" => \$minConnections, ## min number of links required to connect cliques
            "iterateLogBaseArray=s" => \$iterateLogBaseArray, ## array of logBases to iterate on...
            "iterateCliqueSizeArray=s" => \$iterateCliqueSizeArray, ## MUST be same length as $iterateLogBaseArray..
            "iterateCliqueSize=i" => \$iterateCliqueSize, ## size of cliques above which will iterate on..
            "iterateLogBase=f" => \$iterateLogBase2, ## size of cliques above which will iterate on..
            "iterateLogBaseMax=f" => \$iterateLogBaseMax, ## size of cliques above which will iterate on..
            "logBase=f" => \$logBase2,   ## used for scaling the number of links to number in cliques
            "minIterateLogBase=f" => \$minIterateLogBase,   ## used for scaling the number of links to number in cliques in iteration
            "logBaseMax=f" => \$logBaseMax,   ## used for scaling the number of links to number in cliques smooth iteration at each logBase integer between $logBase and $logBaseMax
            "useAllLinks!" => \$useAllLinks,
            "useCloneIds=i" => \$useCloneIds,
            "chimeraFile=s"=>\$chimeraFile,
            "gusConfigFile=s"=>\$gusConfigFile,
            "ignoreFile=s" => \$ignoreFile,
            "consistentEnds!" => \$consistentEnds,
            "reassignLinkNodes!" => \$reassignLinkNodes,
            "debug!" => \$debug,
            "sort!" => \$sort,
            "help!" => \$help,
            "iterateDescending!" => \$iterateDescending,
            "files=s" => \$files );

##NOTE: will always do a linking inititally with logBase2 followed by whatever the logBase that is indicated...
$logBase = 1.5;
$iterateLogBase = 1.5;

if($help || !$files || $lengthCutoff < 10 || $percentCutoff < 50){
  print &getUsage();
  exit;
}

$minIterateLogBase = 2 unless $minIterateLogBase;  ##sets this default..

$minConnections = $minConnections ? $minConnections : 1;  ##at one should act like current closure

my $db;
my $stmt;
if($useCloneIds){
  print STDERR "Establishing dbi login\n" if $verbose;
  my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

  my $db = GUS::ObjRelP::DbiDatabase->new($gusconfig->getDbiDsn(),
					  $gusconfig->getReadOnlyDatabaseLogin(),
					  $gusconfig->getReadOnlyDatabasePassword,
					  $verbose,0,1,
					  $gusconfig->getCoreSchemaName);
  my $dbh = $db->getQueryHandle();
  $stmt = $dbh->prepare("
select /* +RULE */ distinct ls.clone_id 
from dots.assemblysequence a, dots.est est 
where a.assembly_na_sequence_id = ? 
and est.na_sequence_id = a.na_sequence_id");
}

my @files = split(", *",$files);

##iterateArray...
my @itLBArr = split(", *",$iterateLogBaseArray);
my @itCSArr;
my $goodArrays = 1;
if($iterateLogBaseArray){
  if($iterateCliqueSizeArray){
    @itCSArr = split(", *",$iterateCliqueSizeArray);
    if(scalar(@itCSArr) != scalar(@itLBArr)){
      print STDERR "cliquesize and logbase arrays different lengths..using $itCSArr[0] for all clique size limits\n";
      $goodArrays = 0;
    }
    $iterateCliqueSize = $itCSArr[0];
  }elsif(!$iterateCliqueSize){
    die "You must provide --iterateCliqueSize or --iterateCliqueSizeArray in order to iterate\n";
  }
}

print STDERR "Cutoff parameters:\n\tLength: $lengthCutoff\n\tPercent Identity: $percentCutoff\n\tEnds Consistent: ".($consistentEnds ? "1" : "0")."\n\tminConnections: $minConnections\n\tlogBase: $logBase2\n\tlogBaseMax: $logBaseMax\n\titerateCliqueSize: $iterateCliqueSize\n\titerateCliqueSizeArray: '$iterateCliqueSizeArray'\n\titerateDescending: $iterateDescending\n\tminIterateLogBase: $minIterateLogBase\n\titerateLogBase: $iterateLogBase2\n\titerateLogBaseArray: '$iterateLogBaseArray'\n\titerateLogBaseMax: $iterateLogBaseMax\n\treassignLinkNodes: $reassignLinkNodes\n\n";
print STDERR "InputFiles: chimera=$chimeraFile, ignorefile=$ignoreFile\n  matrix: (",join(', ',@files),")\n";


##NOTE: for chimeras, do not want to remove if breakpoint is < $lengthCutoff  bp from either end as will not affect clustering..
## need to get length of assembly sequence from db to determine this...
my %chimera;
if ($chimeraFile) {
  ##need to bet db connections
  my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

  my $db = GUS::ObjRelP::DbiDatabase->new($gusconfig->getDbiDsn(),
					  $gusconfig->getReadOnlyDatabaseLogin(),
					  $gusconfig->getReadOnlyDatabasePassword,
					  $verbose,0,1,
					  $gusconfig->getCoreSchemaName);
  my $dbh = $db->getQueryHandle();
  my $stmt = $dbh->prepare("
select quality_end - quality_start + 1 
from dots.AssemblySequence 
where assembly_sequence_id = ?");
  print STDERR "Opening chimera file $chimeraFile\n" if $verbose;
  open(C,"$chimeraFile") || die "chimera file $chimeraFile not found\n";
	
  ##want to parse more fully to limit chimeras to >2 each side of brkpt...
  my $chim_id;
  while (<C>) {
    if (/^\>(\S+)/) {
      $chim_id = $1;
    } elsif (/Chimera\s\d+:\s1,.*numberLeft=(\d+),\snumRight=(\d+),\sbreakpoint=(\d+)/) {
      my($numLeft,$numRight,$bkpt) = ($1,$2,$3);
      if ($numLeft >= 2 && $numRight >= 2){
        #if $bkpt >= 40 check right end...if both >= 40 then is chimera..
        if($bkpt >= $lengthCutoff){
          #getLength
          $stmt->execute($chim_id);
          my $length = 0;
          while(my($len) = $stmt->fetchrow_array()){
            $length = $len;
          }
          $chimera{$chim_id} = 1 if $length - $bkpt >= $lengthCutoff;
        }
      }
    }
  }
  print STDERR "Chimeras: ",scalar(keys%chimera),"\n" if $verbose;
  close C;
  $db->logout();
}

my %ignore;
if ($ignoreFile) {
  print STDERR "Opening ignore file $ignoreFile\n" if $verbose;
  open(C,"$ignoreFile") || die "ignore file $ignoreFile not found\n";
	
  ##want to parse more fully to limit chimeras to >2 each side of brkpt...
  my $chim_id;
  while (<C>) {
    if (/^\>(\S+)/) {
      $ignore{$1} = 1;
    }
  }
  print STDERR "Ignoring: ",scalar(keys%ignore),"\n" if $verbose;
  close C;
}



##need to have two datastuctures
## 1.  %cliques: key is cluster_id and values array of included ids
## 2.  key is id and value is cluster to which assigned


my $countErr = 0;
my $clusterID = 0;
my %tmp;
my $countLines = 0;
my %graph;  ##contains the graph
my %singletons;
my %numbers;
foreach my $file (@files) {
  if ($file =~ /gz$/) {
    open(F, "gunzip -c $file |");
  } else {
    open(F,"$file");
  }
  while (<F>) {
    if (/^\>(\S+):\s\(*(.*?)\)*$/) {
      my $query = $1;
      $countLines++;
      print STDERR "Processing line $countLines\n" if($verbose && $countLines % 1000 == 0);
      next if ($chimera{$1} || $ignore{$1}); ##ignore this one
      my $noHit = 1;
      foreach my $hit (split(', ',$2)) {
        my($id,$pVal,$length,$percent,$endsCons) = split(':',$hit);
        next if ($chimera{$id} || $ignore{$id}); ##ignore this one
        if ($length >= $lengthCutoff && $percent >= $percentCutoff && (!$consistentEnds || $endsCons)){
          $noHit = 0;
          $graph{$query}->{$id} = 1;
          $graph{$id}->{$query} = 1;
        }
      }
      $singletons{$query} = 1 if $noHit;  ##is a singleton...
    } else {
      $countErr++;
      #    print STDERR "Line in incorrect format: $_";
    }
  }
  close F;
}
##need to remove singletons that have graph node...
foreach my $id (keys%singletons){
  delete $singletons{$id} if exists $graph{$id};
}

print STDERR "\nHave ",scalar(keys%graph)," graph nodes and ",scalar(keys%singletons)," singletons\n" if $verbose;
##now need to print the singletons...and clear up memory
print STDERR "printing cluster for the singletons\n" if $verbose;
foreach my $id (keys%singletons){
  $clusterID++;
  print "Cluster_$clusterID (1 sequences): ($id)\n";
  $numbers{1}++;
}
undef %singletons;  ##free this memory

print STDERR "Building cliques\n" if $verbose;

##Now build the cliques.  Go through in order and put into the largest clique.  
## if have the same number of options, check each to make certain is not  dead end 
#3 or minimally does not form a single clique when other possibilities would do better!
&makeCliques();

##now pass through cliques to see if any can be merged to form larger cliques
##note that this still requires full linkage between all ids...
&growCliques();
#&growCliques(1);
#&growCliques(2);

##try putting nodes from singleton cliques into clique with most commonnodes before linking..
##does not seem to be a good thing..need to merge smallest to largest first...
#&putSingletonCliquesIntoBestClique();

##now have the cliques and connections....connect if apropriate...
##in terms of memory, would be cheapest to just print as go along rather than save and sort
##could save %numbers so get distribution immediately as compromise...
## splitClusters.pl sorts if needed..


print STDERR "\nBuilding Clusters using links...have ",scalar(keys%cliques)," cliques: \n\n" if $verbose;

&runLinkCliques($logBase);  ##links the cliques with the current $logBase...always 1.5 to start..


if($logBaseMax){
  for(my $a = int($logBase) + 1;$a <= $logBaseMax;$a++){
    &runLinkCliques($a);  ##links the cliques with the current $logBase
  }
  if($reassignLinkNodes){ ##just run on the last one..
    &identifyAndMoveNodesCausingSingleEdges();  ##move nodes to better clique if causing a single edge to another clique
    &runLinkCliques($logBaseMax);  ##rerun linking at this logbase as moving nodes may affect it...
  }
}elsif($logBase2){
  &runLinkCliques($logBase2);  ##links the cliques with the current $logBase
  if($reassignLinkNodes){
    &identifyAndMoveNodesCausingSingleEdges();  ##move nodes to better clique if causing a single edge to another clique
    &runLinkCliques($logBase2);  ##rerun linking at this logbase as moving nodes may affect it...
  }
}


##now print the thing out!!
&printClusters($iterateCliqueSize);

&printClusterDistribution() if $verbose;

##now allow to iterate the thing...
##in this case, use the logBase2 value to iterate down from...
if($iterateDescending && $iterateCliqueSize){
  for(my $s = ($logBaseMax ? $logBaseMax : $logBase2) - 1;$s >= $minIterateLogBase + 1;$s--){
    &iterateLinkCliques($iterateCliqueSize,$logBaseMax ? undef : $s,$logBaseMax ? $s : undef);
    &printClusters($iterateCliqueSize);  ##print out..
  }
  &iterateLinkCliques($iterateCliqueSize,$logBaseMax ? undef : $minIterateLogBase,$logBaseMax ? 2 : undef);
  &printClusters();  ##print out..all remaining
}elsif($iterateLogBaseArray){
  for(my $s = 0;$s < scalar(@itLBArr);$s++){
    &iterateLinkCliques($goodArrays ? $itCSArr[$s] : $iterateCliqueSize,$itLBArr[$s],undef);
    &printClusters($s == (scalar(@itLBArr) - 1) ? undef : ($goodArrays ? $itCSArr[$s+1] : $iterateCliqueSize));  ##print out..
  }
}elsif($iterateCliqueSize){
  &iterateLinkCliques($iterateCliqueSize,$iterateLogBase2,$iterateLogBaseMax);
  &printClusters();  ##print out..all remaining
}

##print out the summary

&printClusterDistribution();

sub iterateLinkCliques {
  my($itClSz,$itLgBase2,$itLgBaseMx) = @_;
  
  ##now if $itClSz set will want to print smallest to output, remove them from graph and do larger
  ##delete nodes from graph that have been dealt with...
  &deleteGraphNodes($itClSz);
  print STDERR "\nIterating: Removing cliques smaller than $itClSz and iterating...",scalar(keys%graph)," nodes remaining\n" if $verbose;
  
  ##next ready for new data structures
  undef %cliques;
  undef %idAssign;
  undef %link;
  
  ##make cliques
  print STDERR "Iterating: Building cliques\n" if $verbose;
  &makeCliques();
  &growCliques();
  &printCliqueDistribution() if $verbose;
  
  &runLinkCliques($iterateLogBase); ##link with iterateLogBased which is always 1.5 to start
  
  if($itLgBaseMx){
    for(my $a =int($iterateLogBase)+ 1;$a <= $itLgBaseMx;$a++){
      &runLinkCliques($a);  ##links the cliques with the current $logBase
    }
    if($reassignLinkNodes){
      &identifyAndMoveNodesCausingSingleEdges();  ##move nodes to better clique if causing a single edge to another clique
      &runLinkCliques($itLgBaseMx);  ##rerun linking at this logbase as moving nodes may affect it...
    }
  }elsif($itLgBase2){
    &runLinkCliques($itLgBase2);  ##links the cliques with the current $logBase
    if($reassignLinkNodes){
      &identifyAndMoveNodesCausingSingleEdges();  ##move nodes to better clique if causing a single edge to another clique
      &runLinkCliques($itLgBase2);  ##rerun linking at this logbase as moving nodes may affect it...
    }
  }
  return 1;
}

sub printClusterDistribution {
  my $tot;
  print STDERR "Summary of distribution\n";
  foreach my $num (sort{$b <=> $a}keys%numbers) {
    print STDERR "$num\t$numbers{$num}\n";
    $tot += $num * $numbers{$num};
  }
  print STDERR "TOTAL ids: $tot\n";
}

##delete all nodes less than size
sub deleteGraphNodes {
  my($size) = @_;
  foreach my $c (keys%cliques){
    next if scalar(@{$cliques{$c}}) >= $size;  ##will have already printed
    foreach my $n (@{$cliques{$c}}){
      foreach my $ln (keys%{$graph{$n}}){  ##first delete from other end
        delete $graph{$ln}->{$n};
      }
      delete $graph{$n};
    }
  }
}

##prints all clusters smaller than $size..
sub printClusters {
  my($size) = @_;
  print STDERR "Printing clusters smaller than '$size'\n" if $verbose;
  my $totalPrinted = 0;
  if($sort){
    foreach my $c (sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})}keys%cliques){
      next if $size && scalar(@{$cliques{$c}}) >= $size;
      ##now print... and generate numbers...
      my $num = scalar(@{$cliques{$c}});
      $numbers{$num}++;
      print "Cluster_$c \($num sequences\): \(",join(', ',@{$cliques{$c}}),"\)\n"; 
      $totalPrinted += $num;
    }
  }else{
    foreach my $c (keys%cliques){
      next if $size && scalar(@{$cliques{$c}}) >= $size;
      ##now print... and generate numbers...
      my $num = scalar(@{$cliques{$c}});
      $numbers{$num}++;
      print "Cluster_$c \($num sequences\): \(",join(', ',@{$cliques{$c}}),"\)\n"; 
      $totalPrinted += $num;
    }
  }
  print STDERR "  Printed $totalPrinted identifiers at this size cutoff...\n";
}

sub printCliqueDistribution {
  print STDERR "\nClique distribrution\n";
  my %t;
  foreach my $c (keys%cliques){
    $t{scalar(@{$cliques{$c}})}++;
  }
  my $tot;
  foreach my $c (sort{$b <=> $a}keys%t){
    print STDERR "  $c: $t{$c}\n";
    $tot += ($c * $t{$c});
  }
  print STDERR "Total ids contained: $tot\n";
  
  ##check heree to see if there are duplciated identifiers in cliques..
  if($debug){
    print STDERR "Checking for duplicated ids\n";
    foreach my $c (keys%cliques){
      foreach my $i (@{$cliques{$c}}){
        unless($idAssign{$i} == $c){
          print STDERR "Clique $c: identifier is in a different clique: $idAssign{$i}\n" ;
          print STDERR "  Other click s not valid \n" unless $cliques{$idAssign{$i}}->{$i};
        }
      }
    }
  }
}

sub putSingletonCliquesIntoBestClique {
  print STDERR "\nPutting singleton cliques into the clique with most common nodes\n" if verbose;
  foreach my $c (keys%cliques){
    next if scalar(@{$cliques{$c}}) > 1;
    my $node = $cliques{$c}->[0];
    my $best = &getClusterForNodeWithMostEdges($node);
    &moveNodeToNewCluster($node,$best);
  }
  &printCliqueDistribution() if $verbose;
}

##identify problems by finding ids thatt form single edge...remove from clique and add to better one..
##where better one is the clique which contains most of the nodes it  is connected to..
sub identifyAndMoveNodesCausingSingleEdges {
  print STDERR "Moving nodes causing single edges to cliques...\n" if $verbose;
  foreach my $c (keys%link){
    foreach my $l (keys%{$link{$c}}){
      next unless &getLinkCountForClique($c,$l) <= &getLinkCountForClique($l,$c);
      my @n = &getNodesInLinkForClique($c,$l);
      foreach my $n (@n){
        my $best = &getClusterForNodeWithMostEdges($n);
        &moveNodeToNewCluster($n,$best);
      }
    }
  }
  &printCliqueDistribution() if $verbose;
}

##should look ahead one if have tie for most edges..
sub getClusterForNodeWithMostEdges {
  my($node) = @_;
  my %cl;
  foreach my $n (keys%{$graph{$node}}){
    $cl{$idAssign{$n}}++;
  }
  my @sort = sort{$cl{$b} <=> $cl{$a}}keys%cl;
  ##want to return $node unless there is one with > number of edges
  foreach my $n (@sort){
    last if $cl{$n} < $cl{$sort[0]};
    return $idAssign{$node} if $n eq $idAssign{$node};
  }
  if($verbose){
    print STDERR "getClusterForNodeWithMostEdges ($node->$idAssign{$node}) returns: ";
    foreach my $n (@sort){
      print STDERR "$n ($cl{$n}), ";
    }
    print STDERR "\n";
  }
  return $sort[0];
}

sub moveNodeToNewCluster {
  my($node,$new) = @_;
  my $old = $idAssign{$node};
  return 1 if $old == $new;

  print STDERR "Moving node $node from $old (",scalar(@{$cliques{$old}}),") to $new (",scalar(@{$cliques{$new}}),")\n" if $verbose;

  ##remove node  from old clique..
  my $oldIsValid = &removeNodeFromClique($old,$node); ##$oldIsValid = 0 means has been deleted b/c was singleton..
  print STDERR "\$oldIsValid='$oldIsValid'\n" if $debug;

  ##change $idAssign and put into new clique
  $idAssign{$node} = $new;
  push(@{$cliques{$new}},$node);

  ##now the links
  &reassignLinksForCliques($new,$old);

  ##delete the clique and links for old entirely if !$oldIsValid
  if(!$oldIsValid){ 
    &deleteEmptyClique($old);
  }

  ##take care of clone_ids..NOT simple as don't have mapping from clone_ids to nodes..what to do??
  ##best may be to just remove the old and new cluster from the clone hash so gets recomputed..
  delete $clones{$old};
  delete $clones{$new};
    
  return 1;
}

sub reassignLinksForCliques {
  my(@cs) = @_;
  foreach my $c (@cs){
    &deleteLinksForClique($c);
  }
  foreach my $c (@cs){
    foreach my $n (@{$cliques{$c}}){
      foreach my $e (keys%{$graph{$n}}){
        my $l = $idAssign{$e};
        next if $l == $c;
        $link{$c}->{$l}->{$n} = 1;
        $link{$l}->{$c}->{$n} = 1;
        $link{$c}->{$l}->{$e} = 1;
        $link{$l}->{$c}->{$e} = 1;
      }
    }
  }
}

sub deleteEmptyClique {
  my($c) = @_;
  &deleteLinksForClique($c);
  delete $cliques{$c};
}

sub deleteLinksForClique {
  my($c) = @_;
  foreach my $l (keys%{$link{$c}}){ 
    delete $link{$l}->{$c}; 
    delete $link{$l} if scalar(keys%{$link{$l}}) == 0;
  }
  delete $link{$c}; 
}

##returns number of nodes remaining in the clique..
sub removeNodeFromClique {
  my($c,$n) = @_;
  my @tmp;
  foreach my $id (@{$cliques{$c}}){
    push(@tmp,$id) unless $id eq $n; ##ids may not be ints..
  }
  if(scalar(@tmp) ==  0){
    delete $cliques{$c};
  }else{
    @{$cliques{$c}} = @tmp;
  }
  return scalar(@tmp);
}

my %haveProcessed;
sub runLinkCliques {
  my($lb) = @_;
  $logBase = $lb if $lb;
  my $countLines = 0;
  print STDERR "Linking cliques using logBase $logBase:  " if $verbose;
  my $totMerged = 1;
  my $pass = 0;
  while($totMerged){
    $pass++;
    print STDERR "runLinkCliques: Pass $pass:\n" if $verbose;
    $totMerged = 0;
    undef %haveProcessed;
    foreach my $c (sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})} keys%cliques){
      next if $haveProcessed{$c};
      next unless $cliques{$c};  ##will be deleting after dealt with...
      $totMerged += &linkCliquesNew($c);
      $countLines++;
      print STDERR "Linking cliques $countLines\n" if($debug && $countLines % 1000 == 0);
      $haveProcessed{$c} = 1;
    }
    print STDERR "Merged $totMerged cliques\n" if $verbose;
  }
  if($verbose){
    print STDERR "clique distribution after linking with logBase $logBase:\n";
    &printCliqueDistribution();
  }
}

##method to go through cliques and  grow them...ie join cliques if can make larger valid clique..
##returns the number of merge events...
sub growCliques {
  my($misMatch) = @_;
  $misMatch = $misMatch ? $misMatch : 0;
  my $countLines = 0;
  my $totMerged = 1;
  my $pass = 0;
  print STDERR "growCliques: misMatch=$misMatch\n" if $verbose;
  while($totMerged){
    $pass++;
    print STDERR "  Pass $pass:\n" if $verbose;
    $totMerged = 0;
    foreach my $c (sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})} keys%cliques){
      next unless $cliques{$c};  ##will be deleting after dealt with...
      foreach my $l (sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})}&getLinkedCliques($c)){
        if(&cliquesShouldMerge($c,$l,$misMatch)){
          $totMerged += &mergeCliques($c,$l);
        }
      }
      $countLines++;
      print STDERR "growing cliques $countLines\n" if($debug && $countLines % 1000 == 0);
    }
    print STDERR "Merged $totMerged cliques\n" if $verbose;
  }
  if($verbose){
    print STDERR "clique distribution after growing Cliques:\n";
    &printCliqueDistribution();
  }

}

##return 1 if two cliques can be merged to form larger clique
sub cliquesShouldMerge {
  my($c,$l,$misMatch) = @_;
  my $ct = 0;
  foreach my $a (@{$cliques{$c}}){
    foreach my $b (@{$cliques{$l}}){
      $ct++ if $graph{$a}->{$b};
    }
  }
  print STDERR "cliquesShouldMerge-misMatch=$misMatch: $c(",scalar(@{$cliques{$c}}),") - $l(", scalar(@{$cliques{$l}}),") have $ct edges\n" if $debug;
  return 1 if $ct >= scalar(@{$cliques{$c}}) * scalar(@{$cliques{$l}}) - $misMatch;
}

sub getLinkedCliques {
  my $c = shift;
  return keys%{$link{$c}};
  ##below because was having problems keeping links straight...for testing only!!
  my @tmp;
  foreach my $c (keys%{$link{$c}}){
    push(@tmp,$c) if $cliques{$c};
  }
  return @tmp;
}

##Now build the cliques.  Go through in order and put into the largest clique.  
## if have the same number of options, check each to make certain is not  dead end 
#3 or minimally does not form a single clique when other possibilities would do better!
sub makeCliques {
  my $countLines = 0;
#  foreach my $id (keys%graph){
  foreach my $id (sort{scalar(keys%{$graph{$a}}) <=> scalar(keys%{$graph{$b}})}keys%graph){
#  foreach my $id (sort{scalar(keys%{$graph{$b}}) <=> scalar(keys%{$graph{$a}})}keys%graph){
    print STDERR "Processing $id\n" if $debug;
    $countLines++;
    print STDERR "making cliques $countLines\n" if($verbose && $countLines % 1000 == 0);
    &processSet($id);
  }
  &printCliqueDistribution() if $verbose;
}


##don't do it iteratively and only try to merge to clicks larger than self....will iterate
##at a higher level...calling this method...
sub linkCliquesNew {
  my $c = shift;
  next unless $cliques{$c};
  ##first get th set we are dealing with...need both links and revlinks
  my $ctMerged = 0;
  my @s = sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})} &getLinkedCliques($c);
  for(my $a = 0;$a<scalar(@s);$a++){
#    next if scalar(@{$cliques{$c}}) > scalar(@{$cliques{$s[$a]}});  ##only want to merge into larger cliques than self
    next if $haveProcessed{$s[$a]};  ##think it  is valid to merge heree even if $haveProcessed as the above sort is valid
    if(&getLinkage($c,$s[$a])){  ##if link meets criteria, then merge
      my @good;
      push (@good,$s[$a]);
      my $n;
      for($n = $a+1;$n<scalar(@s);$n++){
        last if scalar(@{$cliques{$s[$a]}}) != scalar(@{$cliques{$s[$n]}}); ##only looking at ones of same size
#        print STDERR "Self: $c, First: $s[$a], current:$s[$n], array(",join(',',@s),")\n";
        push(@good,$s[$n]) if &getLinkage($c,$s[$n]);
      }
      my $best;
      if($useAllLinks){
        $best = &getBestLinkedClique_allLinks($c,\@good);
      }else{
        $best = &getBestLinkedClique($c,\@good);
      }
      &mergeCliques($c,$best);
      $ctMerged++;
      $haveProcessed{$best} = 1;
      last;  ##want to only do one merge / pass for each cluster..
      $a = $n;
    }
  }
  return $ctMerged;
}

sub getBestLinkedClique {
  my($c,$good) = @_;
  return $good->[0] if scalar(@{$good}) == 1;
  print STDERR "Finding best linked clique among ",scalar(@{$good})," cliques of size ",scalar(@{$cliques{$good->[0]}}),", (",join(', ',@{$good}),")\n" if $debug;
  my @all;
  foreach my $l (&getLinkedCliques($c)){
    push(@all,$l) if &getLinkage($c,$l);
  }
  push(@all,$c);  ##put self on list as want to take it into account as well!
  print STDERR " Checking against cliques (",join(', ',@all),")\n" if $debug;
  my %ct;
  foreach my $g (@{$good}){
    foreach my $l (@all){
      next if $l == $g || !$link{$g}->{$l} || !&getLinkage($g,$l);
      $ct{$g} += &getLinkWeight($g,$l);
    }
  }
  my @sort = sort{$ct{$b} <=> $ct{$a}}keys%ct;
  my @val = sort { $b <=> $a } values %ct;
  print STDERR "  returns first one in array (",join(', ',@val),")\n" if $debug;
  return $cliques{$sort[0]} ? $sort[0] : $good->[0];
}

##returns the relative weight of a linkage
##foreach clique the weight = number of ids in link/number of nodes in clique
##return sum of these
sub getLinkWeight {
  my($c,$l) = @_;
  return (&getCliqueLinkWt($c,$l) + &getCliqueLinkWt($l,$c));
}

sub getCliqueLinkWt {
  my($c,$l) = @_;
  return 0 unless $cliques{$c} && $link{$c}->{$l};
  return (&getLinkCountForClique($c,$l) / scalar(@{$cliques{$c}}));
}

##return the clique from the set that has the most total linkages with any clique linked to this one..
sub getBestLinkedClique_allLinks {
  my($c,$good) = @_;
  return $good->[0] if scalar(@{$good}) == 1;
  print STDERR "Finding best linked clique among ",scalar(@{$good})," cliques of size ",scalar(@{$cliques{$good->[0]}}),"\n" if $debug;
  my @all = &getLinkedCliques($c);
  push (@all,$c);
  my %ct;
  foreach my $g (@{$good}){
    foreach my $l (@all){
      next if $l == $g;
      $ct{$g} += &countTotalLinks($g,$l);
    }
  }
  my @sort = sort{$ct{$b} <=> $ct{$a}}keys%ct;
  my @val = sort { $b <=> $a } values %ct;
  print STDERR "  returns first one in array (",join(', ',@val),")\n" if $debug;
  return $cliques{$sort[0]} ? $sort[0] : $good->[0];
}

##link the cliques together..
##go through clique links and put together as meets cutoffsjj
sub linkCliques {
  my $c = shift;
  next unless $cliques{$c};
  my $haveMerge = 0;
  foreach my $l (sort{scalar(@{$cliques{$a}}) <=> scalar(@{$cliques{$b}})}&getLinkedCliques($c)){
    if(&getLinkage($c,$l)){  ##if link meets criteria, then add to list to be processed..
      $haveMerge = 1;
      &mergeCliques($c,$l);
    }
  }
  &linkCliques($c) if $haveMerge;  #if have merged cliques then potentially could now be another clique to link...
}

sub countTotalLinks {
  my($c,$l) = @_;
  return $link{$c}->{$l} ? scalar(keys %{$link{$c}->{$l}}) : 0;
}

#"minConnections=i" => \$minConnections, ## min number of links required to connect cliques
##if this is one then would be just like now...
#"logBase=f" => \$logBase,
##Note: can't require more connections than there are members of the smallest clique
sub getLinkage {
  my($c,$l) = @_;
  my $numLinks = &getMinLinkCount($c,$l);
  if(!$numLinks){
    print STDERR "ERROR: getLinkage....\$numLinks == 0 .... $c -> $l: ",scalar(keys %{$link{$c}->{$l}}),"\n" if $verbose;
    return 0;
  }

  ##make use of clone_ids if not sufficient to meet logBase and logBase >= 2
  my $ret = $numLinks >= &getLogBaseReq($c,$l) ? 1 : 0;
  if($useCloneIds && $logBase >= 2 && !$ret){
    $ret = &getCloneLinkage($c,$l);
  }
  if($logBase){ ##use this as method of choice if set...
    return $numLinks if $ret;
  }else{
    die "logBase must be set to run buildClusters.pl...\n".&getUsage();
  }
}

sub getLogBaseReq {
  my($c,$l) = @_;
  my $a = scalar(@{$cliques{$c}});
  my $b = scalar(@{$cliques{$l}});
  my $minClique = $a <= $b ? $a : $b;
  return int(log($minClique)/log($logBase)) + 1;
}

##need cache here for storing clone_ids with cliques
##%clones;
sub getCloneLinkage {
   my($c,$l) = @_;
   print STDERR "Using clone_ids to link $c and $l: ";# if $debug;
   if(!$clones{$c}){ &getCloneIds($c); }
   if(!$clones{$l}){ &getCloneIds($l); }
   my $ct = 0;
   foreach my $i (keys %{$clones{$c}}){
     $ct++ if $clones{$l}->{$i};
   }
   my $ret = $ct >= &getLogBaseReq($c,$l) ? 1 : 0;
   print STDERR "$ct linkages returns $ret\n";# if $debug;
   return $ret;
}

sub getCloneIds {
  my $c = shift;
  foreach  my $id (@{$cliques{$c}}){
    $stmt->execute($id);
    while(my ($cid) = $stmt->fetchrow_array()){
      $clones{$c}->{$cid} = 1;
    }
  }
}

##gets the set of cliques hit by the ids...
sub processSet {
  my($id) = @_;
  ##first check to see if any of the ids match existing cliques
  ##this serves two purposes, first to determine if this id can be added to any of these cliques
  #3 second, to make links to cliques that don't contain this id.
  my %list;
  foreach my $i (keys%{$graph{$id}}) {
    if ($idAssign{$i}) {
      #      push(@list,$idAssign{$i});
      $list{$idAssign{$i}}++;  ##keeps total of these ids in cliques
    }else{
      my $ncid = &makeNewClique($i);
      $list{$ncid}++;
    }
  }
#  print STDERR "Ids that match query.$id: \(", join(', ', keys%{$graph{$id}}), "\)\n" if $debug == 1;

  ##sorted list of cliques which contains ids hit by this node
  my @sort = sort{$list{$a} <=> $list{$b}}keys%list;
#  my @sort = sort{$list{$b} <=> $list{$a}}keys%list;
  if($debug){
    print STDERR "List of cliques:\n";
    foreach my $l (@sort){
      print STDERR "  $l: $list{$l}, clique=",scalar(@{$cliques{$l}}),"\n";
    }
  }

  ##Go through cliques on @$list and find largest one that contains only ids from @$ids..add $id to this clique
  #3

  my $inClique;
  my $tmpClique;
  ##first...check to see if already in clique..
  ##if it is by itself...then can be added to the next clique it hits....
  if($idAssign{$id}){
    $tmpClique = $idAssign{$id};
  }
  if($tmpClique && scalar(@{$cliques{$tmpClique}}) > 1){
    print STDERR "Already in clique with ",scalar(@{$cliques{$tmpClique}})," members\n" if $debug;
    $inClique = $tmpClique;
  }else{
    for(my $c = 0;$c < scalar(@sort);$c++){
      next if $sort[$c] == $tmpClique;  ##is self...
      if($list{$sort[$c]} == scalar(@{$cliques{$sort[$c]}})){  ##the query should  be added to this clique 
        ##now want to make certain don't choose the worst one if multiple have same number..
        ##first get all with the same number...
        my @good;
        push(@good,$sort[$c]);
        for(my $a = $c+1;$a<scalar(@sort);$a++){
          next if $sort[$a] == $tmpClique;  ##is self...
          last if ($list{$sort[$c]} != $list{$sort[$a]} || scalar(@{$cliques{$sort[$c]}}) !=  scalar(@{$cliques{$sort[$a]}}));
          push(@good,$sort[$a]) if $list{$sort[$a]} == scalar(@{$cliques{$sort[$a]}});
        }
        ##now @good has all possibilities...check that don't get worst one...
        my $best = &getBestClique($id,\@good);
        print STDERR "Adding query to clique $best with ",scalar(@{$cliques{$best}})," members\n" if $debug;
        if($tmpClique){
          &mergeCliques($best,$tmpClique);
          $inClique = $best;
        }else{
          push(@{$cliques{$best}},$id);
          $inClique = $best;  ##which clique this one is in...
          $idAssign{$id} = $best;
        }
        last;
      } 
    } 
    if(!$inClique){  ##was not placed into a clique...create new one...or use tmpClique if it exists
      $inClique = $tmpClique ? $tmpClique : &makeNewClique($id);
    }
  }
  ##then need to build connections...
  foreach my $c (@sort){
    next if $c == $inClique;  ##self....
    $link{$inClique}->{$c}->{$id} = 1;
    $link{$c}->{$inClique}->{$id} = 1;
  }
}

sub getNodesInLinkForClique {
  my($c,$l) = @_;
  my @tmp;
  foreach my $n (keys%{$link{$c}->{$l}}){
    push(@tmp,$n) if $idAssign{$n} == $c;
  }
  return @tmp;
}

sub getLinkCounts {
  my($c,$l) = @_;
  my %ct;
  foreach my $n (keys%{$link{$c}->{$l}}){
    $ct{$idAssign{$n}}++;
  }
  return ($ct{$c},$ct{$l});
}

##note will return the number of links to the first clique
sub getLinkCountForClique {
  my($c,$l) = @_;
  my($a,$b) = &getLinkCounts($c,$l);
  return $a;
}

sub getMinLinkCount {
  my($c,$l) = @_;
  my($a,$b) = &getLinkCounts($c,$l);
  return $a < $b ? $a : $b;
}

sub getBestClique {
  my($id,$good) = @_;
  return $good->[0] if scalar(@{$good}) == 1;
  my @allIds = keys%{$graph{$id}};
  push(@allIds,$id);  ##add self so all come through the following..
  print STDERR "getBestClique: ",scalar(@$good)," cliques, ",scalar(@allIds)," ids\n" if $debug;
  my %num;
  foreach my $c (@{$good}){
    foreach my $i (@{$cliques{$c}}){
      foreach my $h (@allIds){
        $num{$c}++ if $graph{$i}->{$h};
      }
    }
  }
  my @sort = sort{$num{$b} <=> $num{$a}}keys%num;
  print STDERR "getting Best Clique\n" if $debug;
  foreach my $a (@sort){
    print STDERR " $a: $num{$a}\n" if $debug;
  }
  return $sort[0];
}

##add ids from $merge to $keep
##change pointers for  those ids to point to $keep
##change link to reflect the new situation
sub mergeCliques {
  my($keep,$merge) = @_;
  print STDERR "Merging $merge (",scalar(@{$cliques{$merge}})," nodes) into $keep (",scalar(@{$cliques{$keep}})," nodes)\n" if $debug;
  foreach my $l (keys%{$link{$merge}}){  
    ##transfer to new clique
    if($l != $keep){ ##don't link to self!!
      foreach my $i (keys%{$link{$merge}->{$l}}){
        print STDERR "transferring link $i from $merge to $keep\n" if $debug;
        $link{$keep}->{$l}->{$i} = 1;
        $link{$l}->{$keep}->{$i} = 1;
      }
    }
    delete $link{$l}->{$merge};
  }
  delete $link{$merge};

  ##change $idAssign and put into new clique
  foreach my $n (@{$cliques{$merge}}){
    print STDERR "ERROR: \$merge=$merge, \$idAssign.$n = $idAssign{$n}\n" unless $idAssign{$n} == $merge;
    print STDERR "assigning $n from clique $idAssign{$n} to $keep\n"  if $debug;
    $idAssign{$n} = $keep;
    push(@{$cliques{$keep}},$n);
  }

  ##change any  clone_id links if exist
  foreach my $i (keys%{$clones{$merge}}){
    $clones{$keep}->{$i} = 1;
  }
  delete $clones{$merge};

  ##delete from cliques..
  delete $cliques{$merge};
  return 1;
}

sub makeNewClique {
  my($id) = @_;
  print STDERR "makeNewClique..ERROR: $id is already in clique.$idAssign{$id}\n" if $idAssign{$id};
  $clusterID++;
  push(@{$cliques{$clusterID}},$id);
  $idAssign{$id} = $clusterID;
  return $clusterID;
}
  
sub getUsage {
print <<usage;
USAGE:  buildBlastClusters.pl 
  --lengthCutoff=i    ##length cutoff for similarities to create an edge in graph
  --percentCutoff=i   ##percent cutoff for similarities to create an edge in graph
  --consistentEnds!   ##ends must be consistent if true to create an edge in graph
  --chimeraFile=<filename>   ##file containing chimeras identified by blast matrix...will be removed from graph
  --ignoreFile <filename of sequences to ignore>  ##file of identifiers to NOT put into graph
  --verbose!   ##verbose output
  --files '<matrix, files>'  ##matrix files from blast matrix (can be gzipped or not)
  --debug!  ##lots of debugging output
  --logBase=f  ##logBase to use for determining linkages (NOTE: logBase 1.5 is always used initially to merge small
                 clusters prior to using this logBase value in second pass
  --logBaseMax=f  ##like logBase but increments up by ints to the logBaseMax from logBase 2 (ie, if 4 then does 3 then 4).
  --iterateCliqueSize=i  ##cutoff size of cliques to accept in first iteration. Nodes from these cliques will be removed
                           and the merging re-run using iterateLogBase as logBase value.
  --iterateLogBase=f  ## value for logBase to use for large clusters...should be less that logBase!
  --iterateLogBaseArray=s" =>  ## array of logBases to iterate on...
  --iterateCliqueSizeArray=s" =>  ## MUST be same length as $iterateLogBaseArray..
  --iterateLogBaseMax=f  ##like logBaseMax but for iteration
  --iterateDescending!   ##if true then iterates the iteration stating with logBase - 1 down to --minIterateLogBase
  --minIterateLogBase=f [2] ##if --iterateDescending, determines minimum logBase to use
  --reassignLinkNodes!   ##if true then detect nodes involved in linkages less then logBase and reassigns them to the cluster with the most common edges
  --useAllLinks!      ##use total links between cliques for determining best clique to merge into
  --useCloneIds=i      ##uses clone_ids when clusters don't meet logBase but have links
  --sort!  ##sort output by cluster size (note if iterating then will get two sorts..first smaller than iterateCliqueSize
             and then the iterated cliques.
  --help  ##prints this usage statement

usage
}

