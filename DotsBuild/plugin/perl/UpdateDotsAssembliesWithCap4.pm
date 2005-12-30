package DoTS::DotsBuild::Plugin::UpdateDotsAssembliesWithCap4;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use GUS::PluginMgr::Plugin;
use GUS::Model::DoTS::Assembly;
use GUS::Model::DoTS::AssemblySequence;
use GUS::Model::DoTS::MergeSplit;

my $debug = 0;

my $argsDeclaration =
[
    stringArg({
        name => 'clusterfile',
        descr => 'name of cluster file for input',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'cap4Dir',
        descr => 'location of executable cap4',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'directory',
        descr => 'location of working directory accessible from both current machine and cap4_machine',
        default => `pwd`,
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'remote_dir',
        descr => 'working directory on cap4_machine',
        default => "/scratch1/$ENV{USER}/cap4",
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'debug_assem_file',
        descr => 'cap4 output file to  be parse in and iterated on for debugging Assembly',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    booleanArg({
        name => 'assemble_old',
        descr => 'does not assemble clusters with only old ids unless true',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    booleanArg({
        name => 'no_delete',
        descr => 'if true, tags Assembly.description with \'DELETED\' rather than deleting',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    booleanArg({
        name => 'debugPlugin',
        descr => 'if true, turns debugging on specifically in plugin...not relationalrow',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    integerArg({
        name => 'max_iterations',
        descr => 'maximum number of times to iterate',
        default => 1,
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'cap4_machine',
        descr => 'node to rsh to to run cap4',
        default => 'server',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    stringArg({
        name => 'cap4_params',
        descr => 'parameters to pass to cap4',
        default => 'Verbosity=0 -NoRecover -NoPolyBaseMask MinCovRep=500 InOverhang=30 EndOverhang=30 RemOverhang=30 QualSumLim=300 MaxInternalGaps=15 -ESTAssembly MaxOverlaps=50 -KeepDups',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    integerArg({
        name => 'testnumber',
        descr => 'number of iterations for testing',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    booleanArg({
        name => 'reassemble',
        descr => 'Reassembles entirely from AssemblySequences rather than incremental',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    }),
    integerArg({
        name => 'taxon_id',
        descr => 'taxon_id for these assemblies (8=human,14=mouse)',
        constraintFunc => undef,
        reqd => 0,
        isList => 0
    })
];

my $purposeBrief = <<PURPOSEBRIEF;
Incremental update of DOTS Assemblies
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
Incremental update of DOTS Assemblies: reassembles entirely from AssemblySequences if --reassemble
PLUGIN_PURPOSE

#check the documentation for this
my $tablesAffected = [
    ['DoTS::AssemblySequence', '']
];

my $tablesDependedOn = [
    ['DoTS::AssemblySequence', ''],
    ['DoTS::Assembly', ''],
    ['DoTS::MergeSplit', '']
];

my $howToRestart = <<PLUGIN_RESTART;
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
		     cvsRevision => '$Revision$', # cvs fills this in!
		     name => ref($self),
		     argsDeclaration   => $argsDeclaration,
		     documentation     => $documentation
		    });
    
    return $self;
}


$| = 1;

##Global variables##
my $ctx;
my @oldIds;
my @assIds;
my %mapAss;
my %mapAssSeq; 
my $iterateParams;
my $iterationNumber = 0;
my @trimmedAssemblySequences;   ##global list of assembly sequences that get removed entirely by trimming

my $assCache;                   ##note that this will be used to manage the assembly cache....nothing more!!
my $cap4;
my $tmpLib = "tmpLib";
my $count = 0;
my $oldTotal = 0;
my $newTotal = 0;
my $msGroupId = 0;              ##global variable so that all merges and splits from a single "Cluster" are grouped together.
my $algoInvo;                   ##global AlgorithmInvocation used for submits...

sub run {
  my $M   = shift;
  $ctx = shift;
  print $ctx->{cla}->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $ctx->{cla}->{'testnumber'}\n" if $ctx->{cla}->{'testnumber'};
  $cap4 = $ctx->{cla}->{'cap4Dir'}."/cap4";
  if (!(-e "$cap4")) {
    die "$cap4 does not exist";
  }

  ##set no version on if not committing
  $ctx->{self_inv}->setGlobalNoVersion(1) unless $ctx->{cla}->{commit};

  ##NOTE: WANT TO NOT DELETE EVIDENCE OR  SIMILARITY BY DEFAULT TO INCREASE THROUGHPUT...DELETE THESE SEPARATELY
  $ctx->{self_inv}->setGlobalDeleteEvidenceOnDelete(0);
  $ctx->{self_inv}->setGlobalDeleteSimilarityOnDelete(0);

  ##may need more objects as do all in one transaction...
  $ctx->{'self_inv'}->setMaximumNumberOfObjects(300000);
  $algoInvo = $ctx->{self_inv}; 

  chomp $ctx->{cla}->{directory};

  if ((!$ctx->{cla}->{'clusterfile'} && !$ctx->{cla}->{debug_assem_file}) || !$ctx->{cla}->{'taxon_id'} || !$ctx->{cla}->{cap4_machine}) {
    die "You must include --clusterfile --taxon_id and --cap4_machine on the command line\n";
  }

  ##create the Asembly object for the cache..
  $assCache = GUS::Model::DoTS::Assembly->new();
  $assCache->setSetDefaultsOnSubmit(1); ##set the defaults on assemblies that are submitted

  if ($ctx->{cla}->{debugPlugin}) { ##turns on debugging if passed in on cmdline
    $debug = 1;
    $assCache->setAssemblyDebugging(1);
  }

  ##create the remote_dir on the cap4_machine
  if ($ctx->{cla}->{cap4_machine} !~ /^s/i) {
    my @tmp = split(/\//,$ctx->{cla}->{remote_dir});
    my $totPath = "";
    for (my $i = 0;$i< scalar(@tmp);$i++) {
      next unless $tmp[$i];
      system("rsh -n $ctx->{cla}->{cap4_machine} 'mkdir $totPath/$tmp[$i]'") unless $tmp[$i] =~ /scratch/;
      $totPath .= "/$tmp[$i]";
    }
  }

  ##for restarting....
  my %finished;
  if (-e "updateDOTSAssemblies.log") {
    open(LOG, "updateDOTSAssemblies.log"); ##for logging progress...errors go to stderr
    while (<LOG>) {
      ##parse out cluster ids and put into a hash so can ignore...
      if (/^(Cluster_\d+)/) {
        $finished{$1} = 1;
      }
    }
    close LOG;
    print STDERR "Restarting: have finished ".scalar(keys%finished)." clusters\n";
  }

  open(LOG, ">>updateDOTSAssemblies.log");
  select LOG; $| = 1; select STDOUT;
  #  print STDERR "<html>\n<preg>\n\n";

  ##do  the iterate params..and print the params to the log...
  $iterateParams = $ctx->{cla}->{cap4_params};
  $iterateParams =~ s/ MaxOverlaps=\d+//;
  $iterateParams =~ s/ -ClipByBadEnd//;
  ##also don't want to look for chimeras here...as depth too low?
  ##note: not sure I should look for chimeras at all!!
  #  $iterateParams =~ s/-ChimeraOut //;
  print STDERR "IterateParams: $iterateParams\n" if $debug;

  print STDERR "\nCAP4 params: $ctx->{cla}->{cap4_params}\n\n";
  print STDERR  "CAP4 iterate params: $iterateParams\n\n";

  if ($ctx->{cla}->{debug_assem_file}) {
    &runDebugFromFile();
  } else {
    open(F,"$ctx->{cla}->{'clusterfile'}") || die "clusterfile $ctx->{cla}->{'clusterfile'} not found\n";

    while (<F>) {
      ##reset variables...
      undef @oldIds;
      undef @assIds;
      undef %mapAss;
      undef %mapAssSeq;
      undef @trimmedAssemblySequences;
      $iterationNumber = 0;
      $assCache->undefPointerCache(); ##removes circular references for parent/children so can be garbage collected...
      $assCache->flushAssemblyCache(); ##flush the assembly cache...
      $assCache->flushAssemblySequenceCache(); ##flush the assemblySequence cache...
      $algoInvo->removeAllChildren(); ##removes all children
      $msGroupId++;
      my $totalNewIds = 0;
      
      #	print STDERR "Line:$_" if $debug;
      
      ##mechanism to stop clustering and restart so can free ram after very large clusters..
      if (/^EXIT/) {
        print LOG "Exiting to free ram..will restart\n";
        die "Exiting to free RAM...will  restart\n";
      }
      
      if (/^(Cluster_\S+)\s.*:\s\((.*)\)/) {
        my ($countContigs,$countSinglets,$singlets,$align) = (0,0,undef,undef);
        my $cluster = $1;
        next if exists $finished{$cluster}; # && !$debug;
        $count++;
        last if $ctx->{cla}->{'testnumber'} && $count > $ctx->{cla}->{'testnumber'}; 
        foreach my $id (split(', ',$2)) {
          if ($id =~ /^DT\.(\d+)/) {
            push(@oldIds,$1);
          } else {
            push(@assIds,$id);
          }
        }
        
        ##add to the totals...
        $oldTotal += scalar(@oldIds);
        $newTotal += scalar(@assIds);
        
        print STDERR "processing $cluster: oldIds\(".join(', ',@oldIds)."\), newIds\(".join(', ',@assIds)."\)\n" if $debug;
        
        ##first get the sequence entries..returns a new AssemblyCluster object if there are no old ids..
        
        ##if am reassembling....
        my $okToReass = 0;
        if ($ctx->{cla}->{'reassemble'}) {
          $okToReass = &getGusEntriesForReassembly();
        } else {
          $okToReass = &getGusEntries(); ##retrieve all the existing dots_ids and merge genes if necessary
        }
        
        if (!$okToReass) {      ##returns undef if one (or more) of the DOTS entries are not found!!
          print LOG "$cluster ERROR: DOTS entry for $cluster not found\n";
          next;
        }

        print STDERR "Starting $cluster: ",`date` if $debug;

        ##next get the new sequences...cache and get/generate AssemblySequences...
        &getNewSequences();     ##get the new sequences

        $totalNewIds = scalar(@assIds);

        ##NOTE: it is possible that there are no sequences if the assemblysequences couldnot be retrieve
        ##from the db..
        if (scalar(@oldIds) == 0 && $totalNewIds == 0) { ##there are no sequences....don't do anything
          print LOG "$cluster finished: contains NO sequences\n";
          next;
        }

        ##if is singleton then don't need to run cap2...
        ##create new Assembly..add to gene and am finished...
        if (scalar(@oldIds) == 0 && $totalNewIds == 1) { ##it is a singleton....
          print STDERR "Prcessing new singleton\n" if $debug;
          ##first get the cached AssemblySequence
          my $singAssSeq = $assCache->getCachedAssemblySequence($assIds[0]);
          
          ##next create a new Assembly and add it to algoInvo for submitting..
          $algoInvo->addChild(&makeNewAssembly($singAssSeq));
          
          ##submit
          my ($subRet,$reviewed,$delReviewed) = &submitUpdatedAssemblies();
          my($assCt,$assSeqCt) = &countAssembliesAndAssemblySequences();
          if ($subRet) {
            print LOG "$cluster finished: $assCt Assemblies: ",scalar(@oldIds)," DT, $assSeqCt total\n";
            if ($reviewed || $delReviewed) {
              foreach my $r (@$reviewed) {
                print LOG "  REVIEWED: rna.",$r->[0],", DT.",$r->[1],"\n";
              }
              foreach my $r (@$delReviewed) {
                print LOG "  REVIEWED merged: rna.",$r->[0],", DT.",$r->[1],"\n";
              }
            }
          } else {
            print LOG "$cluster ERROR: $cluster NOT SUBMITTED: $assCt Assemblies from $assSeqCt sequences\n"; 
          }
          next;
        }
        
        ##if contains only old DOTS_ids go to next unless $ctx->{cla}->{assemble_old}
        if ($totalNewIds == 0) { ##there are no new sequences...only old (happens in transClosure when bring in cluster info)
          if ( !$ctx->{cla}->{assemble_old} ) { ##is singleton or not assembling only old ones...don't bother to assemble...
            if ($ctx->{cla}->{'reassemble'}) {
              print LOG "$cluster ERROR: unable to retrieve any AssemblySequences\n" ;
              ##NOTE: need to submit here as need to update the assemblysequences that may be low quality...
              my ($subRet,$reviewed,$delReviewed) = &submitUpdatedAssemblies();
              next;
            }
            print LOG "$cluster finished: ",scalar(@oldIds), " DOTS assemblies with no new sequences\n";
            next;
          }
          if (scalar(@oldIds) == 1) {
            print LOG "$cluster finished: 1 DOTS assembly with no new sequences\n";
            next;
          }
          ##NOTE..want to use the iterate parameters for assembling as all are contigs
          ## and thus want to limit trimming...won't  iterate  at end...
          &iterateAssembly(scalar(@oldIds),0);

        } else {                ##else $totalNewIds > 0


          ##if reassembling and one oldId and one NewId then is singleton...
          if ($ctx->{cla}->{reassemble} && scalar(@oldIds) == 1 && $totalNewIds == 1) {
            print STDERR "$cluster Processing Old singleton assembly\n" if $debug;
            my $oldSingAss = $algoInvo->getFromDbCache('GUS::Model::DoTS::Assembly',$oldIds[0]);
            ##gappedSequence for the singleton  should be set..
            ##need to make certain that the single AssemblySequence is child of this sequence..
            ##if not then add it....
            if (!$oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')) {
              $oldSingAss->addChild($assCache->getCachedAssemblySequence($assIds[0]));
            }
            $oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->setGappedSequence($oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->getSequence()) unless $oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->getGappedSequence() eq $oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->getSequence();
            $oldSingAss->setGappedConsensus($oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->getSequence()) unless  $oldSingAss->getGappedConsensus() eq $oldSingAss->getChild('GUS::Model::DoTS::AssemblySequence')->getSequence();
            $assCache->cacheAssembly($oldSingAss);
            $algoInvo->addChild($oldSingAss);
           
          } else {

            ##now do the cap4....incrementally if more than 500 new sequences ??
            ##think about incrementally...cap4 seems able to do large numbers fairly linearly in time

            ##needs to be done incrementally if > 500 ass seqs
            ##assemble just the ass seqs first then use the iterate params to do with
            ##any cached assemblies...essentialy just iterate...pehaps should iterate twice??
            ##NOTE: in any case, want to only use -ClipByBadEnd if there are no consensus
            ##sequences as clips them too much...thus should incrementally assemble the
            ##new seuqences and then do an iteration with the cached (and incoming) assemblies

            print STDERR "Processing $cluster ($totalNewIds new sequences)\n" if $debug;

            for (my $i = 0;$i < $totalNewIds;$i += 500) {
              open(T,">$tmpLib");
              print T "<CAML>\n\n";
              ##first print the cached assemblies...but only if singletons...otherwise wait  for iteration
              ##this will be trying singletons with each 500...if  assemble will no longer be singletons
              foreach my $ass ($assCache->getAllCachedAssemblies()) {
                if (!$ass->isMarkedDeleted() && $ass->getNumberOfContainedSequences() == 1) { ##assemblies that are merged get marked deleted...
                  print T $ass->toCAML(1)."\n";
                  print STDERR "printing sequence for D".$ass->getCacheId()."\n" if $debug;
                }
              }

              ##don't want to assemble fewer than 100 if iterating...
              my $amt = ($totalNewIds - $i > 499 && $totalNewIds - $i < 599) ? 599 : 499;
              foreach my $assid (@assIds[$i..$i+$amt]) {
                last unless $assid;
                print STDERR "printing sequence for $assid\n" if $debug;
                my $a = $assCache->getCachedAssemblySequence($assid);
                print T $a->toCAML(1)."\n"; 
              }
              print T "</CAML>\n";
              close T;
              $i += 100 if $amt == 599;
            
              ##now run cap4 and parse output...
              my($t_countContigs,$t_countSinglets,$t_singlets,$t_align) = &processCap4();
              if (!defined $t_align && !defined $t_singlets) { ##didn't get any valid cap4 output
                next;
              }
              &processAlignment($t_align);

              ##need to now total  the things above...
              $countContigs += $t_countContigs;
              $countSinglets += $t_countSinglets;
          
            }

            if ($countSinglets == 0 && $countContigs == 0) {
              print STDERR "ERROR: cap4 did not produce any valid output..\n";
              print LOG "$cluster ERROR: cap4 did not produce any valid output, skipping\n";
              next;
            }

            ##It's party time....iterate with cap4 until freezes down  to stable set of assemblies
            ##  need to determine the rules for this very carefully!!
            ##no! just iterate once but with the MaxOverlaps not set so does full alignment
          
            &iterateAssembly($countContigs,$countSinglets);

            ##perhaps should assemble trimmed AssSeqs here rather than in iterateAssembly sub


            ##Lastly, need to assign identifiers if am reassembling
            if ($ctx->{cla}->{'reassemble'}) {
              if (! &assignIdsForReassembly()) {
                ##reassignment failed....print error msg to log and do not submit...
                print LOG "$cluster ERROR, NOT SUBMITTED: reassigning ids failed\n";
                next;
              }
            } 
          } 
        } 

        
        ##now do the submits.... 
        print STDERR "Submitting non-singleton\n" if $debug;
        my ($subRet,$reviewed,$delReviewed) = &submitUpdatedAssemblies();
        my($assCt,$assSeqCt) = &countAssembliesAndAssemblySequences();
        ##debugging information...
        if (!$ctx->{cla}->{commit} || $debug) {
          &printDebuggingInfo($debug);
        } 
        if ($subRet) {
          
          #          print G "$cluster finished: $assCt Assemblies from $assSeqCt sequences\n";
          print LOG "$cluster finished: $assCt Assemblies: ",scalar(@oldIds)," DT, $assSeqCt total\n";
          if ($reviewed || $delReviewed) {
            foreach my $r (@$reviewed) {
              print LOG "  REVIEWED: rna.",$r->[0],", DT.",$r->[1],"\n";
            }
            foreach my $r (@$delReviewed) {
              print LOG "  REVIEWED merged: rna.",$r->[0],", DT.",$r->[1],"\n";
            }
          }
          if ($ctx->{cla}->{'reassemble'} && $totalNewIds != $assSeqCt) {
            print LOG "  Number of AssemblySequences ($assSeqCt) does not match total input for reassembly ($totalNewIds)\n";
          }
        } else {
          print LOG "$cluster ERROR, NOT SUBMITTED: $assCt Assemblies from $assSeqCt sequences\n"; 
        }
        
      }
    }
  }

  #  unlink "$tmpLib";

  ############################################################
  ###  put an informative summary in the results variable
  ############################################################
  my $results = "Processed $count clusters: $oldTotal DT.na_sequence_ids and $newTotal new sequences. Inserted ".($ctx->{self_inv}->getTotalInserts() - 1).", Updated ".$ctx->{self_inv}->getTotalUpdates()." and Deleted ".$ctx->{self_inv}->getTotalDeletes();
  ###updates the AlgorithnInvocation table with results and time complete

  print LOG "\n$results\n";

  #  print STDERR "\n</pre>\n</html>\n";

  return $results;
}

##################################################
# Subroutines...
##################################################

sub printDebuggingInfo {
  my($level) = shift;
  print STDERR "\nDebugging information...\n\n" if $level;
  foreach my $ass ($algoInvo->getChildren('GUS::Model::DoTS::Assembly')) {
    print STDERR "\nFinished Assembly: DT.",$ass->getId(),", cacheId=",$ass->getCacheId(),"\nAlignment:\n" if $level; 
    print STDERR $ass->getCap2Alignment() if $level;
    print STDERR $ass->toXML(0,1) if $level;
    #              print STDERR $ass->toCAML(1);
    #              $ass->findSNPs(4,0.1);
    foreach my $a ($ass->getChildren('GUS::Model::DoTS::AssemblySequence')) {
      my $strim = $a->getSequenceStart() - $a->getQualityStart();
      my $etrim = $a->getQualityEnd() - $a->getSequenceEnd();
      print STDERR "  AssemblySequence.",$a->getParent('GUS::Model::DoTS::ExternalNASequence')->getSourceId(),": length=",($a->getSequenceEnd() - $a->getSequenceStart() + 1),", trimmed, start=$strim, end=$etrim\n";
      #                if($strim + $etrim > 60){
      #                  print STDERR "    trimmed: ",$a->getSequence(),"\n";
      #                  $a->resetAssemblySequence();
      #                  $a->setSequenceStart($a->getQualityStart());
      #                  $a->setSequenceEnd($a->getQualityEnd());
      #                  print STDERR "Not trimmed: ",$a->getSequence(),"\n";
      #                }
    }
  } 
}

##purpose here is to iterate until somehow freezes and get no changes...just how to do this is the question!!
sub iterateAssembly {
  my($countContigs,$countSinglets) = @_;
  print STDERR "Iterating assembly($countContigs contigs, $countSinglets singlets)\n" if $debug;
  $iterationNumber++;
  if ($iterationNumber > $ctx->{cla}->{max_iterations} || ($countContigs == 1 && $countSinglets == 0)) {
    ##need to assemble the trimmed assemblysequences here if there are any... 
    &processTrimmedAssemblySequences(); 
    &makeSingletonAssemblies($countSinglets);
    return $countSinglets; 
  }

  ##first print sequences...
  open(T,">$tmpLib");
  print T "<CAML>\n";
  ##first print the cached assemblies...unless reassembling...
  foreach my $ass ($algoInvo->getChildren('GUS::Model::DoTS::Assembly')) {
    if (!$ass->isMarkedDeleted()) { ##assemblies that are merged get marked deleted...
      print T $ass->toCAML(1);
      print STDERR "printing sequence for D".$ass->getCacheId()."\n" if $debug;
    }
  }
  ##next print the new ids 
  if ($countSinglets > 0) {     ##have some singlets
    foreach my $a ($assCache->getAllCachedAssemblySequences()) {
      #    foreach my $assid (@assIds){
      #      next unless $assid;
      #      my $a = $assCache->getCachedAssemblySequence($assid);
      next unless $a->isSinglet();
      #      print STDERR "printing sequence for singleton $assid\n";# if $debug;
      print T $a->toCAML(1); 
      $a->setSinglet(0);
    }
  }
  print T "</CAML>\n";
  close T;
  #  exit;
  ##then run cap4..
  my ($newCountContigs,$newCountSinglets,$newSinglets,$newAlign) = &processCap4(1);
  ##then parse and check to see if have fewer assemblies + singletons...
  ##if fewer,process the alignments else submit what I have...
  if ($newCountContigs && ($newCountSinglets < $countSinglets || $newCountContigs < $countContigs)) {
    ##need to process and then redo iteration
    &processAlignment($newAlign);
    &iterateAssembly($newCountContigs,$newCountSinglets);
  } else {
    &processTrimmedAssemblySequences(); 
    &makeSingletonAssemblies($newCountSinglets);
  }
  return $newCountSinglets;     ##what to return!!..number of singlets...
}

sub makeSingletonAssemblies {
  my($numSinglets) = @_;
  
  if ($numSinglets > 0) {
    foreach my $a ($assCache->getAllCachedAssemblySequences()) {
      next unless $a->isSinglet();
      $algoInvo->addChild(&makeNewAssembly($a));
      $a->setSinglet(0);
    }
  }
}

sub makeNewAssembly {
  my(@aseq) = @_;
  my $singAss = GUS::Model::DoTS::Assembly->new({'sequence_type_id' => 5, ##all Assemblies are of type mRNA
                               'taxon_id' => $ctx->{cla}->{'taxon_id'},
                               'subclass_view' => 'Assembly',
                               'description' => 'Not Annotated'} );
  
  print STDERR "Adding AssemblySequence child\n" if $debug == 1;
  $singAss->addChildren(@aseq);

  $singAss->getGappedConsensus();
  $singAss->getGappedLength();
  $singAss->setNumberOfContainedSequences(scalar(@aseq));
  
  $singAss->createFeaturesAndAASequence();

  $assCache->cacheAssembly($singAss);

  ##want to set the sequence gaps for the assemblysequences...should  only be 1..
  foreach my $c (@aseq) {
    $c->setGappedSequence($c->getSequence()) if scalar(@aseq) == 1;
  }
  
  #  print "Checking singleTonAssembly\n",$singAss->toXML(0,1) if $debug;
  return $singAss;
}

##need submit method here because want to deal appropriately with assemblies...
##am I submitting assemblies twice??...
##note  that need to submit assemblies that have been deleted last!!
sub submitUpdatedAssemblies {
  my %submitted;
  my @reviewed;
  my @delReviewed;              ##THose that have been deleted...but man rev

  ##need to check to see that all assemblies that are not marked deleted have assemblysequences..
  #3this  can happen if an assembly gets trimmed entirely...it's
  #  foreach my $id (@oldIds){
  #    my $ass = $algoInvo->getFromDbCache('GUS::Model::DoTS::Assembly',$id);
  #    if(!$ass->isMarkedDeleted() && scalar($ass->getChildren('GUS::Model::DoTS::AssemblySequence')) == 0){
  #      &markAssemblyAndChildrenDeleted($ass);
  #      $algoInvo->addChild($ass) unless $ass->getParent('GUS::Model::DoTS::AlgorithmInvocation');
  #    }elsif($ass->isMarkedDeleted() && $ass->getChildren('GUS::Model::DoTS::AssemblySequence')){
  #      print STDERR "ERROR: $id is marked deleted but still has AssemblySequence children\n";
  #      &printDebuggingInfo(1);
  #    }
  #  }

  $algoInvo->manageTransaction(undef,'begin');

  ##need to first submit any assemblysequence children (these have not been assembled and may have foreign key
  ## refereences to an assembly that is getting deleted..
  foreach my $a ($algoInvo->getChildren('GUS::Model::DoTS::AssemblySequence')) {
    $a->submit(1,1);            ##submit not_deep and no_trans
  }

  foreach my $c (sort{$a->isMarkedDeleted() <=> $b->isMarkedDeleted()}$algoInvo->getAllChildren(0,1)) {
    next if $c->getClassName() eq 'GUS::Model::DoTS::AssemblySequence'; ##already submitted above
    ##first don't submit if is is asembly and marked deleted and does not have an id...
    if ($c->getClassName() eq 'GUS::Model::DoTS::Assembly') {
      next if exists $submitted{$c->getCacheId()};
      $submitted{$c->getCacheId()} = 1;
      next if ( !$c->getId() && ($c->isMarkedDeleted() || scalar($c->getChildren('GUS::Model::DoTS::AssemblySequence')) < 1 ));
      #      $c->markDeleted() if scalar($c->getChildren('GUS::Model::DoTS::AssemblySequence')) < 1;
      if ($c->getId()) {        ##want to get the manually reviewed things...
        my $rna = $c->getRNA(1,1);
        if ($rna && $rna->getReviewStatusId() == 1 ) {
          if ($c->isMarkedDeleted()) {
            push(@delReviewed,[$rna->getId(),$c->getId()]);
          } else {
            push(@reviewed,[$rna->getId(),$c->getId()]);
          }
        }
      }
      ##note: if not deleting and  is assembly that is marked deleted need to catch and unmark deleted
      if ($c->isMarkedDeleted() && $ctx->{cla}->{no_delete}) {
        $c->markUnDeleted();
      }
      ###check to see if still has assemblyseuqencechildren and deal with if  so...
      ##is hack but can fix before dumping sequences so leave in until find problem
      if ($c->isMarkedDeleted()) {
        $c->retrieveChildrenFromDB('GUS::Model::DoTS::AssemblySequence');
        if (scalar($c->getChildren('GUS::Model::DoTS::AssemblySequence')) > 0) {
          $c->markUnDeleted();
          $c->removeAllChildren();
          $c->undefSubmitList();
          $c->setLength(0);
          $c->setGappedConsensus('');
          $c->setSequence('');
          $c->setDescription('ERROR: Needs to be reassembled');
          print LOG "ERROR: Assembly.",$c->getId()," still has AssembySequence children\n";
        }
      }
    }                           # else{ next; } $c->submit(1,1); ##uncomment here and comment next also return 1 in Assembly-submit.
    ## to to assembly so can output the alignments for testing without submitting.
    ##remove algoinvo parent so id does not get set if no other changed atts..
    $c->removeParent($algoInvo);
    $c->submit(undef,1);
    $c->setParent($algoInvo);
  }
  ##sanity check...some assemblies have no sequence...check here and print to STDERR if not good then roll back so
  ##can test....later
  my $ctZero = 0;
  foreach my $a ($algoInvo->getChildren('GUS::Model::DoTS::Assembly')) {
    if ($a->getLength() == 0) {
      print LOG "ERROR: Assembly.",$a->getId()," sequence length = 0\n";
      $ctZero++;
      ##Do something here!!
    }
  }
  ##  $algoInvo->setRollBack(1) if $ctZero;  ##2will roll back when manage transaction;
  return ($algoInvo->manageTransaction(undef,'commit'),\@reviewed,\@delReviewed);
}

sub countAssembliesAndAssemblySequences {
  my @tmp = $algoInvo->getChildren('GUS::Model::DoTS::Assembly');
  my $ct = 0;
  #  my @reviewed;
  #  my @delReviewed;  ##THose that have been deleted...but man rev
  my $ctAssem = 0;
  foreach my $a (@tmp) {
    $ct += $a->getNumberOfContainedSequences();
    $ctAssem++ if $a->getNumberOfContainedSequences();
    #    my $r = $a->getRNA(1);
    #    if($r && $r->getManuallyReviewed()){
    #      push(@reviewed,[$r->getId(),$a->getId()]);
    #      print STDERR "ManuallyReviewed: rna.",$r->getId(),", DT.",$a->getId(),"\n" if $debug;
    #    }
  }

  #  foreach my $rna ($algoInvo->getChildren('GUS::Model::DoTS::RNA')){
  #    if($rna->getManuallyReviewed()){
  #      push(@delReviewed,[$rna->getId(),$rna->getNASequence('GUS::Model::DoTS::Assembly',undef,1)->getId()]);
  #      print STDERR "DeletedManuallyReviewed: rna.",$rna->getId(),", DT.",$rna->getNASequence('GUS::Model::DoTS::Assembly',undef,1)->getId(),"\n";
  #    }
  #  }

  ##need to also add in the  AssemblySequences that were  marked as chimera..
  foreach my $a ($algoInvo->getChildren('GUS::Model::DoTS::AssemblySequence')) {
    $ct++ if $a->getProcessedCategory() =~ /CHIMERA/;
  }
  
  #  return (scalar(@tmp),$ct,\@reviewed,\@delReviewed);
  return ($ctAssem,$ct);
}


##assigns the identifiers to the new assemblies based upon contained AssemblySequences
#my %mapAss; contains array of assemblysequence ids
#my %mapAssSeq; keys assemblysequence id values the assemblyto which assigned
sub assignIdsForReassembly {
  my %mapNew;
  my %numNew;
  print STDERR "assignIdsForReassembly gene \n" if $debug;
  foreach my $a ($algoInvo->getChildren('GUS::Model::DoTS::Assembly')) {
    if ($a->getId()) {
      ##is one of old ones...already have the old mapping and these should not have any assemblies
      if (scalar($a->getChildren('GUS::Model::DoTS::AssemblySequence')) > 0) {
        print STDERR "ERROR - assignIdsForReassembly: old assembly ".$a->getId()." still has AssemblySequence children (";
        foreach my $c ($a->getChildren('GUS::Model::DoTS::AssemblySequence')) {
          print $c->getId().", ";
        }
        print ")\n";
        ##what should happen in this case???
        return undef;  
      }
      next;
    }
    foreach my $aseq ($a->getChildren('GUS::Model::DoTS::AssemblySequence')) {
      $mapNew{$a->getCacheId()}->{$mapAssSeq{$aseq->getId()}}++ if exists $mapAssSeq{$aseq->getId()}; #$mapNew{newAssemblyCacheId}->{oldAssemblyId}++;
      $numNew{$a->getCacheId()}++;
    }
  }
  ##now what...foreach newAss assign children to most likely old one...don't forget merge split
  ##if old ones left ovet then get the mergeSplit.  Any new assemblies without old ids stay untouched
  ##need to remove the new ids all the way up to RNA...since are not from db should be able to just RNA->removeParent(TU).
  ##note that if any ids are merged into other assemblies that are larger then this id will be merged and any
  ##other assembly that contains some of these ids will be considered new

  my %dealtWith;
  foreach my $new_id (sort { $numNew{$b} <=> $numNew{$a} } keys %numNew) { ##go through the new assemblies..
    my $newAss = $assCache->getCachedAssembly($new_id);
    my @sort = sort { $mapNew{$new_id}->{$b} <=> $mapNew{$new_id}->{$a} } keys %{$mapNew{$new_id}}; ##sort by number of identifiers contained

    my $haveAssigned = 0;
    my $oldAss;
    foreach my $old_id (@sort) {
      next if exists $dealtWith{$old_id};
      if (!$haveAssigned) {     ##have not yet made an assignment
        
        print STDERR "Assigning $numNew{$new_id} children of $new_id to $old_id\n" if $debug;
        $oldAss = $assCache->getFromDbCache('GUS::Model::DoTS::Assembly',$old_id); ##will be in the DbCache as came from database

        $oldAss->addChildren($newAss->getChildren('GUS::Model::DoTS::AssemblySequence')); ##give all the AssemblySequence children to the old id am keeping
        ##so can reverse complement if need be
        print STDERR "assignIds..transferring values\n" if $debug;
        $oldAss->setGappedConsensus($newAss->getGappedConsensus());
        $oldAss->setGappedLength($newAss->getGappedLength());  
        $oldAss->setQualityValues($newAss->getQualityValues());
        $oldAss->setSequence($oldAss->makeConsensus());
        $oldAss->setNumberOfContainedSequences(scalar($oldAss->getChildren('GUS::Model::DoTS::AssemblySequence')));
        print STDERR "assignIds..transferring values complete\n" if $debug;
        $haveAssigned = 1;
        $algoInvo->addChild($oldAss);
        $newAss->markDeleted();
      } else {
        ##create mergesplits here... and markDeleted!!!
        print STDERR "creating mergesplit for $old_id to ".$oldAss->getId()."\n" if $debug;
        my $delAss = $newAss->getFromDbCache('GUS::Model::DoTS::Assembly',$old_id);
        ##add to algoInvo so gets submitted
        $algoInvo->addChild($delAss);
        ##mergesplit...
        my $delRNA = $delAss->getRNA(1);
        my $oldRNA = $oldAss->getRNA(1);
        ##NOTE...need to do merge split for the Assembly not the RNA
        $delAss->addToSubmitList(&createMergeSplit($delRNA,$oldRNA,1)) if $delRNA && $oldRNA;
        $delAss->addToSubmitList(&createMergeSplit($delAss,$oldAss,1));
        ##now mark deleted!!
        &markAssemblyAndChildrenDeleted($delAss);
      }
      $dealtWith{$old_id} = 1;
    }
  }	
  return 1;
}

sub markAssemblyAndChildrenDeleted {
  my($delAss) = @_;
  if ($ctx->{cla}->{no_delete}) { ##just tag that has been deleted...delete later
    $delAss->setDescription('DELETED');
    $delAss->setNumberOfContainedSequences(0);
    $delAss->set('sequence','NULL');
    $delAss->setGappedConsensus('NULL');
    $delAss->setQualityValues('NULL');
    $delAss->markDeleted();
    ##    print STDERR "Tagging Assembly ",$delAss->getId(), " deleted\n";
    return;
  }
  foreach my $ts ( $delAss->getTranslatedAASequences(1,1)) {
    $delAss->markFromAASequenceToSelfDeleted($ts);
    $delAss->addToSubmitList($ts);
    #    $algoInvo->addChild($ts);
  }
  my $rna = $delAss->getRNA(1,1);
  if ($rna) {
    $delAss->markFromRNAToSelfDeleted();
    $delAss->addToSubmitList($rna);
    #    $algoInvo->addChild($rna);  ##so gets submitted in the submit transaction
  }
  $delAss->retrieveAllChildrenExceptAssSeqsFromDB();
  $delAss->markDeleted(1);
}

##want to assemble all trimmed sequences into assembly..sould expect that
##all will assemble into a single assembly or at most two if came from two
##different ends of assembly(s)
sub processTrimmedAssemblySequences {
  return if scalar(@trimmedAssemblySequences) < 1;
  print STDERR "processTrimmedAssemblySequences: starting\n" if $debug;
  open(T,">$tmpLib");
  print T "<CAML>\n\n";
  foreach my $a (@trimmedAssemblySequences) {
    $a->resetAssemblySequence();
    $a->setSinglet(0);
    print T $a->toCAML(1)."\n";
  }
  print T "</CAML>\n";
  close T;
  undef @trimmedAssemblySequences; ##have processed these so don't need....undef just to make sure don't do again.
  my($countContigs,$countSinglets,$singlets,$align) = &processCap4();
  &processAlignment($align);
  print STDERR "processTrimmedAssemblySequences: complete\n" if $debug;
}

sub processAlignment{
  my($align) = @_;
  my $ct = 1;
  foreach my $al (@$align) {
    print STDERR "Processing Alignment ".$ct++.":\n\n$al\n\n" if $debug; ##$al\n" if $debug == 1;
    #				next;
    ##create assembly...parse
    my $newAss = GUS::Model::DoTS::Assembly->new();
    
    push(@trimmedAssemblySequences, $newAss->parseCap4Caml($al));

    ##note following have changed....push removed assemblyseqs above and then asemble together
    ##at end
    ##if AssemblySequences get removed entirely from assembly...just make a singleton assembly
    #    foreach my $das (@delAssSeqs){
    #      $algoInvo->addChild(&makeNewAssembly($das));
    #    }
    #    if(scalar(@delAssSeqs) > 0){
    #      $algoInvo->addChildren(@delAssSeqs);
    #    }

    if ($debug) {
      print STDERR "Raw alignment:\n".$newAss->getCap2Alignment();
      print STDERR $newAss->toXML(0,1);
    }

    ##what if there is something wrong with the parse and there aren't any AssemblySequences?
    ##or is conntig that did not assemble in next iteration
    ##want to just go on to the next one...
    if (! $newAss->getChildren('GUS::Model::DoTS::AssemblySequence')) {
      #      print STDERR "ERROR: Unable to parse cap4 alignment\n$al\n";
      #      die "testing...so dying...\n";
      next;
    }
    
    #				print STDERR "Parsed alignment\n".$newAss->getCap2Alignment()."\n";
    
    ##note that at this point, if there is only a single old ID then should be finished...The
    ##old assembly is still cached so there is no reason to mess with this one further.
    if ($newAss->countContainedAssemblies() == 1 && scalar($newAss->getChildren('GUS::Model::DoTS::AssemblySequence') == 1)) {
      print STDERR "\nAlignment contains only a single old ID and no new ones..moving on to next\n\n" if $debug;
      next;
    }

    
    ##build complete alignment
    print STDERR "Building the alignment...\n" if $debug;
    $newAss->buildCompleteAlignment();
    
    print STDERR "CompleteAlignment:\n" . $newAss->getCap2Alignment() . "\nAlignment finished\n\n" if $debug == 1;
    ##assign id and mark others for deletion...need to delete RNASequence...etc..
    ##cache new one ... will replace old one in cache if inherits old identifier
    ##if new assembly (no old identifiers) then add to $gene..
    if ($newAss->countContainedAssemblies()) {
      my $markDel;
      print STDERR "Assigning identifier to new assembly\n" if $debug;
      ($newAss,$markDel) = $newAss->assignIdentifier();
      print STDERR "Keeping Assembly: ".$newAss->getCacheId()."\n" if $debug;
      ##NOTE:  when marking things for deleteion need to deal with all children with fk constraints.
      ##in this case, the RNA can have a protein child(ren) that is dependent....others will crop up.
      foreach my $d (@{$markDel}) {
        ##need to create mergeSplit...only if both have valid ids..are from db!
        $algoInvo->addChild(&createMergeSplit($d,$newAss,1)) if ($d->getId() && $newAss->getId());
        &markAssemblyAndChildrenDeleted($d);
        
        $d->setParent($algoInvo);
      }
    } else {
      ##is new assembly...create features and protein...
      $newAss->createFeaturesAndAASequence();
      $newAss->setParent($algoInvo); ##so will be submitted...
    }
    ##set last things with this assembly and add to cache
    ##set parent ids relevant to this program..
    ##need these in order to create new RNA if new assembly
    $newAss->setTaxonId($ctx->{cla}->{'taxon_id'});
    $newAss->setSubclassView('Assembly');
    $newAss->setDescription('Not Annotated');
    
    $newAss->cacheAssembly($newAss); ##note that will need to be added as will replace the one currently there..
  }
}

##this needs to be modified extensively in face of manual annotation...
# sub markAssemblyToRNADeleted {
#   my $assembly = shift;
#   print STDERR "Deleting Assembly: ".$assembly->getCacheId()."\n" if $debug;
#   my $r = $assembly->getRNA(1,1);
#   $algoInvo->addChild($r);  ##do this so gets submitted in submit...
#   foreach my $p ($assembly->getRNA()->getChildren('GUS::Model::DoTS::Protein',1)){
#     print "Deleting $p from old $assembly\n" if $debug;
#     $p->retrieveAllChildrenFromDB(1);  ##retrieves all children recursively
#     $p->markDeleted(1);  ##marks self and all children deleted recursively
#   }
#   ##need to mark whole tree deleted up to RNA
#   $assembly->markFromRNAToSelfDeleted();
#   print STDERR "RNA $r ".$assembly->getRNA(undef,1)->getId()." delete state '".$assembly->getRNA(undef,1)->isMarkedDeleted()."' for assembly ".$assembly->getCacheId()."\n" if $debug;
# }

#3runs cap4 on $tmpLib and returns an array of cap4 alignments....
##should also count the number of contigs and singlets to determine if need to rerun...
##also need to return the singlets...
sub processCap4 {
  my($iterate) = @_;
  print STDERR "Processing cap4 alignment\n" if $debug == 1;
  my $countRedo = 0;
  ##first unlink output file so can detect failures..
  unlink "$tmpLib.assem.caml";
  ##deal  with debugging file...
  if ($ctx->{cla}->{debug_assem_file} && !$iterate) {
    die "--debug_assem_file...comment line 749 and following else stmtnt etc to use\n";
    open(C, "$ctx->{cla}->{debug_assem_file}") || die "debug_assem_file: $ctx->{cla}->{debug_assem_file} not found\n";
  }                             # else{
  ##for running on the server...having problem with nodes and want to finish
  #    my $cmd = "$cap4 $tmpLib ".($iterate ? $iterateParams : $ctx->{cla}->{cap4_params});
  ##for rning on the nodes...
  #    my $cmd = "rsh -n $ctx->{cla}->{cap4_machine} 'cd $ctx->{cla}->{directory}; $cap4 $tmpLib ".($iterate ? $iterateParams : $ctx->{cla}->{cap4_params})."'";

  my $cmd;
  ##first should unlink the output file so that will not  be an error  if cap4 fails..
  if ($ctx->{cla}->{cap4_machine} =~ /^s/i) { ##running on server
    $cmd = "$cap4 $tmpLib ".($iterate ? $iterateParams : $ctx->{cla}->{cap4_params});
  } else {
    my $rmCmd = "rsh -n $ctx->{cla}->{cap4_machine} 'cd $ctx->{cla}->{remote_dir}; /bin/rm $tmpLib"."*;'";
    system($rmCmd);
    system("rcp $tmpLib $ctx->{cla}->{cap4_machine}:$ctx->{cla}->{remote_dir}");
    $cmd = "rsh -n $ctx->{cla}->{cap4_machine} 'cd $ctx->{cla}->{remote_dir}; $cap4 $tmpLib ".($iterate ? $iterateParams : $ctx->{cla}->{cap4_params}).($ctx->{cla}->{cap4_machine} =~ /^s/i ? "" : "'");
    $cmd .= "; rcp $ctx->{cla}->{cap4_machine}:$ctx->{cla}->{remote_dir}/$tmpLib.assem.caml .";
  }
  print STDERR "$cmd\n" if $debug;

 REDO:

  system("$cmd");
  print STDERR "parsing cap4 output\n" if $debug == 1;
  open(C, "$tmpLib.assem.caml");
  #  }
  my @align;
  my $singlets;
  my $caml;
  my $countContigs = 0;
  my $countSinglets = 0;
  my $good = 0;
  while (<C>) {
    $good = 1 if /\<CAML\>/;
    $caml .= $_;
    $countSinglets++ if /\<SINGLET/; ## && !/\<CHIMERA/;
    $singlets .= $_ if ($countSinglets && $_ !~ /CAML/);
    $countContigs ++ if /\<CONTIG ID/;
  }
  if (!$good) {
    ##don't have valid output...
    $countRedo++;
    return undef if $countRedo > 2; #try twice for three total times...
    sleep 10;
    goto REDO;
  }
  $caml =~ s/\n//g;
  while ($caml =~ m/(\<CONTIG\s.*?\<\/CONTIG\>)/g) {
    push(@align,$1);
  } 
  ##get the singlets
  while ($caml =~ m/(\<SINGLET.*?\<\/SINGLET\>)/g) {
    my $string = $1;
    if ($string =~ /NAME=\"*\'*(\w+)/) {
      my $id = $1;
      next if $id =~ /^D/;      ##is a singleton Assembly which am already dealing with
      my $as;
      if ($assCache->assemblySequenceIsCached($id)) {
        $as = $assCache->getCachedAssemblySequence($id);
        print STDERR "Getting AssemblySequence $id from cache\n" if $debug;
      } else {
        print STDERR "Getting new Singleton Assembly Sequence $id from db....\n" if $debug;
        $as = GUS::Model::DoTS::AssemblySequence->new({'assembly_sequence_id' => $id});
        $as->retrieveFromDB();
        $assCache->cacheAssemblySequence($as);
      }
      ##what  about if is chimera....
      if ($string =~ /CHIMERA=\"*\'*(\w+)/) {
        $as->resetAssemblySequence();
        $as->setProcessedCategory("CHIMERA=$1");
        $as->removeParent($as->getParent('GUS::Model::DoTS::Assembly')) if $as->getParent('GUS::Model::DoTS::Assembly');
        $algoInvo->addChild($as);
        print LOG "  ".$as->getId()." CHIMERA=$1, iterate='$iterate'\n";
      } else {
        $as->setSinglet(1);
      }
    }
  }
  #  unlink "$tmpLib.assem.caml";
  print STDERR "\nprocessCap4 returns(",scalar(@align),",$countSinglets,\$singlets,\@align)\n\n" if $debug;
  return (scalar(@align),$countSinglets,$singlets,\@align);
}

##sequences could be either from dbEST or GenBank...check if dbEST first then GenBank.
##all sequences should be in AssemblySequence
sub getNewSequences {
  print STDERR "getting ".scalar(@assIds)." new sequences: \(",join(', ',@assIds),"\)\n" if $debug == 1;
  my %del;
  foreach my $id (@assIds) {
    print STDERR "getNewSequences: $id\n" if $debug;
    my $as = $assCache->getFromDbCache('GUS::Model::DoTS::AssemblySequence',$id);
    if (!$as) {
      $as = GUS::Model::DoTS::AssemblySequence->new({'assembly_sequence_id' => $id});
      if (!$as->retrieveFromDB()) {
        ##has been deleted for some reason...
        ##want to remove from the list...
        print STDERR "ERROR: Unable to retrieve AssemblySequence.$id\n";
        $del{$id} = 1;
        next;
      }
    } else {
      print STDERR "  $as: $id retrieved from dbCache\n" if $debug;
    }
    ##mark the sequenceGaps deleted...note that this is important if reassembling and
    ## is a good check if new sequence as could potentially have some sequence gaps that
    ## could screw up the waterworks
    $as->resetAssemblySequence();

    ##don't want to assemble sequences that are shorter than 50 bp...results from setting quality_end = lenssequence.quality_stop
    ##NOTE: here need to add this to the $algoInnvo for submitting and remove the Assembly parent..
    if ($as->getLength() < 50) {
      #      print STDERR "$id ",$as->getParent('GUS::Model::DoTS::ExternalNASequence',1)->getSourceId(),": length too short...\n";
      $del{$id} = 1;
      $algoInvo->addChild($as);
      $as->removeParent($as->getParent('GUS::Model::DoTS::Assembly')) if $as->getParent('GUS::Model::DoTS::Assembly');
      $as->setProcessedCategory('low_quality') unless $as->getProcessedCategory() eq 'low_quality';
    }

    ##lastly, add the AssemblySequences to the cache...
    $assCache->cacheAssemblySequence($as);
  }
  if (scalar(keys%del) > 0) {
    my @tmp = @assIds;
    undef @assIds;
    foreach my $id (@tmp) {
      push(@assIds,$id) unless $del{$id};
    }
  }
}

##this will be used if reassembling...
sub getGusEntriesForReassembly {
  print STDERR "retrieving gus entries \(".join(', ',@oldIds)."\)\n" if $debug;
  if (scalar(@oldIds) == 0) {   ##is a new gene....contains no old dots ids...
    return 1;
  }
  foreach my $na_seq_id (@oldIds) {
    if (! exists $mapAss{$na_seq_id}) { ##new assembly...
      #my $a  = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $na_seq_id,'taxon_id' => $ctx->{cla}->{taxon_id} }); alered by DP 3/08/03-see line below
      my $a  = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $na_seq_id });
      my $haveAssSeqs = 0;
      if ($a->retrieveFromDB()) {
        #        $assCache->cacheAssembly($a);
        ##check this....should use the where hash ref...need to modify these Gene and RNA methods..
        #        print STDERR "getGusEntriesForReassembly new Assembly: ".$a->toXML() if $debug;
        
        ##now want to get the AssemblySequence ids....push onto @assIds 
        ##need to be able at end to track these and assign identifiers correctly so need datastructure..
        foreach my $aseq ($a->getChildren('GUS::Model::DoTS::AssemblySequence',1)) {
          push(@assIds,$aseq->getId());
          $mapAssSeq{$aseq->getId()} = $a->getId();
          push(@{$mapAss{$a->getId()}},$aseq->getId());
          $haveAssSeqs = 1;
        }
        
      }
      if (! $haveAssSeqs) {
        print STDERR "ERROR:  Unable to retrieve Assembly or AssemblySequences for '$na_seq_id': removing from cluster\n";
        ##want to remove the missing id and  continue
        my @tmp = @oldIds;
        undef @oldIds;
        foreach my $id (@tmp) {
          push(@oldIds,$id) unless $id == $na_seq_id;
        }
      }
    } else {
      print STDERR "ERROR: getGusEntriesForReassembly: Already have entry for Assembly $na_seq_id\n" if $debug;
    }
  }
  return 1;
}

##this method retrieves retrieves and caches all assemblies
sub getGusEntries {
  #  my(@oldIds) = @_; ##is global variable
  print STDERR "retrieving gus entries \(".join(', ',@oldIds)."\)\n" if $debug;
  if (scalar(@oldIds) == 0) {   ##is a new gene....contains no old dots ids...
    ##create a dots gene, add a TU child and return it...
    return 1;
  }
  foreach my $na_seq_id (@oldIds) {
    if (! $assCache->assemblyIsCached($na_seq_id)) { ##new assembly...
      ##check this....should use the where hash ref...need to modify these Gene and RNA methods..
      #my $a  = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $na_seq_id,'taxon_id' => $ctx->{cla}->{taxon_id} });altered by DP 3/08/03-see next line
      my $a  = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $na_seq_id});
      if ($a->retrieveFromDB()) {
        #        print STDERR "getGusEntries: ".$a->toXML() if $debug;
        $assCache->cacheAssembly($a);
        print STDERR "Caching: ".$a->getCacheId()."\n" if $debug;
        ##need to add to algoInfo as that is how am currenty submitting and all assemblies will need to 
        #3 be submitted.
        $algoInvo->addChild($a);
        #        $a->retrieveChildrenFromDB('GUS::Model::DoTS::AssemblySequence');  ##may not need to do this but shouldn't impact the spped unless this assembly remains a singleton...
      } else {
        print STDERR "ERROR:  Unable to retrieve Assembly for '$na_seq_id'\n";
        ##want to remove the missing id and  continue
        my @tmp = @oldIds;
        undef @oldIds;
        foreach my $id (@tmp) {
          push(@oldIds,$id) unless $id == $na_seq_id;
        }
      }
    }
  }
  return 1;
}

sub createMergeSplit {
  my($o,$n,$is_merge) = @_;
  print STDERR "Creating MergeSplit: Old:".$o->getId().", New:".$n->getId().", is_merge:'$is_merge'\n" if $debug;
  my $ms = GUS::Model::DoTS::MergeSplit->new({'old_id' => $o->getId(),
                            'new_id' => $n->getId(),
                            'table_id' => $o->getTableIdFromTableName($o->getClassName()),
                            'merge_split_group_id' => $msGroupId,
                            'is_merge' => $is_merge });
  return $ms;
}

sub runDebugFromFile {
  print STDERR "Parsing cap4 file for debugging\n";
  my ($countContigs,$countSinglets,$singlets,$align) = &processCap4();
  if (scalar(@$align) == 0) {   ##didn't get any valid cap2 output
    print STDERR "ERROR: cap4 did not produce any valid output....exiting\n";
    print LOG "cluster_DEBUG: ERROR cap4 did not produce any valid output\n";
    next;
  }
  &processAlignment($align);

  ##It's party time....iterate with cap4 until freezes down  to stable set of assemblies
  ##  need to determine the rules for this very carefully!!
  ##no! just iterate once but with the MaxOverlaps not set so does full alignment
  
  &iterateAssembly($countContigs,$countSinglets) if $ctx->{cla}->{max_iterations};
  ##now do the submits.... gene will submit all RNAs and also genes that have been merged with it
  #  my ($subRet,$reviewed,$delReviewed) = &submitUpdatedAssemblies();
  #  my($assCt,$assSeqCt) = &countAssembliesAndAssemblySequences();
  #  if($subRet){
    
  #    print LOG "cluster_DEBUG finished: $assCt Assemblies from $assSeqCt sequences\n";
  #    if($reviewed || $delReviewed){
  #      foreach my $r (@$reviewed){ print LOG "  REVIEWED: rna.",$r->[0],", DT.",$r->[1],"\n"; }
  #      foreach my $r (@$delReviewed){ print LOG "  REVIEWED deleted: rna.",$r->[0],", DT.",$r->[1],"\n"; }
  #    }
  #    if(!$ctx->{cla}->{commit} || $debug){
  #      print "\ncluster_DEBUG finished: $assCt Assemblies from $assSeqCt sequences\n";
  #      if($reviewed || $delReviewed){
  #        foreach my $r (@$reviewed){ print "  REVIEWED: rna.",$r->[0],", DT.",$r->[1],"\n"; }
  #        foreach my $r (@$delReviewed){ print "  REVIEWED deleted: rna.",$r->[0],", DT.",$r->[1],"\n"; }
  #      }
  foreach my $ass ($algoInvo->getChildren('GUS::Model::DoTS::Assembly')) {
    #            print STDERR $ass->toXML(0,1);
    #            print STDERR $ass->toCAML(1);
    print STDERR "\nFinished Assembly: DT.",$ass->getId(),"\n",$ass->getCap2Alignment();
  }
  #    }
  #  } else {
  #    print LOG "cluster_DEBUG ERROR, NOT SUBMITTED: $assCt Assemblies from $assSeqCt sequences\n"; 
  #  }
  
}

1;

