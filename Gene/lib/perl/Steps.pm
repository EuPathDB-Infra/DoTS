use strict;

sub createDotsGenePipelineDir {
  my ($mgr) = @_;

  my $propertySet = $mgr->{propertySet};

  my $buildName = $mgr->{buildName};

  my $dotsGeneBuildDir = $propertySet->getProp('dotsGeneBuildDir');
  my $serverPath = $propertySet->getProp('serverPath');
  my $nodePath = $propertySet->getProp('nodePath');

  return if (-e "$dotsGeneBuildDir/$buildName/seqfiles");

  $mgr->runCmd("mkdir -p $dotsGeneBuildDir/$buildName/seqfiles");

  $mgr->runCmd("chmod -R g+w $dotsGeneBuildDir/$buildName");
}

sub createGenomeDir {
  my ($mgr) = @_;

  my $propertySet = $mgr->{propertySet};
  my $buildName = $mgr->{buildName};

  my $dotsGeneBuildDir = $propertySet->getProp('dotsGeneBuildDir');
  my $serverPath = $propertySet->getProp('serverPath');
  my $nodePath = $propertySet->getProp('nodePath');
  my $gaTaskSize = $propertySet->getProp('genome.taskSize');
  my $gaPath = $propertySet->getProp('genome.path');
  my $gaOptions = $propertySet->getProp('genome.options');
  my $genomeVer = 'goldenpath/' . $propertySet->getProp('genomeVersion');
  my $extGDir = $propertySet->getProp('externalDbDir') . '/' . $genomeVer;
  my $srvGDir = $propertySet->getProp('serverExternalDbDir') . '/'. $genomeVer;

  &makeGenomeDir("dots", "genome", $buildName, $dotsGeneBuildDir, $serverPath,
		   $nodePath, $gaTaskSize, $gaOptions, $gaPath, $extGDir, $srvGDir);

  $mgr->runCmd("chmod -R g+w $dotsBuildDir/$buildName/");
}

sub makeAssemblyDir {
  my ($name, $buildName, $localDir, $mgr) = @_;

  $mgr->runCmd("mkdir -p $localDir/$buildName/assembly/$name/big");
  $mgr->runCmd("mkdir -p $localDir/$buildName/assembly/$name/small");
}

# $name is prevDots or intermedDots
sub extractDots {
  my ($name, $DTprefix, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  $name = "${name}Dots";
  my $signal = "${name}Extract";

  return if $mgr->startStep("Extracting $name assemblies from GUS", $signal);

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');
  my $taxonId = $propertySet->getProp('taxonId');
  my $species = $propertySet->getProp('speciesFullname');

  my $seqFile = "$mgr->{pipelineDir}/seqfiles/$name.fsa";
  my $logFile = "$mgr->{pipelineDir}/logs/${name}Extract.log";

  my $sql = "select $DTprefix na_sequence_id,'[$species]',description,'('||number_of_contained_sequences||' sequences)','length='||length,sequence from dots.Assembly where taxon_id = $taxonId";

  my $cmd = "dumpSequencesFromTable.pl --outputFile $seqFile --gusConfigFile $gusConfigFile  --verbose --idSQL \"$sql\" 2>>  $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub startGenomicAlignmentOnLiniac {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $liniacServer = $propertySet->getProp('liniacServer');
  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = "release".$propertySet->getProp('dotsRelease')."/".$propertySet->getProp('speciesNickname');
  my $signal = "rungenomealign";

  return if $mgr->startStep("Starting genomic alignment", $signal);

  $mgr->endStep($signal);
  my $liniacCmdMsg = "submitPipelineJob runGenomeAlign $serverPath/$buildName NUMBER_OF_NODES";
  my $liniacLogMsg = "monitor $serverPath/$buildName/logs/*.log and xxxxx.xxxx.stdout";

  $mgr->exitToLiniac($liniacCmdMsg, $liniacLogMsg, 1);
}

sub copyGenomeDotsFromLiniac {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $serverPath = $propertySet->getProp('serverPath');
  my $liniacServer = $propertySet->getProp('liniacServer');
  my $liniacUser = $propertySet->getProp('liniacUser');
  my $buildName = $mgr->{'buildName'};
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $signal = "genomeAssemSeqsFromLiniac";
  return if $mgr->startStep("Copying genome alignment of DoTS from $liniacServer", $signal);
    
  $mgr->copyFromLiniac($liniacServer,
		       "$serverPath/$buildName/genome/assemSeqs-genome/master/mainresult",
		       "per-chr",
		       "$pipelineDir/genome/assemSeqs-genome",
		       $liniacUser);
  
  $mgr->runCmd("mkdir -p $pipelineDir/repeatmask/assemSeqs/master/mainresult");
  $mgr->copyFromLiniac($liniacServer,
		       "$serverPath/$buildName/repeatmask/assemSeqs/master/mainresult",
		       "blocked.seq",
		       "$pipelineDir/repeatmask/assemSeqs/master/mainresult",
		       $liniacUser);
    
  $mgr->endStep($signal);
}

sub loadGenomeAlignments {
  my ($mgr, $queryName, $targetName) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "$queryName-$targetName" . 'Load';

  return if $mgr->startStep("Loading $queryName-$targetName alignments", $signal);

  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  my $gapTabSpace = $propertySet->getProp('genomeGapLogin');
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $pslDir = "$pipelineDir/genome/$queryName-$targetName/per-chr";

  my $qFile = "$pipelineDir/repeatmask/$queryName/master/mainresult/blocked.seq";
  my $qTabId = ($queryName =~ /dots/i ? 56 : 57);

  my $args = "--blat_dir $pslDir --query_file $qFile --gap_table_space $gapTabSpace --keep_best 2 "
    . "--query_table_id $qTabId --query_taxon_id $taxonId "
      . "--target_table_id 245 --target_db_rel_id $genomeId --target_taxon_id $taxonId "
	. "--max_query_gap 5 --min_pct_id 95 max_end_mismatch 10 "
	  . "--end_gap_factor 10 --min_gap_pct 90 "
	    . "--ok_internal_gap 15 --ok_end_gap 50 --min_query_pct 10";
  if ($qTabId == 57) {
    my $gb_db_rel_id = $propertySet->getProp('gb_db_rel_id');
    $args .= " --query_db_rel_id $gb_db_rel_id";
  }

  $mgr->runPlugin("LoadBLATAlignments", 
			  "GUS::Common::Plugin::LoadBLATAlignments",
			  $args, "loading genomic alignments of $queryName vs $targetName");
  $mgr->endStep($signal);
}

sub copyDotsToLiniac {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');
  my $liniacServer = $propertySet->getProp('liniacServer');

  my $seqfilesDir = "$mgr->{pipelineDir}/seqfiles";
  my $f = "${name}Dots.fsa";

  my $signal = "${name}Dots2liniac";
  return if $mgr->startStep("Copying $seqfilesDir/$f to $serverPath/$mgr->{buildName}/seqfiles on $liniacServer", $signal);

  $mgr->copyToLiniac($seqfilesDir, $f, $liniacServer, "$serverPath/$mgr->{buildName}/seqfiles");

  $mgr->endStep($signal);
}

sub makeBuildName {
  my ($nickName, $release) = @_;

  return makeBuildNameRls($nickName, $release);
}

sub makeBuildNameRls {
  my ($nickName, $release) = @_;

  return "release${release}/" . $nickName;
}

sub usage {
  print STDERR "usage:  dotsgenebuild propertiesfile\n";
  exit 1;
}

##################################### main ####################################

1;
