use strict;

#
# NOTE: some of the Dots Gene pipeline steps are duplicated
#       instead of imported from DotsBuild pipeline steps.
#       This is because the "mgr"s each assume its own distinct buildDir.
#       May consider combining the prop files and make Dots Gene pipeline
#       part of the DotsBuild pipeline in the future.
#

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

  $mgr->runCmd("chmod -R g+w $dotsGeneBuildDir/$buildName/");
}

sub copyPipelineDirToLiniac {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $nickName = $propertySet->getProp('speciesNickname');
  my $dotsRelease = "release".$propertySet->getProp('dotsRelease');
  my $serverPath = $propertySet->getProp('serverPath') . "/$dotsRelease";
  my $liniacServer = $propertySet->getProp('liniacServer');
  my $liniacUser = $propertySet->getProp('liniacUser');
  my $fromDir =   $propertySet->getProp('dotsGeneBuildDir') . "/$dotsRelease";
  my $signal = "dir2liniac";
  return if $mgr->startStep("Copying $fromDir to $serverPath on $liniacServer", $signal);

  $mgr->copyToLiniac($fromDir, $nickName, $liniacServer, $serverPath, $liniacUser);

  $mgr->endStep($signal);
}

sub extractDots {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $name = "finalDots";
  my $signal = "${name}Extract";

  return if $mgr->startStep("Extracting final Dots seqs from GUS", $signal);

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');
  my $taxonId = $propertySet->getProp('taxonId');
  my $species = $propertySet->getProp('speciesFullname');

  my $seqFile = "$mgr->{pipelineDir}/seqfiles/$name.fsa";
  my $logFile = "$mgr->{pipelineDir}/logs/${name}Extract.log";

  my $sql = "select 'DT.' || na_sequence_id, 'length='||length, sequence from dots.Assembly where taxon_id = $taxonId";

  my $cmd = "dumpSequencesFromTable.pl --outputFile $seqFile --gusConfigFile $gusConfigFile  --verbose --idSQL \"$sql\" 2>>  $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub copyDotsToLiniac {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');
  my $liniacServer = $propertySet->getProp('liniacServer');
  my $liniacUser = $propertySet->getProp('liniacUser');

  my $seqfilesDir = "$mgr->{pipelineDir}/seqfiles";
  my $f = "${name}Dots.fsa";

  my $signal = "${name}Dots2liniac";
  return if $mgr->startStep("Copying $seqfilesDir/$f to $serverPath/$mgr->{buildName}/seqfiles on $liniacServer", $signal);

  $mgr->copyToLiniac($seqfilesDir, $f, $liniacServer, "$serverPath/$mgr->{buildName}/seqfiles", $liniacUser);

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
  my $liniacCmdMsg = "submitPipelineJob runDotsGenomeAlign $serverPath/$buildName NUMBER_OF_NODES";
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
  my $signal = "genomeDotsFromLiniac";
  return if $mgr->startStep("Copying genome alignment of DoTS from $liniacServer", $signal);
    
  $mgr->copyFromLiniac($liniacServer,
		       "$serverPath/$buildName/genome/dots-genome/master/mainresult",
		       "per-chr",
		       "$pipelineDir/genome/dots-genome",
		       $liniacUser);
      
  $mgr->endStep($signal);
}

sub deleteBlatAlignment {
    my ($mgr) = @_;
    my $propertySet = $mgr->{propertySet};
 
    my $signal = "deleteBlatAlignment";
 
    return if $mgr->startStep("Deleting Dots BLAT alignments from GUS", $signal);
    
    my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');

    my $sql = "select /*+ RULE */ blat_alignment_id from dots.BlatAlignment where query_taxon_id = $taxonId and query_table_id = 56 and target_table_id = 245 and target_taxon_id = $taxonId and target_external_db_release_id = $genomeId";

    my $cmd = "deleteEntries.pl --table DoTS::BlatAlignment --idSQL \"$sql\" --verbose 2>> $logFile";
    
    $mgr->runCmd($cmd);

    $mgr->endStep($signal);
}

sub cacheEstClonePairs {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};

    my $signal = 'cacheEstClonePairs';

    return if $mgr->startStep("Caching EST clone pairs", $signal);

    my $taxonId = $propertySet->getProp('taxonId');
    my $tempLogin = $propertySet->getProp('tempLogin');

    my @joinCols = ('washu_name', 'image_id');
    foreach my $jc (@joinCols) {
	my $args = "--taxon_id $taxonId --temp_login $tempLogin --join_column $jc";
	$mgr->runPlugin("EstClonePairs", 
			"DoTS::Gene::Plugin::EstClonePairs",
			$args, "chacing EST clone pairs joined by $jc");
    }
    $mgr->endStep($signal);
}

sub loadGenomeAlignments {
  my ($mgr, $queryName, $targetName) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "$queryName-$targetName" . 'LoadBlatALignment';

  return if $mgr->startStep("Loading $queryName-$targetName alignments", $signal);

  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  my $gapTabSpace = $propertySet->getProp('genomeGapLogin');
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $pslDir = "$pipelineDir/genome/$queryName-$targetName/per-chr";

  my $qFile = "$pipelineDir/repeatmask/$queryName/master/mainresult/blocked.seq";
  $qFile = "$pipelineDir/seqfiles/finalDots.fsa" if $queryName =~ /dots/i;

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

  $mgr->runPlugin($signal, "GUS::Common::Plugin::LoadBLATAlignments",
		  $args, "loading genomic alignments of $queryName vs $targetName");
}

sub findGenomicSignals {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  my $tmpLogin = $propertySet->getProp('tempLogin');

  my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin";

  $mgr->runPlugin("FindGenomicSignals", "DoTS::Gene::Plugin::FindGenomicSignals",
		  $args, "searching for splice/polyA signals etc in genomic context of DT alignments");
}

sub createGenomeDotsGene {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --est_pair_cache EstClonePair";

    $mgr->runPlugin('CreateGenomeDotsGene', "DoTS::Gene::Plugin::CreateGenomeDotsGene",
		    $args, "create genome-based Dots Gene from DT alignments");
}

sub computeQualityScore {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --blat_signals_cache BlatAlignmentSignals --subset_selector 1";

    $mgr->runPlugin('ComputeGenomeDotsGeneScore',
		    "DoTS::Gene::Plugin::ComputeGenomeDotsGeneScore",
		    $args, "compute heuristic quality score for genome DoTS genes");
}

sub markAntisense {
    my ($mgr) = @_;
    die "to do";
}

sub mapToSimilarityDotsGene {
    my ($mgr) = @_;
    die "to do";
}

sub integrateWithGus {
    my ($mgr) = @_;
    die "to do";
}

sub makeReleaseFiles {
    my ($mgr) = @_;
    die "to do";
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
  my $prog = `basename $0`;
  chomp $prog;
  print STDERR "usage: $prog propertiesfile\n";
  exit 1;
}

##################################### main ####################################

1;
