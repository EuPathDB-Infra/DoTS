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

sub downloadGenome {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGenome";
  my $doDo = 'isNewGenome';

  return if $mgr->startStep("Downloading genome", $signal, $doDo);

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $downloadSubDir = "$externalDbDir/goldenpath/$genomeVer";

  my $logfile = "$mgr->{pipelineDir}/logs/downloadGenome.log";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $gUrl = $propertySet->getProp('genome_download_url');
  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -P $downloadSubDir $gUrl/chromFa.zip";
  $mgr->runCmd($cmd);
  $mgr->runCmd("unzip $downloadSubDir/chromAgp.zip -d $downloadSubDir");
  $mgr->runCmd("rm $downloadSubDir/chromFa.zip");

  $mgr->endStep($signal);
}

sub createGenomeDbRelease {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $pipelineDir = $mgr->{pipelineDir};

  my $signal = 'createGenomeDbRelease';
  my $doDo = 'isNewGenome';

  # return if $mgr->startStep("Creating ExternalDatabaseRelease in GUS", $signal, 'isNewGenome');

  my $dbrXml = "$pipelineDir/genome/genome_release.xml";
  my $gDbId = $propertySet->getProp('genome_db_id');
  my $gRlsDate = $propertySet->getProp('genome_db_rls_date');
  my $gVer = $propertySet->getProp('genomeVersion');
  my $gVer2 = $propertySet->getProp('genome_db_rls_note');
  my $gUrl = $propertySet->getProp('genome_download_url');
  &makeGenomeReleaseXml ($dbrXml, $gDbId, $gRlsDate, "$gVer: $gVer2", $gUrl);
  
  my $args = "--comment \"add this ext db rel for $gVer ($gVer2)\" --filename $dbrXml";
  $mgr->runPlugin($signal, "GUS::Common::Plugin::UpdateGusFromXML",
		  $args, "Making genome release", $doDo);

  # $mgr->endStep($signal);
}

sub useNewGenomeDbRelease {
  my ($mgr) = @_;

  my $signal = "useNewGenomeDbRelease";
  return if $mgr->startStep("updating prop file to use new genome release", $signal, 'isNewGenome');
  # $mgr->endStep($signal);

  my $propertySet = $mgr->{propertySet};
  my $propFile = $mgr->{propertiesFile};
  my $gDbId = $propertySet->getProp('genome_db_id');
  my $gVer = $propertySet->getProp('genomeVersion');
  my $gVer2 = $propertySet->getProp('genome_db_rls_note');
  my $sql = "SELECT external_database_release_id from SRES.ExternalDatabaseRelease"
      . " WHERE external_database_id = $gDbId AND version = '$gVer: $gVer2'";
  my $task = "1, find the new genome release id ($sql)\n"
      . "2, update genome_db_rls_id in $propFile";
  $mgr->waitForAction($task, $signal);

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

  my $sql = "select 'DT.' || na_sequence_id, 'length='||length, sequence"
      . " from dots.Assembly where taxon_id = $taxonId";

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

sub copyGenomeToLiniac {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet}; 
  my $signal = "genome2liniac";
  my $doDo = 'isNewGenome';

  my $gVer = $propertySet->getProp('genomeVersion');
  my $fromDir = $propertySet->getProp('externalDbDir') . '/goldenpath';
  my $serverPath = $propertySet->getProp('serverExternalDbDir') . '/goldenpath';
  my $liniacServer = $propertySet->getProp('liniacServer');
  my $liniacUser = $propertySet->getProp('liniacUser');
  return if $mgr->startStep("Copying $fromDir/$gVer to $serverPath on $liniacServer", $signal, $doDo);

  $mgr->copyToLiniac($fromDir, $gVer, $liniacServer, $serverPath, $liniacUser);

  $mgr->endStep($signal);
}

sub prepareGenomeAlignmentOnLiniac {
  my ($mgr) = @_;
  my $signal = "prepGenomeAlign";

  my $propertySet = $mgr->{propertySet};
  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = $mgr->{buildName};
  my $chr_lst = "$serverPath/$buildName/genome/dots-genome/input/target.lst";

  my $gaPath = $propertySet->getProp('genome.path');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $srvGDir = $propertySet->getProp('serverExternalDbDir') . '/goldenpath/'. $genomeVer;

  my $liniacCmdMsg = "$gaPath $chr_lst x.fa x.psl -makeOoc=$srvGDir/11.ooc";
  my $liniacLogMsg = "wait a little while for it to finish ";

  return if $mgr->startStep("Making 11.ooc file for over-represented words on liniac", $signal);
  $mgr->endStep($signal);

  $mgr->exitToLiniac($liniacCmdMsg, $liniacLogMsg, 1);
}

sub startGenomicAlignmentOnLiniac {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $liniacServer = $propertySet->getProp('liniacServer');
  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = $mgr->{buildName};
  my $signal = "runGenomeAlign";

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

sub createGenomeVirtualSequence {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "createGenomeVirtualSequences";
  my $doDo = 'isNewGenome';

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeDir = "$externalDbDir/goldenpath/$genomeVer";
  my $gd = CBIL::Util::GenomeDir->new($genomeDir);

  my $taxon = $propertySet->getProp('taxonId');
  my $genome_rel = $propertySet->getProp('genome_db_rls_id');
  my $xml = $mgr->{pipelineDir} . "/genome/genome_virtual_sequence.xml";
  $gd->makeGusVirtualSequenceXml($taxon, $genome_rel, $xml);

  my $args = "--comment \"create Dots.VirtualSequence entries for $genome_rel\" --filename $xml";
  $mgr->runPlugin($signal, "GUS::Common::Plugin::UpdateGusFromXML", $args,
		  "Creating virtual seqs for chroms in new genome", $doDo);
}

sub downloadGenomeGaps {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGenomeGaps";
  my $doDo = 'isNewGenome';

  return if $mgr->startStep("Downloading genome gaps", $signal, $doDo);

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');

  my $downloadSubDir = "$externalDbDir/goldenpath/$genomeVer";
  my $downloadGapDir = "$externalDbDir/goldenpath/$genomeVer" . 'gaps';

  my $logfile = "$mgr->{pipelineDir}/logs/downloadGenomeGaps.log";
  $mgr->runCmd("mkdir -p $downloadGapDir");

  my $gd = CBIL::Util::GenomeDir->new($downloadSubDir);
  my @chrs = $gd->getChromosomes;
  foreach my $chr (@chrs) {
      my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -P $downloadGapDir "
	  . "http://genome.ucsc.edu/goldenPath/$genomeVer/database/chr${chr}_gap.txt.gz";
      $mgr->runCmd($cmd);
  }
  $mgr->runCmd("gunzip $downloadGapDir/*.gz -d $downloadGapDir");
  
  $mgr->endStep($signal);
}

sub insertGenome {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "insertGenome";
  my $doDo = "isNewGenome";

  return if $mgr->startStep("Loading genomic sequences", $signal, $doDo);

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeDir = "$externalDbDir/goldenpath/$genomeVer";
  my $gd = CBIL::Util::GenomeDir->new($genomeDir);
  my @chr_files = $gd->getChromosomeFiles();

  my $gRlsId = $propertySet->getProp('genome_db_rls_id');
  foreach my $cf (@chr_files) {
      my $chr = $1 if $cf =~ /chr(\S+)\.fa/;
      my $args = "--comment \"load genomic seqs for $genomeVer\" "
	  . "--fasta_files $cf --external_database_release_id $gRlsId";
      $mgr->runPlugin("${signal}Chr$chr", "GUS::Common::Plugin::UpdateNASequences",
		      $args, "Loading genomic sequences");
  }

  $mgr->endStep($signal);
}

sub loadGenomeGaps {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "loadGenomeGaps";
  my $doDo = 'isNewGenome';

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $dbi_str = $propertySet->getProp('dbi_str');
  my $genomeDir = "$externalDbDir/goldenpath/$genomeVer";
  my $gapDir = "$externalDbDir/goldenpath/$genomeVer" . 'gaps';
  my $temp_login = $propertySet->getProp('tempLogin');  
  my $temp_password = $propertySet->getProp('tempPassword'); 

  my $args = "--tempLogin \"$temp_login\" --tempPassword \"$temp_password\" "
      . "--dbiStr \"$dbi_str\" --gapDir $gapDir";
  $mgr->runPlugin($signal, "GUS::Common::Plugin::LoadGenomeGaps", $args, "Loading genome gaps", $doDo);
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
	$mgr->runPlugin("${jc}EstClonePairs", 
			"DoTS::Gene::Plugin::EstClonePairs",
			$args, "chacing EST clone pairs joined by $jc");
    }
    $mgr->endStep($signal);
}

sub loadGenomeAlignments {
  my ($mgr, $queryName, $targetName) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "$queryName-$targetName" . 'LoadBlat';

  return if $mgr->startStep("Loading $queryName vs $targetName BLAT", $signal);

  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  my $gapTabSpace = $propertySet->getProp('genomeGapLogin');
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $pslDir = "$pipelineDir/genome/$queryName-$targetName/per-chr";
  my $pslReader = CBIL::Bio::BLAT::PSLDir->new($pslDir);
  my @pslFiles = $pslReader->getPSLFiles();

  my $qFile = "$pipelineDir/repeatmask/$queryName/master/mainresult/blocked.seq";
  $qFile = "$pipelineDir/seqfiles/finalDots.fsa" if $queryName =~ /dots/i;
  my $qTabId = ($queryName =~ /dots/i ? 56 : 57);

  my $args = "--query_file $qFile --keep_best 2 "
      . "--query_table_id $qTabId --query_taxon_id $taxonId "
      . "--target_table_id 245 --target_db_rel_id $genomeId --target_taxon_id $taxonId "
      . "--max_query_gap 5 --min_pct_id 95 max_end_mismatch 10 "
      . "--end_gap_factor 10 --min_gap_pct 90 "
      . "--ok_internal_gap 15 --ok_end_gap 50 --min_query_pct 10";
  unless ($gapTabSpace eq 'n/a') { $args . " --gap_table_space $gapTabSpace"; }
  if ($qTabId == 57) {
      my $gb_db_rel_id = $propertySet->getProp('gb_db_rel_id');
      $args .= " --query_db_rel_id $gb_db_rel_id";
  }

  foreach my $pf (@pslFiles) {
      my $chr = $2 if $pf =~ /(.*\/).*chr(\S+?)\./i;
      $mgr->runPlugin($signal . "Chr$chr" . "Strip", "GUS::Common::Plugin::LoadBLATAlignments",
		      "--blat_files $pf $args --action strip", "strip extra headers in $pf");
      $mgr->runPlugin($signal . "Chr$chr" . "Load", "GUS::Common::Plugin::LoadBLATAlignments",
		      "--blat_files $pf $args --action load", "loading BLAT alignments in $pf");      
  }

  $mgr->runPlugin($signal . "SetBest", "GUS::Common::Plugin::LoadBLATAlignments",
		  "$args --action setbest", "set top alignment status, delete unwanted alignments");      


  $mgr->endStep($signal);
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
    # $args .= ' --skip_chrs 1,2';

    $mgr->runPlugin('CreateGenomeDotsGene', "DoTS::Gene::Plugin::CreateGenomeDotsGene",
		    $args, "create genome-based Dots Gene from DT alignments");
}

sub computeQualityScore {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --blat_signals_cache BlatAlignmentSignals --subset_selector 0";

    # restart
    # $args .= " --skip_chrs 1,2,3";

    $mgr->runPlugin('ComputeGenomeDotsGeneScore',
		    "DoTS::Gene::Plugin::ComputeGenomeDotsGeneScore",
		    $args, "compute heuristic quality score for genome DoTS genes");
}

sub markAntisense {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --genome_dots_gene_cache GenomeDotsGene";

    $mgr->runPlugin('MarkAntisenseGenomeDotsGene',
		    "DoTS::Gene::Plugin::MarkAntisenseGenomeDotsGene",
		    $args, "look for antisense genome DoTS genes and mark db entries");
}

sub mapToSimilarityDotsGene {
    my ($mgr) = @_;

    my $signal = 'mapToSimilarityDotsGene';

    return if $mgr->startStep("mapping genome dots genes to similairity dots genes", $signal);

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args0 = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin "
	. "--temp_map_table GenomeAndSimilarityDotsGeneMap "
	. "--genome_dots_gene_cache GenomeDotsGene --genome_dots_transcript_cache GenomeDotsTranscript ";

    my @modes = ('prep', 'init', 'oneone', 'onemore', 'transfer');
    foreach my $mode (@modes) {
	my $args = $args0 . " --mode $mode";
	$mgr->runPlugin("${mode}MapGenomeAndSimilarityDotsGene",
				"DoTS::Gene::Plugin::MapGenomeAndSimilarityDotsGene",
				$args, "map genome-based and similarity based dots genes via shared DTs");
    }

    $mgr->endStep($signal);
}

sub moveToAllgenes {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --genome_dots_gene_cache GenomeDotsGene --genome_dots_transcript_cache GenomeDotsTranscript";

    $mgr->runPlugin('MoveGenomeDotsGeneToAllgenes',
		    "DoTS::Gene::Plugin::MoveGenomeDotsGeneToAllgenes",
		    $args, "move genome dots gene temp tables to Allgenes schema");
}

sub integrateWithGus {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $args = "--taxon_id $taxonId --genome_db_rls_id $genomeId --temp_login $tmpLogin --genome_dots_gene_cache GenomeDotsGene --genome_dots_transcript_cache GenomeDotsTranscript";

    $mgr->runPlugin('IntegrateGenomeDotsGeneWithGus',
		    "DoTS::Gene::Plugin::IntegrateGenomeDotsGeneWithGus",
		    $args, "integrate genome dots gene info into GUS");
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
