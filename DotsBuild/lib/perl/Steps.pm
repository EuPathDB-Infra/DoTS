use strict;

#use lib "$ENV{GUS_HOME}/lib/perl";
#use GUS::Pipeline::Manager;
#use GUS::Pipeline::MakeTaskDirs;
#use CBIL::Util::PropertySet;
#use File::Basename;
use Bio::SeqIO;
use CBIL::Util::GenomeDir;

sub createDotsPipelineDir {
  my ($mgr) = @_;

  my $propertySet = $mgr->{propertySet};

  my $buildName = $mgr->{buildName};

  my $dotsBuildDir = $propertySet->getProp('dotsBuildDir');
  my $serverPath = $propertySet->getProp('serverPath');
  my $nodePath = $propertySet->getProp('nodePath');
  my $rmTaskSize = $propertySet->getProp('repeatmask.taskSize');
  my $rmPath = $propertySet->getProp('repeatmask.path');
  my $rmOptions = $propertySet->getProp('repeatmask.options');
  my $dangleMax = $propertySet->getProp('repeatmask.dangleMax');
  my $bmTaskSize = $propertySet->getProp('blastmatrix.taskSize');
  my $bsTaskSize = $propertySet->getProp('blastsimilarity.taskSize');
  my $bsBparam = $propertySet->getProp('blastsimilarity.Bparam');
  my $bsVparam = $propertySet->getProp('blastsimilarity.Vparam');
  my $bsEparam = $propertySet->getProp('blastsimilarity.Eparam');
  my $wuBlastBinPathCluster = $propertySet->getProp('wuBlastBinPathCluster');
  my $ncbiBlastBinPathCluster = $propertySet->getProp('ncbiBlastBinPathCluster');

  return if (-e "$dotsBuildDir/$buildName/seqfiles");

  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/seqfiles");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/epcr");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/misc");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/genetrap");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/blastSite");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/downloadSite");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/cluster/initial");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/cluster/intermed");
  $mgr->runCmd("mkdir -p $dotsBuildDir/$buildName/cluster/gene");

  &makeRMDir("assemSeqs", $buildName, $dotsBuildDir,
	     $serverPath, $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);
  &makeRMDir("prevDots", $buildName, $dotsBuildDir,
	     $serverPath,  $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);
  &makeRMDir("unalignedAssemSeqs", $buildName, $dotsBuildDir,
	     $serverPath, $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);
  &makeRMDir("alignedDots", $buildName, $dotsBuildDir,
	     $serverPath,  $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);
  &makeRMDir("intermedDots", $buildName, $dotsBuildDir,
	     $serverPath,  $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);
  &makeRMDir("finalDots", $buildName, $dotsBuildDir,
	     $serverPath,  $nodePath, $rmTaskSize, $rmOptions, $dangleMax, $rmPath);

  &makeMatrixDir("assemSeqs", "assemSeqs", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("prevDots", "assemSeqs", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("prevDots", "prevDots", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("unalignedAssemSeqs", "unalignedAssemSeqs", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("alignedDots", "unalignedAssemSeqs", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("alignedDots", "alignedDots", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("intermedDots", "intermedDots", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);
  &makeMatrixDir("finalDots", "finalDots", $buildName, $dotsBuildDir,
		 $serverPath, $nodePath, $bmTaskSize, $wuBlastBinPathCluster);

  # These parm lists seem to have problems
  #  &makeSimilarityDir("finalDots", "$serverPath/$buildName/seqfiles","nrdb", $buildName, $dotsBuildDir,
  #		     $serverPath, $nodePath, $bsTaskSize,
  #		     $wuBlastBinPathCluster,
  #		     "nrdb.fsa", ,'finalDots.fsa','(\d+)', 'blastx',
  #		     "-wordmask=seg+xnu W=3 T=1000 B=$bsBparam V=$bsVparam E=$bsEparam");
  #  &makeSimilarityDir("finalDots", "$serverPath/$buildName/seqfiles", "prodom", $buildName, $dotsBuildDir,
  #		     $serverPath, $nodePath, $bsTaskSize,
  #		     $wuBlastBinPathCluster,
  #		     "prodom.fsa", ,'finalDots.fsa','(\S+)', 'blastx',
  #		     "-wordmask=seg+xnu W=3 T=1000 B=$bsBparam V=$bsVparam E=$bsEparam");
  #  &makeSimilarityDir("finalDots", "$serverPath/$buildName/seqfiles", "cdd", $buildName, $dotsBuildDir,
  #		     $serverPath, $nodePath, $bsTaskSize,
  #		     $ncbiBlastBinPathCluster,
  #		     "cdd/All",  ,'finalDots.fsa','\w+\|\w+\|\d+\s+(\w+)', 'rpsblast',
  #		     "-a 2 -e .1 -p F");
  &makeSimilarityDir("finalDots", "nrdb", $buildName, $dotsBuildDir,
		     $serverPath, $nodePath, $bsTaskSize,
		     $wuBlastBinPathCluster,
		     "nrdb.fsa", "$serverPath/$buildName/seqfiles", 'finalDots.fsa', '(\d+)', 'blastx',
		     "-wordmask=seg+xnu W=3 T=1000 B=$bsBparam V=$bsVparam E=$bsEparam");
  &makeSimilarityDir("finalDots", "prodom", $buildName, $dotsBuildDir,
		     $serverPath, $nodePath, $bsTaskSize,
		     $wuBlastBinPathCluster,
		     "prodom.fsa", "$serverPath/$buildName/seqfiles", 'finalDots.fsa', '(\S+)', 'blastx',
		     "-wordmask=seg+xnu W=3 T=1000 B=$bsBparam V=$bsVparam E=$bsEparam");
  &makeSimilarityDir("finalDots", "cdd", $buildName, $dotsBuildDir,
		     $serverPath, $nodePath, $bsTaskSize,
		     $ncbiBlastBinPathCluster,
		     "cdd/All", "$serverPath/$buildName/seqfiles", 'finalDots.fsa',
		     '\w+\|\w+\|\d+\s+(\w+)', 'rpsblast', "-a 2 -e .1 -p F");

  &makeAssemblyDir("initial", $buildName, $dotsBuildDir, $mgr);
  &makeAssemblyDir("intermed", $buildName, $dotsBuildDir, $mgr);

  $mgr->runCmd("chmod -R g+w $dotsBuildDir/$buildName");
}

sub createGenomeDir {
  my ($mgr) = @_;
  my $signal = "createGenomeDir";
  return if $mgr->startStep("Creating genome dir", $signal);

  my $propertySet = $mgr->{propertySet};
  my $buildName = $mgr->{buildName};
  my $dotsBuildDir = $propertySet->getProp('dotsBuildDir');
  my $serverPath = $propertySet->getProp('serverPath');
  my $nodePath = $propertySet->getProp('nodePath');
  my $gaTaskSize = $propertySet->getProp('genome.taskSize');
  my $gaPath = $propertySet->getProp('genome.path');
  my $gaOptions = $propertySet->getProp('genome.options');
  my $genomeVer = $propertySet->getProp('genomeDir') . '/'
                  . $propertySet->getProp('genomeVersion');
  my $extGDir = $propertySet->getProp('externalDbDir') . '/' . $genomeVer;
  my $srvGDir = $propertySet->getProp('serverExternalDbDir') . '/'. $genomeVer;

  &makeGenomeDir("assemSeqs", "genome", $buildName, $dotsBuildDir, $serverPath,
		 $nodePath, $gaTaskSize, $gaOptions, $gaPath, $extGDir, $srvGDir);

  $mgr->runCmd("chmod -R g+w $dotsBuildDir/$buildName/");
  $mgr->endStep($signal);
}

sub makeAssemblyDir {
  my ($name, $buildName, $localDir, $mgr) = @_;

  $mgr->runCmd("mkdir -p $localDir/$buildName/assembly/$name/big");
  $mgr->runCmd("mkdir -p $localDir/$buildName/assembly/$name/small");
}

sub downloadGenbank {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGenbank";

  return if $mgr->startStep("Downloading Genbank", $signal,'downloadGenbank');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genbankRel = $propertySet->getProp('genbankRel');

  my $downloadSubDir = "$externalDbDir/genbank/$genbankRel";

  my $logfile = "$mgr->{pipelineDir}/logs/downloadGenbank.log";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $rejectFiles = $propertySet->getProp('gbRejectFiles');
  my $acceptFiles = $propertySet->getProp('gbAcceptFiles');

  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -R \"$rejectFiles\" -A \"$acceptFiles\" -P $downloadSubDir ftp://ftp.ncbi.nih.gov/genbank/";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub downloadRefSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadRefSeq";

  return if $mgr->startStep("Downloading RefSeq", $signal,'downloadGenbank');

  my $logfile = "$mgr->{pipelineDir}/logs/downloadRefSeq.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/refseq/$date";
    
  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $subDir = $propertySet->getProp('speciesNickname') eq 'hum' ? 'H_sapiens' : 'M_musculus';

  my $ftpsite = "ftp://ftp.ncbi.nih.gov/refseq/$subDir/mRNA_Prot/";
  my $ftpfile = "*.gbff.gz";

  my $cmd = "wget -t5 -o $logfile  -m -np -nd -nH --cut-dirs=3 -A \"$ftpfile\" -P $downloadSubDir $ftpsite";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub downloadTaxon {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadTaxon";

  return if $mgr->startStep("Downloading Taxon", $signal,'downloadTaxon');

  my $logfile = "$mgr->{pipelineDir}/logs/downloadTaxon.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/taxonomy/$date";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=2 -A \"gi_taxid_nucl.dmp.gz,gi_taxid_prot.dmp.gz,taxdump.tar.gz\" -P $downloadSubDir  ftp://ftp.ncbi.nih.gov/pub/taxonomy/";
  $mgr->runCmd($cmd);

  $mgr->runCmd("gunzip $downloadSubDir/taxdump.tar.gz");

  $mgr->runCmd("tar --extract --file $downloadSubDir/taxdump.tar -C $downloadSubDir");

  $mgr->runCmd("gunzip $downloadSubDir/gi_taxid_nucl.dmp.gz");
  $mgr->runCmd("gunzip $downloadSubDir/gi_taxid_prot.dmp.gz");

  $mgr->endStep($signal);

}

sub downloadNRDB {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadNRDB";

  return if $mgr->startStep("Downloading NRDB", $signal, 'downloadNRDB');

  my $logfile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/nrdb/$date";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $cmd = "wget -t5 -m -np -nd -nH -o $logfile --cut-dirs=4 -A \"nr.gz\"  -P $downloadSubDir  ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/;gunzip $downloadSubDir/nr.gz";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub downloadGenome {
  my ($mgr, $genomeSrcIdRegEx) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGenome";

  return if $mgr->startStep("Downloading genome", $signal, 'downloadGenome');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeDir = $propertySet->getProp('genomeDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeUrlPath = $propertySet->getProp('genomeUrlPath');
  my $genomeFile = $propertySet->getProp('genomeFile');

  my $downloadSubDir = "$externalDbDir/$genomeDir/$genomeVer";

  my $logfile = "$mgr->{pipelineDir}/logs/downloadGenome.log";

  $mgr->runCmd("mkdir -p $downloadSubDir");

#  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -P $downloadSubDir "
#	. "http://genome.ucsc.edu/$genomeDir/$genomeVer/bigZips/chromFa.zip";
  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -P $downloadSubDir "
	. $genomeUrlPath . $genomeFile;

  $mgr->runCmd($cmd);

  if ( $genomeFile =~ '\.zip$') {
    $mgr->runCmd("unzip $downloadSubDir/$genomeFile -d $downloadSubDir");
    $mgr->runCmd("rm $downloadSubDir/$genomeFile");
  } elsif ( $genomeFile =~ '\.gz$') {
    $mgr->runCmd("gunzip -f $downloadSubDir/$genomeFile");
  }

  my $gd = CBIL::Util::GenomeDir->new($downloadSubDir, $genomeSrcIdRegEx);

  if ($genomeSrcIdRegEx) {
    $gd->splitUp();
  }

  $mgr->endStep($signal);
}


sub downloadGenomeGaps {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGenomeGaps";

  return if $mgr->startStep("Downloading genome gaps", $signal, 'downloadGenomeGaps');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeDir = $propertySet->getProp('genomeDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeUrlPath = $propertySet->getProp('genomeUrlPath');
  my $genomeFile = $propertySet->getProp('genomeFile');

  my $downloadSubDir = "$externalDbDir/$genomeDir/$genomeVer";
  my $downloadGapDir = "$externalDbDir/$genomeDir/$genomeVer" . 'gaps';

  my $logfile = "$mgr->{pipelineDir}/logs/downloadGenome.log";

  $mgr->runCmd("mkdir -p $downloadGapDir");

  my $gd = CBIL::Util::GenomeDir->new($downloadSubDir);
  my @chrs = $gd->getChromosomes;
  foreach my $chr (@chrs) {
    my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=1 -P $downloadGapDir "
      . "http://genome.ucsc.edu/$genomeDir/$genomeVer/database/chr${chr}_gap.txt.gz";
    $mgr->runCmd($cmd);
  }

  $mgr->runCmd("gunzip $downloadGapDir/*.gz -d $downloadGapDir");
  $mgr->endStep($signal);
}


sub insertTaxon {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/taxonomy/$date";

  my $names = "$downloadSubDir/names.dmp";
  my $nodes = "$downloadSubDir/nodes.dmp";
  my $gencode = "$downloadSubDir/gencode.dmp";
  my $merged = "$downloadSubDir/merged.dmp";
  my $restart = $propertySet->getProp('insertTaxonRestart');

  my $args = "--names $names --nodes $nodes --gencode $gencode --merged $merged --restart $restart";

  $mgr->runPlugin("loadTaxon", "GUS::Common::Plugin::LoadTaxon", $args,
		  "Loading Taxon tables",'downloadTaxon');
}

sub parseGenbank {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseGenBank";

  my $genbankRel = $propertySet->getProp('genbankRel');

  my $gb_db_rel_id = $propertySet->getProp('gb_db_rel_id');

  my @gbFiles = split (/,/, $propertySet->getProp('gbFiles'));

  my $externalDbDir = $propertySet->getProp('externalDbDir');


  foreach my $file (@gbFiles) {

    my $dirAndFile = "$externalDbDir/genbank/$genbankRel/$file";

    my $args = "--gbRel $genbankRel --file $dirAndFile  --db_rel_id $gb_db_rel_id";

    my $signal = "gbParse_${file}";

    $mgr->runPlugin($signal, "GUS::Common::Plugin::GBParser", $args,
		    "Loading GenBank files into GUS", 'insertGenbank');

  }

  foreach my $file (@gbFiles) {

    my $subDir = "gbParse_".$file;

    my $failFiles = "$mgr->{pipelineDir}/plugins/$subDir/gbparserFailures/*.gb";

    my @fileArr = <$failFiles>;

    if ((scalar @fileArr) >= 1) {
      die "There are GenBank entry failures - evaluate and run GBParser manually - then restart dotsbuild\n";
    }
  }

  $mgr->endStep($signal);

}

sub parseRefSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseRefSeq";

  my $refseqRel = $propertySet->getProp('refseqRel');
  
  my $refseq_rel_id = $propertySet->getProp('refseq_rel_id');
  
  my $refseqFile = $propertySet->getProp('refseqFile');
  
  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $dirAndFile = "$externalDbDir/refseq/$date/$refseqFile";
    
  my $args = "--gbRel $refseqRel --file $dirAndFile --db_rel_id $refseq_rel_id";
    
  my $signal = "refseq_${refseqFile}";
    
  $mgr->runPlugin($signal, "GUS::Common::Plugin::GBParser", $args, "Loading RefSeq files into GUS");

  my $failFiles = "$mgr->{pipelineDir}/plugins/$refseqFile/gbparserFailures/*.gb";
  
  my @fileArr = <$failFiles>;
    
  if ((scalar @fileArr) >= 1) {
    die "There are RefSeq entry failures - evaluate and run GBParser manually - then restart dotsbuild\n";
  }

  $mgr->endStep($signal);

}

sub parsedbEST {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $restart = $propertySet->getProp('dbESTRestart');

  my $taxonIdList = &getTaxonIdList($mgr);

  my $args = "--log $mgr->{pipelineDir}/logs/dbest.log --fullupdate --span 500 --project 'dbEST Parser' --taxon_id_list '$taxonIdList' --restart_number $restart";

  $mgr->runPlugin("loadDbEst", "GUS::Common::Plugin::dbEST", $args,
		  "Loading dbEST files into GUS", 'loadDbEst');
}

sub makeAssemSeqs {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $file = $propertySet->getProp('fileOfRepeats');

  my $taxonIdList = &getTaxonIdList($mgr);

  my $repeatFile = "$externalDbDir/repeat/$file";

  my $phrapDir = $propertySet->getProp('phrapDir');

  my $args = "--taxon_id_list '$taxonIdList' --repeatFile $repeatFile --phrapDir $phrapDir";

  $mgr->runPlugin("makeAssemSeqs",
		  "DoTS::DotsBuild::Plugin::MakeAssemblySequences", $args,
		  "Making assembly table sequences");
}

sub extractAssemSeqs {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonIdList = &getTaxonIdList($mgr);

  my $outputFile = "$mgr->{pipelineDir}/seqfiles/assemSeqs.fsa";
  my $args = "--taxon_id_list '$taxonIdList' --outputfile $outputFile --extractonly";

  $mgr->runPlugin("extractAssemSeqs",
		  "DoTS::DotsBuild::Plugin::ExtractAndBlockAssemblySequences",
		  $args, "Extracting assembly table sequences");

}

sub extractUnalignedAssemSeqs {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonIdList = &getTaxonIdList($mgr);

  my $outputFile = "$mgr->{pipelineDir}/seqfiles/unalignedAssemSeqs.fsa";
  my $args = "--taxon_id_list '$taxonIdList' --outputfile $outputFile --extractonly";

  $mgr->runPlugin("extractUnalignedAssemSeqs",
		  "DoTS::DotsBuild::Plugin::ExtractAndBlockAssemblySequences",
		  $args, "Extracting unaligned assembly sequences");

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

sub copyGenomeToCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet}; 
  my $signal = "genome2cluster";

  my $gVer = $propertySet->getProp('genomeVersion');
  my $gDir = $propertySet->getProp('genomeDir');
  my $fromDir = $propertySet->getProp('externalDbDir') . "/$gDir";
  my $serverPath = $propertySet->getProp('serverExternalDbDir') . "/$gDir";
  return if $mgr->startStep("Copying $fromDir/$gVer to $serverPath on cluster",
			      $signal, 'copyGenomeToCluster');

  $mgr->{cluster}->copyTo($fromDir, $gVer, $serverPath);

  $mgr->endStep($signal);
}



sub copyPipelineDirToCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $nickName = $propertySet->getProp('speciesNickname');
  my $dotsRelease = "release".$propertySet->getProp('dotsRelease');
  my $serverPath = $propertySet->getProp('serverPath') . "/$dotsRelease";
  my $fromDir =   $propertySet->getProp('dotsBuildDir') . "/$dotsRelease";
  my $signal = "dir2cluster";
  return if $mgr->startStep("Copying $fromDir to $serverPath on cluster", $signal);

  $mgr->{cluster}->copyTo($fromDir, $nickName, $serverPath);

  $mgr->endStep($signal);
}

sub prepareGenomeAlignmentOnCluster {
  my ($mgr, $queryName, $targetName) = @_;
  my $signal = "prepGenomeAlign";
  my $doDo = 'downloadGenome';
  my $propertySet = $mgr->{propertySet};
  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = $mgr->{buildName};
  my $seq_lst = "$serverPath/$buildName/genome/$queryName-$targetName/input/target.lst";

  my $gaPath = $propertySet->getProp('genome.path');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeDir = $propertySet->getProp('genomeDir');
  my $srvGDir = $propertySet->getProp('serverExternalDbDir') . '/' . $genomeDir . '/' . $genomeVer;

  my $clusterCmdMsg = "$gaPath $seq_lst x.fa x.psl -makeOoc=$srvGDir/11.ooc";
  my $clusterLogMsg = "wait a little while for it to finish ";

  return if $mgr->startStep("Making 11.ooc file for over-represented words on cluster", $signal, $doDo);
  $mgr->endStep($signal);

  $mgr->exitToCluster($clusterCmdMsg, $clusterLogMsg, 1);
}

sub startGenomicAlignmentOnCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = "release".$propertySet->getProp('dotsRelease')."/".$propertySet->getProp('speciesNickname');
  my $signal = "rungenomealign";

  return if $mgr->startStep("Starting genomic alignment", $signal);

  $mgr->endStep($signal);
  my $clusterCmdMsg = "submitPipelineJob runGenomeAlign $serverPath/$buildName NUMBER_OF_NODES";
  my $clusterLogMsg = "monitor $serverPath/$buildName/logs/*.log and xxxxx.xxxx.stdout";

  $mgr->exitToCluster($clusterCmdMsg, $clusterLogMsg, 1);
}

sub insertGenome {
  my ($mgr, $genomeSrcIdRegEx) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "insertGenome";
print "before startStep(genome)\n";
  return if $mgr->startStep("Inserting genome sequences into GUS", $signal);
print "after startStep(genome)\n";
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeDir = $propertySet->getProp('genomeDir');
  my $dbi_str = $propertySet->getProp('dbi_str');
  my $genomeDir = "$externalDbDir/$genomeDir/$genomeVer";
#  my $gapDir = "$externalDbDir/$genomeDir/$genomeVer" . 'gaps';  huh?
  my $temp_login = $propertySet->getProp('tempLogin');  
  my $temp_password = $propertySet->getProp('tempPassword'); 
  my $gd = CBIL::Util::GenomeDir->new($genomeDir, $genomeSrcIdRegEx);
print "after genomedir->new()\n";
  my $taxon = $propertySet->getProp('taxonId');
  my $genome_rel = $propertySet->getProp('genome_db_rls_id');
  my $dbrXml = "$genomeDir/gusExtDbRel.xml";

  $gd->makeGusVirtualSequenceXml($taxon, $genome_rel, $dbrXml, $genomeSrcIdRegEx);

  my $args = "--comment \"add this ext db rel for $genomeVer\" --filename $dbrXml";
  $mgr->runPlugin("populateVirtualSequence", "GUS::Common::Plugin::UpdateGusFromXML",
		    $args, "Populating Virtual Sequence");


  my $gRlsId = $propertySet->getProp('genome_db_rls_id');

  if ($genomeSrcIdRegEx) {
    my @seqs = $gd->getSequenceFiles();
    foreach my $seqfile (@seqs) {
      my $seq = $1 if $seqfile =~ /([^\/\s]+)\.fa/;
      my $args = "--comment \"load genomic seqs for $genomeVer\" "
	. "--fasta_files $seqfile --external_database_release_id $gRlsId";
      $mgr->runPlugin("${signal}sequence_$seq", "GUS::Common::Plugin::UpdateNASequences",
		      $args, "Loading genomic sequences");
    }
  } else {
    my @chr_files = $gd->getChromosomes();
    foreach my $cf (@chr_files) {
      my $chr = $1 if $cf =~ /chr(\S+)\.fa/;
      my $args = "--comment \"load genomic seqs for $genomeVer\" "
	. "--fasta_files $cf --external_database_release_id $gRlsId";
      $mgr->runPlugin("${signal}Chr$chr", "GUS::Common::Plugin::UpdateNASequences",
		      $args, "Loading genomic sequences");
    }
  }

  $mgr->endStep($signal);
}

sub insertGenomeGaps {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "insertGenomeGaps";
  return if $mgr->startStep("Inserting genome sequences into GUS", $signal, 'insertGenome');
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $genomeVer = $propertySet->getProp('genomeVersion');
  my $genomeDir = $propertySet->getProp('genomeDir');
  my $dbi_str = $propertySet->getProp('dbi_str');
  my $genomeDir = "$externalDbDir/$genomeDir/$genomeVer";
  my $gapDir = "$externalDbDir/$genomeDir/$genomeVer" . 'gaps';
  my $temp_login = $propertySet->getProp('tempLogin');  
  my $temp_password = $propertySet->getProp('tempPassword'); 
  my $gd = CBIL::Util::GenomeDir->new($genomeDir);
  my $taxon = $propertySet->getProp('taxonId');
  my $genome_rel = $propertySet->getProp('genome_db_rls_id');

  my $args = "--tempLogin \"$temp_login\" --tempPassword \"$temp_password\" "
    . "--dbiStr \"$dbi_str\" --gapDir $gapDir";
  $mgr->runPlugin("loadGenomeGaps", "GUS::Common::Plugin::LoadGenomeGaps",
		  $args, "Loading genome gaps", 'insertGenome');

  $mgr->endStep($signal);
}

sub copyGenomeAssemSeqsFromCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $serverPath = $propertySet->getProp('serverPath');
  my $buildName = $mgr->{'buildName'};
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $signal = "genomeAssemSeqsFromCluster";
  return if $mgr->startStep("Copying genome alignment of assemSeqs from cluster", $signal);
    
  $mgr->{cluster}->copyFrom("$serverPath/$buildName/genome/assemSeqs-genome/master/mainresult",
		       "per-seq",
		       "$pipelineDir/genome/assemSeqs-genome");
  
  $mgr->runCmd("mkdir -p $pipelineDir/repeatmask/assemSeqs/master/mainresult");
  $mgr->{cluster}->copyFrom("$serverPath/$buildName/repeatmask/assemSeqs/master/mainresult",
		       "blocked.seq",
		       "$pipelineDir/repeatmask/assemSeqs/master/mainresult");
    
  $mgr->endStep($signal);
}
sub deleteGenomeAlignments {
  my ($mgr,$queryName) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "delete${queryName}GenomeAlignments";
  return if $mgr->startStep("Deleting entries from BLATAlignments",$signal);
  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  my $qTabId = ($queryName =~ /dots/i ? 56 : 57);
  $qTabId = 89 if ($queryName =~ /geneTag/i);
  my $qRelId; 
  if ($qTabId == 57) {
    $qRelId = $propertySet->getProp('gb_db_rel_id');
  }
  if ($qTabId == 89) {
    my @DBs = split(/,/, $propertySet->getProp('geneTrapDbRls'));
    foreach my $db (@DBs) {
      my ($name, $id) = split(/:/, $db);
      $qRelId .= "id,";
    }
    $qRelId =~ s/,$//;
  }
  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";
  my $sql = "select blat_alignment_id from dots.blatalignment where target_table_id = 245 and target_external_db_release_id = $genomeId and query_external_db_release_id in ($qRelId) and query_table_id = $qTabId and query_taxon_id=$taxonId and target_taxon_id = $taxonId";
  my $cmd = "deleteEntries.pl --table DoTS::BLATAlignment --idSQL \"$sql\" --verbose 2>> $logFile";
  $mgr->runCmd($cmd);
  $mgr->endStep($signal);
}

sub loadGenomeAlignments {
  my ($mgr, $queryName, $targetName) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');
  my $genomeId = $propertySet->getProp('genome_db_rls_id');
  #my $gapTabSpace = $propertySet->getProp('genomeGapLogin');
  my $pipelineDir = $mgr->{'pipelineDir'};
  my $pslDir = "$pipelineDir/genome/$queryName-$targetName/per-seq";

  my $qFile;

  if ($queryName =~ /dots/i) {
    $qFile = "$pipelineDir/seqfiles/finalDots.fsa" if $queryName =~ /dots/i;
  } else {
    # copy qFile to /tmp directory to work around a bug in the
    # LoadBLATAlignments plugin's call to FastaIndex
    my $qDir = "/tmp/" . $propertySet->getProp('speciesNickname');
    $mgr->runCmd("mkdir $qDir") if ! -d $qDir;
    $qFile = $qDir . "/blocked.seq";
    $mgr->runCmd("cp $pipelineDir/repeatmask/$queryName/master/mainresult/blocked.seq $qFile");
  }

  my $qTabId = ($queryName =~ /dots/i ? 56 : 57);
  #--gap_table_space $gapTabSpace
  my $args = "--blat_dir $pslDir --query_file $qFile --keep_best 2 "
    . "--query_table_id $qTabId --query_taxon_id $taxonId "
      . "--target_table_id 245 --target_db_rel_id $genomeId --target_taxon_id $taxonId "
	. "--max_query_gap 5 --min_pct_id 95 max_end_mismatch 10 "
	  . "--end_gap_factor 10 --min_gap_pct 90 "
	    . "--ok_internal_gap 15 --ok_end_gap 50 --min_query_pct 10";
  if ($qTabId == 57) {
    my $gb_db_rel_id = $propertySet->getProp('gb_db_rel_id');
    $args .= " --query_db_rel_id $gb_db_rel_id";
  }

  $mgr->runPluginNoCommit("LoadBLATAlignments",
			  "GUS::Common::Plugin::LoadBLATAlignments",
			  $args, "loading genomic alignments of $queryName vs $targetName");
}

sub clusterByGenome {
    my ($mgr, $name) = @_;
    my $propertySet = $mgr->{propertySet};

    my $pipelineDir = $mgr->{'pipelineDir'};
    my $taxonId = $propertySet->getProp("taxonId");
    my $extDbRelId = $propertySet->getProp("genome_db_rls_id");
    my $gb_db_rel_id = $propertySet->getProp('gb_db_rel_id');

    my $outputFile = "$pipelineDir/cluster/$name/cluster.out";
    my $logFile = "$pipelineDir/logs/${name}Cluster.log";

    my $args = "--stage $name --taxon_id $taxonId --query_db_rel_id $gb_db_rel_id "
	. "--target_db_rel_id $extDbRelId --out $outputFile --sort 1";
    # $args .= " --test_chr 5";

    $mgr->runPlugin("ClusterByGenome", 
		    "DoTS::DotsBuild::Plugin::ClusterByGenome",
		    $args, "$name clustering by genome alignment");

}


sub qualityStart {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $date = $propertySet->getProp('qualityStartDate');
  my $taxonIdList = &getTaxonIdList($mgr);
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $file = $propertySet->getProp('fileOfRepeats');
  my $repeatFile = "$externalDbDir/repeat/$file";
  my $phrapDir = $propertySet->getProp('phrapDir');
  
  my $sql = "select a.assembly_sequence_id from dots.assemblysequence a, dotsver.externalnasequencever v where a.na_sequence_id = v.na_sequence_id and v.version_date > $date and v.taxon_id in ($taxonIdList)";

  my $args = "--idSQL \"$sql\" --repeatFile $repeatFile --phrapDir $phrapDir";

  $mgr->runPlugin("setQualityStart", 
		   "DoTS::DotsBuild::Plugin::SetAssSeqQualStartStop",
		   $args, "setting quality start in AssemblySequence",'runQualityStart');

}


sub startBlastMatricesOnCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');
  my $signal = "runmatrices";
  return if $mgr->startStep("Starting blast matrices", $signal);

  $mgr->endStep($signal);

  my $script = "runInitialMatrices";
  if ($propertySet->getProp('firstTime') eq "yes") {
    $script = "runAssemAssemMatrices";
  }

  my $clusterCmdMsg = "submitPipelineJob $script $serverPath/$mgr->{buildName} NUMBER_OF_NODES";
  my $clusterLogMsg = "monitor $serverPath/$mgr->{buildName}/logs/*.log and xxxxx.xxxx.stdout";

  $mgr->exitToCluster($clusterCmdMsg, $clusterLogMsg, 1);
}

sub copyBlastMatricesFromCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "matricesFromCluster";
  return if $mgr->startStep("Copying matrices from cluster", $signal);

  $mgr->{cluster}->copyFrom("$serverPath/$mgr->{buildName}/repeatmask/assemSeqs/master/mainresult",
			    "blocked.err",
			    "$mgr->{pipelineDir}/repeatmask/assemSeqs");

  my @names = ("assemSeqs-assemSeqs", "prevDots-assemSeqs", "prevDots-prevDots");
  if ($propertySet->getProp('firstTime') eq 'yes') {
    @names = ("assemSeqs-assemSeqs");
  }

  foreach my $name (@names) {
    $mgr->{cluster}->copyFrom("$serverPath/$mgr->{buildName}/matrix/$name/master/mainresult",
			 "blastMatrix.out.gz",
			 "$mgr->{pipelineDir}/matrix/$name");
  }

  $mgr->endStep($signal);
}

sub initialCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  if ($propertySet->getProp('firstTime') eq 'yes') {
    &cluster($mgr, "initial", "assemSeqs-assemSeqs");
  } else {
    &cluster($mgr, "initial", "prevdots-prevDots",
	     "assemSeqs-assemSeqs", "prevDots-assemSeqs");
  }
}

sub cluster {
  my ($mgr, $name, @matrices) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "${name}Cluster";

  return if $mgr->startStep("Clustering $name", $signal);

  my $length = $propertySet->getProp("$signal.length");
  my $percent = $propertySet->getProp("$signal.percent");
  my $logbase = $propertySet->getProp("$signal.logbase");
  my $consistentEnds = $propertySet->getProp("$signal.consistentEnds");
  my $cliqueSzArray = $propertySet->getProp("$signal.cliqueSzArray");
  my $logbaseArray = $propertySet->getProp("$signal.logbaseArray");

  my @matrixFileArray;
  foreach my $matrix (@matrices) {
    push(@matrixFileArray,
	 "$mgr->{pipelineDir}/matrix/$matrix/blastMatrix.out.gz");
  }
  my $matrixFiles = join(",", @matrixFileArray);

  my $ceflag = ($consistentEnds eq "yes")? "--consistentEnds" : "";

  my $outputFile = "$mgr->{pipelineDir}/cluster/$name/cluster.out";
  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $cmd = "buildBlastClusters.pl --lengthCutoff $length --percentCutoff $percent --verbose --files '$matrixFiles' --logBase $logbase --iterateCliqueSizeArray $cliqueSzArray $ceflag --iterateLogBaseArray $logbaseArray --sort > $outputFile 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub splitCluster {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "${name}SplitCluster";

  return if $mgr->startStep("SplitCluster $name", $signal);

  my $clusterFile = "$mgr->{pipelineDir}/cluster/$name/cluster.out";
  my $splitCmd = "splitClusterFile $clusterFile";

  $mgr->runCmd($splitCmd);
  $mgr->endStep($signal);
}

sub deleteAlignmentInfo {
    my ($mgr) = @_;

    my $propertySet = $mgr->{propertySet};
    my $taxonId = $propertySet->getProp('taxonId');
    my $genomeId = $propertySet->getProp('genome_db_rls_id');
    my $tmpLogin = $propertySet->getProp('tempLogin');

    my $signal = "deleteAlignmentinfo";

    return if $mgr->startStep("integrate genome dots gene info into GUS", $signal);

    my $args = "--taxon_id $taxonId --only_delete --genome_db_rls_id $genomeId --temp_login $tmpLogin --genome_dots_gene_cache GenomeDotsGene --genome_dots_transcript_cache GenomeDotsTranscript";

    $mgr->runPlugin($signal . 'Delete', "DoTS::Gene::Plugin::IntegrateGenomeDotsGeneWithGus",
		    $args, "delete alignment info from RNAFeature and children");
}

sub assemble {
  my ($old, $reassemble, $name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "${name}Assemble";

  return if $mgr->startStep("Assemble $name", $signal);

  my $clusterFile = "$mgr->{pipelineDir}/cluster/$name/cluster.out";

  &runAssemblePlugin($clusterFile, "big", $name, $old, $reassemble, $mgr);
  &runAssemblePlugin($clusterFile, "small", $name, $old, $reassemble, $mgr);
  $mgr->endStep($signal);
  my $msg =
    "EXITING.... PLEASE DO THE FOLLOWING:
 1. check for errors in assemble.errLog and sql failures in updateDOTSAssemblies.log
 2. resume when assembly completes (validly) by re-runnning 'dotsbuild $mgr->{propertiesFile}'
";
  print STDERR $msg;
  print $msg;
  $mgr->goodbye($msg);
}

sub runAssemblePlugin {
  my ($file, $suffix, $name, $assembleOld, $reassemble, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');
  my $cap4Dir = $propertySet->getProp('cap4Dir');

  my $reass = $reassemble eq "yes"? "--reassemble" : "";
  my $args = "--clusterfile $file.$suffix $assembleOld $reass --taxon_id $taxonId --cap4Dir $cap4Dir";
  my $pluginCmd = "ga DoTS::DotsBuild::Plugin::UpdateDotsAssembliesWithCap4 --commit $args --comment '$args'";

  my $logfile = "$mgr->{pipelineDir}/logs/${name}Assemble.$suffix.log";

  my $assemDir = "$mgr->{pipelineDir}/assembly/$name/$suffix";
  $mgr->runCmd("mkdir -p $assemDir");
  chdir $assemDir || die "Can't chdir to $assemDir";

  my $cmd = "runUpdateAssembliesPlugin --clusterFile $file.$suffix --pluginCmd \"$pluginCmd\" 2>> $logfile";
  $mgr->runCmdInBackground($cmd);
}

sub reassemble {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "${name}Reassemble";

  return if $mgr->startStep("Reassemble $name", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $sql = "select na_sequence_id from dots.assembly where taxon_id = $taxonId  and (assembly_consistency < 90 or length < 50 or length is null or description = 'ERROR: Needs to be reassembled')";

  my $clusterFile = "$mgr->{pipelineDir}/cluster/$name/cluster.out";

  my $suffix = "reassemble";

  my $old = "";

  my $reassemble = "yes";

  my $cmd = "makeClusterFile --idSQL \"$sql\" --clusterFile $clusterFile.$suffix";

  $mgr->runCmd($cmd);

  &runAssemblePlugin($clusterFile, $suffix, $name, $old, $reassemble, $mgr);

  $mgr->endStep($signal);
  my $msg =
    "EXITING.... PLEASE DO THE FOLLOWING:
 1. resume when reassembly completes (validly) by re-runnning 'dotsbuild $mgr->{propertiesFile}'
";
  print STDERR $msg;
  print $msg;
  $mgr->goodbye($msg);
}

sub deleteAssembliesWithNoAssemblySequences {
  my ($mgr, $name) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $args = "--taxon_id $taxonId";
  
  $mgr->runPlugin("${name}deleteAssembliesWithNoAssSeq", 
		  "DoTS::DotsBuild::Plugin::DeleteAssembliesWithNoAssemblySequences",
		  $args, "Deleting assemblies with no assemblysequences");
  
}


sub matrix {
  my ($name, $mgr) = @_;

  &copyDotsToCluster($name, $mgr);

  &startDotsMatrixOnCluster($name, $mgr);

  #&copyDotsMatrixFromCluster($name, $mgr);
}

sub copyDotsToCluster {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $seqfilesDir = "$mgr->{pipelineDir}/seqfiles";
  my $f = "${name}Dots.fsa";

  my $signal = "${name}Dots2cluster";
  return if $mgr->startStep("Copying $seqfilesDir/$f to $serverPath/$mgr->{buildName}/seqfiles on cluster", $signal);

  $mgr->{cluster}->copyTo($seqfilesDir, $f, "$serverPath/$mgr->{buildName}/seqfiles");

  $mgr->endStep($signal);
}

sub copySeqFileToCluster {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $seqfilesDir = "$mgr->{pipelineDir}/seqfiles";
  my $f = "${name}.fsa";

  my $signal = "${name}2cluster";
  return if $mgr->startStep("Copying $seqfilesDir/$f to $serverPath/$mgr->{buildName}/seqfiles on cluster", $signal);

  $mgr->{cluster}->copyTo($seqfilesDir, $f, "$serverPath/$mgr->{buildName}/seqfiles");

  $mgr->endStep($signal);
}

sub startDotsMatrixOnCluster {
  my ($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "${name}DotsMatrix";
  return if $mgr->startStep("Starting ${name}Dots matrix", $signal);

  $mgr->endStep($signal);

  my $cmd = "run" . ucfirst($signal);

  my $clusterCmdMsg = "submitPipelineJob $cmd $serverPath/$mgr->{buildName} NUMBER_OF_NODES";
  my $clusterLogMsg = "monitor $serverPath/$mgr->{buildName}/logs/*.log and xxxxx.xxxx.stdout";

  $mgr->exitToCluster($clusterCmdMsg, $clusterLogMsg, 0);
}

sub copyDotsMatrixFromCluster {
  my($name, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "${name}MatrixFromCluster";

  return if $mgr->startStep("Copying ${name} matrix from cluster",
			    $signal);

  $mgr->{cluster}->copyFrom("$serverPath/$mgr->{buildName}/matrix/$name/master/mainresult",
			    "blastMatrix.out.gz",
			    "$mgr->{pipelineDir}/matrix/$name");

  $mgr->endStep($signal);
}

sub sortClusters {
  my ($mgr) = @_;

  my $signal = "sortClusters";

  return if $mgr->startStep("Sort clusters into descending order", $signal);

  my $in = "$mgr->{pipelineDir}/cluster/gene/cluster.out";
  my $out = "$mgr->{pipelineDir}/cluster/gene/cluster.out.descending";

  $mgr->runCmd("sort --numeric-sort --field-separator ' ' --key 2.2 --reverse $in > $out");

  $mgr->endStep($signal);
}

sub loadRNAClusters {
  my ($mgr) = @_;
  my $file = "$mgr->{pipelineDir}/cluster/gene/cluster.out.descending";
  my $args = "--sort_desc --clusterfile $file";

  $mgr->runPlugin("loadRNAClusters", "DoTS::DotsBuild::Plugin::MakeRNAClustersForAssemblies", $args,
		  "Loading RNA Clusters");
}

sub deleteGenesWithNoRNA {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $user = $propertySet->getProp('userId');

  my $args = "--user_id $user";

  $mgr->runPlugin("deleteGenesWithNoRNA", "DoTS::DotsBuild::Plugin::DeleteGenesWithNoRNA",$args,"Deleting genes with no rna");

}
sub versionGeneTrapAssembly {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "versionGeneTrapAssembly";

  return if $mgr->startStep("Versioning entries from GeneTrapAssembly", $signal,'loadGeneTrapAssembly');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $taxonId = $propertySet->getProp('taxonId');

  my $userId = $propertySet->getProp('userId');

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $sql = "select gene_trap_assembly_id from dots.genetrapassembly g, dots.assembly a where a.na_sequence_id=g.assembly_na_sequence_id and a.taxon_id=$taxonId";

  my $cmd = "versionEntries.pl --table DoTS.GeneTrapAssembly --idSQL \"$sql\" --tablePK 'gene_trap_assembly_id' --userId $userId --gusConfigFile $gusConfigFile --verbose 2>> $logFile";
    
  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}


sub deleteGeneTrapAssembly {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteGeneTrapAssembly";

  return if $mgr->startStep("Deleting entries from GeneTrapAssembly",
			    $signal, 'loadGeneTrapAssembly');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $taxonId = $propertySet->getProp('taxonId');

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $sql = "select gene_trap_assembly_id from dots.genetrapassembly g, dots.assembly a where a.na_sequence_id=g.assembly_na_sequence_id and a.taxon_id=$taxonId";

  my $cmd = "deleteEntries.pl --table DoTS::GeneTrapAssembly --idSQL \"$sql\" --verbose 2>> $logFile";
    
  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}


sub downloadGeneTrapTags {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadGeneTrapTags";

  return if $mgr->startStep("Download gene trap tags manually", $signal,'loadGeneTrapAssembly');

  $mgr->manualTask("download gene trap tags using bulk download on Entrez",$signal);

  $mgr->endStep($signal);
}



sub  updateGeneTrapTags {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateGeneTrapTags";

  return if $mgr->startStep("Updating gene trap tags in GUS", $signal,'loadGeneTrapAssembly');

  my $genbankRel = $propertySet->getProp('genbankRel');

  my @DBs = split(/,/, $propertySet->getProp('geneTrapDbRls'));

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  foreach my $db (@DBs) {
    my ($file, $id) = split(/:/, $db);

    my $dirAndFile = "$externalDbDir/genetags/genbank/$genbankRel/$file";

    my $args = "--gbRel $genbankRel --file $dirAndFile  --db_rel_id $id";

    my $indivSignal = "gbParse_${file}";

    $mgr->runPlugin($indivSignal, "GUS::Common::Plugin::GBParser", $args,
		    "Loading GenBank gene trap files into GUS", 'insertGenbank');

  }

  foreach my $db (@DBs) {
    my ($file, $id) = split(/:/, $db);

    my $subDir = "gbParse_".$file;

    my $failFiles = "$mgr->{pipelineDir}/plugins/$subDir/gbparserFailures/*.gb";

    my @fileArr = <$failFiles>;

    if ((scalar @fileArr) >= 1) {
      die "There are Genbank gene trap tag  entry failures - evaluate and run GBParser manually - then restart dotsbuild\n";
    }
  }

  $mgr->endStep($signal);

}



sub extractGeneTrapTags {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "extractGeneTags";

  return if $mgr->startStep("Extracting gene trap tags from GUS", $signal, 'loadGeneTrapAssembly');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $taxonId = $propertySet->getProp('taxonId');

  my @DBs = split(/,/, $propertySet->getProp('geneTrapDbRls'));

  foreach my $db (@DBs) {
    my ($name, $id) = split(/:/, $db);
    my $seqFile = "$mgr->{pipelineDir}/genetrap/${name}.fsa";
    my $logFile = "$mgr->{pipelineDir}/logs/geneTrapTag${name}.log";

    my $sql = "select na_sequence_id,sequence from dots.ExternalNASequence where taxon_id = $taxonId and external_database_release_id = $id";

    my $cmd = "dumpSequencesFromTable.pl --outputFile $seqFile --verbose --gusConfigFile $gusConfigFile  --idSQL \"$sql\" 2>>  $logFile";
    
    $mgr->runCmd($cmd);
  }
  
  $mgr->endStep($signal);
}

sub formatFinalDots {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "formatFinalDots";

  return if $mgr->startStep("Formatting finalDots.fsa for blast", $signal, 'loadGeneTrapAssembly');

  my $blastBinDir = $propertySet->getProp('wuBlastBinPath');

  my $outputFile1  = "$mgr->{pipelineDir}/seqfiles/finalDots.fsa";
  my $fastalink1 = "$mgr->{pipelineDir}/blastSite/finalDots";

  $mgr->runCmd("ln -s $outputFile1 $fastalink1");
  $mgr->runCmd("$blastBinDir/xdformat -n $fastalink1");
  $mgr->runCmd("rm -rf $fastalink1");
  $mgr->endStep($signal);
}


sub blastGeneTrapTags {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "blastGeneTrapTags";

  return if $mgr->startStep("Blasting gene trap tags vs final mouse DoTS", $signal, 'loadGeneTrapAssembly');

  my $dotsFile = "$mgr->{pipelineDir}/blastSite/finalDots";

  my $blastBinDir = $propertySet->getProp('wuBlastBinPath');

  my $blastn = "${blastBinDir}/blastn";

  my @DB = split (/,/, $propertySet->getProp('geneTrapDbRls'));

  foreach my $db (@DB) {

    my ($name, $id) = split(/:/, $db);

    my $tagFile = "$mgr->{pipelineDir}/genetrap/${name}.fsa";

    my $dotsRelease = $propertySet->getProp('dotsRelease');

    my $outputDir = "$mgr->{pipelineDir}/genetrap/$name";

    my $logFile = "$mgr->{pipelineDir}/logs/${name}Blast.log";

    my $mkdir = "mkdir $outputDir";

    $mgr->runCmd($mkdir) unless (-e "$outputDir");

    my $cmd = "blastAll.pl --blastn $blastn --seqfile $tagFile --musdots $dotsFile --targetdirlogin $outputDir 2>> $logFile";

    $mgr->runCmd($cmd);
  }

  $mgr->endStep($signal);
}

sub loadGeneTrapAssembly {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "loadGeneTrapAssembly";

  return if $mgr->startStep("Loading gene trap tags vs final mouse DoTS", $signal, 'loadGeneTrapAssembly');

  my @DB = split (/,/, $propertySet->getProp('geneTrapDbRls'));
  
  foreach my $db (@DB) {
    my ($name, $id) = split(/:/, $db);
    my $blastDir = "$mgr->{pipelineDir}/genetrap/$name";
    my $log = "$mgr->{pipelineDir}/logs/load${name}GeneTrapBlast.out";
    my $args = "--external_db_release $id --blast_dir $blastDir --logfile $log";
    $mgr->runPlugin("load${name}GeneTrapBlast", "DoTS::DotsBuild::Plugin::CalculateGeneTrapLinks", $args, "loading blast results for $name gene trap tags",'loadGeneTrapAssembly');
  } 

  $mgr->endStep($signal);
}


sub makeFrameFinder {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
    
  my $taxonId = $propertySet->getProp('taxonId');

  my $idSQL = "select na_sequence_id from dots.assembly where taxon_id = $taxonId";

  my $wordfile = $propertySet->getProp('wordfile'); 

  my $restart = $propertySet->getProp('frameFinderRestart');

  my $ffDir = $propertySet->getProp('frameFinderDir');
  my $dianaDir = $propertySet->getProp('dianaDir');

  my $args = "--wordfile $wordfile --restart $restart --ffdir $ffDir --dianadir $dianaDir --idSQL \"$idSQL\" ";

  $mgr->runPlugin("makeFramefinder", 
		  "DoTS::DotsBuild::Plugin::FrameFinder",
		  $args, "running FrameFinder plugin");
}

sub insertNRDB {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $dbi_str = $propertySet->getProp('dbi_str');

  my $sourceDB = $propertySet->getProp('sourceDB');

  my $nrdb_maketemp = $propertySet->getProp('nrdb_maketemp');

  my $nrdb_plugin = $propertySet->getProp('nrdb_plugin');

  my $nrdb_delete = $propertySet->getProp('nrdb_delete');

  my $maketemp = $nrdb_maketemp eq "yes" ? "--maketemp" : "";

  my $plugin = $nrdb_plugin eq "yes" ? "--plugin" : "";

  my $delete = $nrdb_delete eq "yes" ? "--delete" : "";

  my $gitax = "$externalDbDir/taxonomy/$date/gi_taxid_prot.dmp";

  my $nrdb = "$externalDbDir/nrdb/$date/nr";

  my $temp_login = $propertySet->getProp('tempLogin');

  my $temp_password = $propertySet->getProp('tempPassword');

  my $restart = $propertySet->getProp('nrdbRestart') > 1 ? "--restart " . $propertySet->getProp('nrdbRestart') : "";

  my $nrdbReleaseId = $propertySet->getProp('nrdb_db_rls_id');

  my $args = "--temp_login \"$temp_login\" --sourceDB $sourceDB --temp_password \"$temp_password\" --dbi_str \"$dbi_str\" $restart --gitax $gitax --nrdb $nrdb --extDbRelId $nrdbReleaseId $maketemp $plugin $delete";

  $mgr->runPlugin("loadNRDB", "GUS::Common::Plugin::LoadNRDB", $args, "Loading NRDB", 'downloadNRDB');
}

sub extractMarkers {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "extractMarkers";

  return if $mgr->startStep("Extracting Markers", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $markerFile = "$mgr->{pipelineDir}/epcr/rh.sts";

  my $cmd = "extractMarkers --taxon_id $taxonId --outputFile $markerFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}
sub runEPCR {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
 
  my $signal = "runEPCR";

  return if $mgr->startStep("Running e-PCR", $signal);

  my $ePCRinPath = $propertySet->getProp('ePCRinPath');

  my $seqFile = "$mgr->{pipelineDir}/seqfiles/finalDots.fsa";

  my $markersFile = "$mgr->{pipelineDir}/epcr/rh.sts";
  my $epcrFile = "$mgr->{pipelineDir}/epcr/finalDots.epcr";
  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $cmd = "$ePCRinPath/e-PCR $markersFile $seqFile > $epcrFile 2>> $logFile";
    
  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub deleteEPCR {
 my ($mgr) = @_;
 my $propertySet = $mgr->{propertySet};
 
 my $signal = "deleteEPCR";
 
 return if $mgr->startStep("Deleting EPCR from GUS", $signal);
    
 my $taxonId = $propertySet->getProp('taxonId');

 my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

 my $sql = "select /*+ RULE */ e.epcr_id from dots.epcr e, dots.nasequenceimp n where n.taxon_id = $taxonId and e.na_sequence_id=n.na_sequence_id";

 my $cmd = "deleteEntries.pl --table DoTS::EPCR --idSQL \"$sql\" --verbose 2>> $logFile";
    
 $mgr->runCmd($cmd);

 $mgr->endStep($signal);
}


sub insertEPCR {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "insertEPCR";

  my $epcrFile = "$mgr->{pipelineDir}/epcr/finalDots.epcr";
  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";
  
  my $args = "--idcol string1 --file $epcrFile --dir $mgr->{pipelineDir}/plugins/$signal --log $logFile --maptableid 2782 --seqsubclassview Assembly";

  $mgr->runPlugin($signal, "DoTS::DotsBuild::Plugin::LoadEPCR", 
		  $args, "Inserting EPCR");
}

sub deleteAnatomyPercent {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteAnatomyPercent";

  return if $mgr->startStep("Deleting AssemblyAnatomyPercent entries from GUS", $signal);
    
  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $sql = "select assembly_anatomy_percent_id from dots.assemblyanatomypercent where na_sequence_id in (select  na_sequence_id from dots.assembly where taxon_id = $taxonId)";

  my $cmd = "deleteEntries.pl --table DoTS::AssemblyAnatomyPercent --idSQL \"$sql\" --verbose 2>> $logFile";
    
  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub insertAnatomyPercent {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
    
  my $taxonId = $propertySet->getProp('taxonId');

  my $idSQL = "select na_sequence_id from dots.Assembly where taxon_id = $taxonId";

  my $args = "--idSQL \"$idSQL\" --taxon_id $taxonId --restart";

  $mgr->runPlugin("insertAnatomyPercent", 
		  "DoTS::DotsBuild::Plugin::UpdateAssemblyAnatomyPercent", 
		  $args, 
		  "mapping assemblies onto anatomy_id in AssemblyAnatomypercent");
} 


sub extractNRDB {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $nrdbReleaseId = $propertySet->getProp('nrdb_db_rls_id');

  my $sql = "select aa_sequence_id,'source_id='||source_id,'secondary_identifier='||secondary_identifier,description,'length='||length,sequence from dots.ExternalAASequence where external_database_release_id = $nrdbReleaseId";

  &extractProteinSeqs("nrdb", $sql, $mgr, 'copyNRDBToCluster');
}

sub extractProteinSeqs {
  my ($name, $sql, $mgr, $doItProp) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "${name}Extract";

  return if $mgr->startStep("Extracting $name protein sequences from GUS", $signal, $doItProp);

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $seqFile = "$mgr->{pipelineDir}/seqfiles/$name.fsa";
  my $logFile = "$mgr->{pipelineDir}/logs/${name}Extract.log";

  my $cmd = "dumpSequencesFromTable.pl --gusConfigFile $gusConfigFile  --outputFile $seqFile --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub downloadCDD {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  # produces seqfiles/cdd.tar

  my $signal = "downloadCDD";

  return if $mgr->startStep("Downloading CDD", $signal, 'downloadCDD');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('cddDate');

  my $downloadSubDir = "$externalDbDir/cdd/$date";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $logfile = "$mgr->{pipelineDir}/logs/$signal.log";

  $mgr->runCmd("/bin/rm -fr $downloadSubDir/*");

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $cmd = "wget -t5 -m -np -nd -nH -o $logfile --cut-dirs=3 -A \"cdd.tar.gz\"  -P $downloadSubDir ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub unpackCDD {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "unpackCDD";
  return if $mgr->startStep("Unpacking CDD files", $signal, 'downloadCDD');

  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('cddDate');

  my $downloadSubDir = "$externalDbDir/cdd/$date";

  $mgr->runCmd("gunzip $downloadSubDir/cdd.tar.gz") unless (-e "$downloadSubDir/cdd.tar");

  $mgr->runCmd("tar --extract --file $downloadSubDir/cdd.tar -C $downloadSubDir");

  if (-e "$downloadSubDir/LOAD") {
    $mgr->runCmd("rm -f $downloadSubDir/LOAD");
  }

  $mgr->runCmd("catFiles --fileGlob  '$downloadSubDir/LOAD*.csq' --outFile $downloadSubDir/LOAD");

  if (-e "$downloadSubDir/PFAM") {
    $mgr->runCmd("rm -f $downloadSubDir/PFAM");
  }

  $mgr->runCmd("catFiles --fileGlob '$downloadSubDir/pfam*.csq' --outFile $downloadSubDir/PFAM");

  if (-e "$downloadSubDir/SMART") {
    $mgr->runCmd("rm -f $downloadSubDir/SMART");
  }

  $mgr->runCmd("catFiles --fileGlob '$downloadSubDir/smart*.csq' --outFile $downloadSubDir/SMART");

  if (-e "$downloadSubDir/COG") {
    $mgr->runCmd("rm -f $downloadSubDir/COG");
  }

  $mgr->runCmd("catFiles --fileGlob '$downloadSubDir/COG*.csq' --outFile $downloadSubDir/COG");

  if (-e "$downloadSubDir/KOG") {
    $mgr->runCmd("rm -f $downloadSubDir/COG");
  }

  $mgr->runCmd("catFiles --fileGlob '$downloadSubDir/KOG*.csq' --outFile $downloadSubDir/KOG");

  if (-e "$downloadSubDir/CD") {
    $mgr->runCmd("rm -f $downloadSubDir/CD");
  }

  $mgr->runCmd("catFiles --fileGlob '$downloadSubDir/cd*.csq' --outFile $downloadSubDir/CD");

  $mgr->runCmd("cd $downloadSubDir; cat LOAD PFAM SMART COG KOG CD > All");

  $mgr->runCmd("rm -f $downloadSubDir/cdd.tar");

  $mgr->endStep($signal);
}

sub insertCDD {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('cddDate');

  my $downloadSubDir = "$externalDbDir/cdd/$date";

  my $loadDB = $propertySet->getProp('load_db_rls_id');

  # >gnl|CDD|3794 LOAD_ACT, ACT, small ligand binding domain
  my $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(LOAD[\\_\\w]+)\\,\\s";
  my $regex_name = "^\\S+\\sLOAD[\\_\\w]+\\,\\s([\\w\\_\\s]+)\\,\\s";
  my $regex_desc = "^\\S+\\sLOAD[\\_\\w]+\\,\\s(.*)";

  my $loadArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/LOAD' --external_database_release_id $loadDB --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertLoad",
		  "GUS::Common::Plugin::InsertNewExternalSequences",
		  $loadArgs,
		  'Inserting Load', 'downloadCDD');

  my $pfamDB = $propertySet->getProp('pfam_db_rls_id');

  #>gnl|Pfam|pfam00291 PALP, Pyridoxal-phosphate dependent enzyme
  $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(pfam\\w+)";
  $regex_name = "^\\S+\\spfam\\w+\\,\\s([\\w\\_\\s]+)\\,\\s";
  $regex_desc = "^\\S+\\spfam\\w+\\,\\s(.*)";

  my $pfamArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/PFAM' --external_database_release_id $pfamDB --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertPfam", "GUS::Common::Plugin::InsertNewExternalSequences", $pfamArgs,
		  "Inserting Pfam", 'downloadCDD');

  my $smartDB = $propertySet->getProp('smart_db_rls_id');

  #>gnl|Smart|smart00238 BIR, Baculoviral inhibition of apoptosis protein repeat; Domain found in inhibitor of apoptosis proteins (IAPs) and other proteins. Acts as a direct inhibitor of caspase enzymes
  $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(smart\\w+)";
  $regex_name = "^\\S+\\ssmart\\w+\\,\\s([\\w\\_\\s]+)\\,\\s";
  $regex_desc = "^\\S+\\ssmart\\w+\\,\\s(.*)";

  my $smartArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/SMART' --external_database_release_id $smartDB  --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertSmart",
		  "GUS::Common::Plugin::InsertNewExternalSequences",
		  $smartArgs,
		  "Inserting Smart", 'downloadCDD');

  my $cogDB = $propertySet->getProp('cog_db_rls_id');

  #>gnl|CDD|9880 COG0004, AmtB, Ammonia permease [Inorganic ion transport and metabolism]

  $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(COG\\w+)";
  $regex_name = "^\\S+\\sCOG\\w+\\,\\s([\\w\\_\\s]+)\\,\\s";
  $regex_desc = "^\\S+\\sCOG\\w+\\,\\s(.*)";

  my $cogArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/COG' --external_database_release_id $cogDB  --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertCOG",
		  "GUS::Common::Plugin::InsertNewExternalSequences",
		  $cogArgs,
		  "Inserting COG", 'downloadCDD');

  my $kogDB = $propertySet->getProp('kog_db_rls_id');

  #>gnl|CDD|17806 KOG0007, KOG0007, Splicing factor 3a, subunit 1 [RNA processing and modification]

  $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(KOG\\w+)";
  $regex_name = "^\\S+\\sKOG\\w+\\,\\sKOG\\w+\\,\\s([\\w\\,\\s]+)\\[";
  $regex_desc = "^\\S+\\sKOG\\w+\\,\\sKOG\\w+\\,\\s[\\w\\,\\s]+\\[([\\w\\,\\s]+)\\]";

  my $kogArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/KOG' --external_database_release_id $kogDB  --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertKOG",
		  "GUS::Common::Plugin::InsertNewExternalSequences",
		  $kogArgs,
		  "Inserting KOG", 'downloadCDD');

  my $cdDB = $propertySet->getProp('cd_db_rls_id');

  $regex_src_id = "^\\>\\w+\\|\\w+\\|\\w+\\s(cd\\w+)";
  $regex_name = "^\\S+\\scd\\w+\\,\\s([\\w\\_\\s]+)\\,\\s";
  $regex_desc = "^\\S+\\scd\\w+\\,\\s(.*)";

  my $cdArgs = "--verbose --table_name 'DoTS::MotifAASequence' --sequencefile '$downloadSubDir/CD' --external_database_release_id $cdDB  --regex_source_id \'$regex_src_id\' --regex_name \'$regex_name\'  --regex_desc \'$regex_desc\'";

  $mgr->runPlugin("insertCD",
		  "GUS::Common::Plugin::InsertNewExternalSequences",
		  $cdArgs,
		  "Inserting CD", 'downloadCDD');

}

sub downloadProdom {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadProdom";

  return if $mgr->startStep("Downloading Prodom", $signal, 'downloadProdom');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/prodom/$date";

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $logfile = "$mgr->{pipelineDir}/logs/$signal.log";

  $mgr->runCmd("/bin/rm -f $downloadSubDir/*");

  $mgr->runCmd("mkdir -p $downloadSubDir");

  my $cmd = "wget -t5 -m -np -nd -nH -o $logfile  --cut-dirs=3 -A \"prodom.cons.gz\"  -P $downloadSubDir ftp://ftp.toulouse.inra.fr/pub/prodom/current_release/";

  $mgr->runCmd($cmd);
  $mgr->runCmd("gunzip $downloadSubDir/prodom.cons.gz");

  $mgr->endStep($signal);
}

sub insertProdom {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadProdom";

  return if $mgr->startStep("Inserting Prodom", $signal, 'insertProdom');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $subdir = $propertySet->getProp('prodomRelease');

  my $downloadSubDir = "$externalDbDir/prodom/$subdir";

  my $prodomFile = "$downloadSubDir/prodom.cons";

  $mgr->runCmd ("gunzip ${prodomFile}.gz") unless -e "$downloadSubDir/prodom.cons" ; 

  my $logFile = "$mgr->{pipelineDir}/logs/insertProdom.log";

  my $prodomDB = $propertySet->getProp('prodom_db_rls_id');

  #>PD000001 p2001.3 (5600) KIT(16) KPC1(13) CDC2(13) // KINASE TRANSFERASE ATP-BINDING SERINE/THREONINE-PROTEIN RECEPTOR TYROSINE-PROTEIN 2.7.1.-PHOSPHORYLATION PRECURSOR
  #>CONSENSUS#PD000001 | 380 | pd_PD000001; | (5600)  KINASE TRANSFERASE ATP-BINDING SERINE/THREONINE-PROTEIN RECEPTOR TYROSINE-PROTEIN 2.7.1.- PHOSPHORYLATION DOMAIN SERINE/THREONINE  as of 2003.1
  #my $regex_src_id = "^\\>(\\w+)";
  my $regex_src_id = "^\\>CONSENSUS#(\\w+)";
  my $regex_contained_seqs = "\\((\\d+)";
  #my $regex_name = "(\\(.+)\\s\\/\\/";
  #my $regex_desc = "\\/\\/\\s(.+)";
  my $regex_desc = "\\)\\s+(.+)";

  my $args = "--verbose --table_name DoTS::MotifAASequence --sequencefile $prodomFile --external_database_release_id $prodomDB --regex_source_id \"$regex_src_id\"  --regex_contained_seqs \"$regex_contained_seqs\" --regex_desc \"$regex_desc\"";

  $mgr->runPlugin("insertProdom",
		  "GUS::Common::Plugin::InsertNewExternalSequences", $args,
		  "Inserting Prodom", 'insertProdom');

  $mgr->endStep($signal);
}

sub extractProdom {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $prodomDB = $propertySet->getProp('prodom_db_rls_id');

  my $sql = "select aa_sequence_id,'source_id='||source_id,'secondary_identifier='||secondary_identifier,description,'length='||length,sequence from dots.MotifAASequence where external_database_release_id = $prodomDB";

  &extractProteinSeqs("prodom", $sql, $mgr, 'copyProdomToCluster');
}

sub copyProteinDBsToCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "proteinDBs2Cluster";
  return if $mgr->startStep("Copying NRDB, CDD and Prodom to $serverPath/$mgr->{buildName}/seqfiles on cluster", $signal);

  my $release = "release" . $propertySet->getProp('dotsRelease');
  my $proteinRelease;
  my $proteinDir;
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $dotsBuildDir = $propertySet->getProp('dotsBuildDir');
  my $seqfilesDir = "$dotsBuildDir/$release/$speciesNickname/seqfiles";

  my $copyNRDBToCluster = $propertySet->getProp('copyNRDBToCluster');
  my $f = "nrdb.fsa";
  if ($copyNRDBToCluster eq 'yes') {
    $mgr->{cluster}->copyTo($seqfilesDir, $f,
		       "$serverPath/$mgr->{buildName}/seqfiles");
  } else {
    $proteinRelease = "release" . $propertySet->getProp('proteinRelease');
    $proteinDir = $propertySet->getProp('proteinDir');
    my $linkCmd = "ln $serverPath/$proteinRelease/$proteinDir/seqfiles/$f $serverPath/$release/$speciesNickname/seqfiles/$f";
    $mgr->{cluster}->runCmdOnCluster($linkCmd);
  }

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('cddDate');
  my $downloadSubDir = "$externalDbDir/cdd";
  my $tmpCddDir = "$downloadSubDir/cdd";
  die "$tmpCddDir exists but it shouldn't" if -e $tmpCddDir;
  my $copyCDDToCluster = $propertySet->getProp('copyCDD');
  $mgr->runCmd("mv $downloadSubDir/$date $tmpCddDir") if ($copyNRDBToCluster eq 'yes');
  $f = "cdd";
  if ($copyNRDBToCluster eq 'yes') {
    $mgr->{cluster}->copyTo($downloadSubDir, "cdd",
		       "$serverPath/$mgr->{buildName}/seqfiles");
  }else {
    my $linkCmd = "ln -s $serverPath/$proteinRelease/$proteinDir/seqfiles/$f $serverPath/$release/$speciesNickname/seqfiles/$f";
    $mgr->{cluster}->runCmdOnCluster($linkCmd);
  }
  $mgr->runCmd("mv $tmpCddDir $downloadSubDir/$date") if ($copyNRDBToCluster eq 'yes');

  my $copyProdomToCluster = $propertySet->getProp('copyProdomToCluster');
  $f = "prodom.fsa";
  if ($copyProdomToCluster eq 'yes') {
    $mgr->{cluster}->copyTo($seqfilesDir, $f,
		       "$serverPath/$mgr->{buildName}/seqfiles");
  } else {
    my $linkCmd = "ln $serverPath/$proteinRelease/$proteinDir/seqfiles/$f $serverPath/$release/$speciesNickname/seqfiles/$f";
    $mgr->{cluster}->runCmdOnCluster($linkCmd);
  }

  $mgr->endStep($signal);
}

sub startSimilaritiesOnCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "findSimilarities";
  return if $mgr->startStep("Starting NRDB, CDD and Prodom similarites on cluster", $signal);

  $mgr->endStep($signal);

  my $clusterCmdMsg = "submitPipelineJob runSimilarities $serverPath/$mgr->{buildName} NUMBER_OF_NODES";
  my $clusterLogMsg = "monitor $serverPath/$mgr->{buildName}/logs/*.log and xxxxx.xxxx.stdout";

  $mgr->exitToCluster($clusterCmdMsg, $clusterLogMsg, 1);
}

sub copySimilaritiesFromCluster {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $serverPath = $propertySet->getProp('serverPath');

  my $signal = "copySimilaritiesFromCluster";
  return if $mgr->startStep("Copying protein similarities from cluster",
			    $signal);
  my @names = ("finalDots-nrdb","finalDots-prodom", "finalDots-cdd");
  foreach my $name (@names) {
    $mgr->{cluster}->copyFrom("$serverPath/$mgr->{buildName}/similarity/$name/master/mainresult",
			      "blastSimilarity.out.gz",
			      "$mgr->{pipelineDir}/similarity/$name");
  }
  $mgr->endStep($signal);
}

sub deleteOldSimilarities {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');
  my $nrdbRel = $propertySet->getProp('nrdb_db_rls_id');
  my $prodomRel = $propertySet->getProp('prodom_db_rls_id');
  my $smartRel = $propertySet->getProp('smart_db_rls_id');
  my $pfamRel = $propertySet->getProp('pfam_db_rls_id');
  my $loadRel = $propertySet->getProp('load_db_rls_id');
  my $cogRel = $propertySet->getProp('cog_db_rls_id');
  my $cdRel = $propertySet->getProp('cd_db_rls_id');
  my $kogRel = $propertySet->getProp('kog_db_rls_id');

  my $sql = "select /*+ RULE */ similarity_id from dots.similarity s, dots.assembly a, dots.aasequenceimp aas where s.query_table_id = 56 and s.query_id = a.na_sequence_id and a.taxon_id = $taxonId  and s.subject_table_id in (83,84,277) and s.subject_id = aas.aa_sequence_id and aas.external_database_release_id in ($nrdbRel,$prodomRel,$smartRel,$pfamRel,$loadRel,$cogRel,$cdRel,$kogRel)";

  my $args = "--idSQL \"$sql\" ";

  $mgr->runPlugin("deleteProteinSims",
		  "GUS::Common::Plugin::DeleteSimilarities", $args,
		  "Deleting old protein similarities");
}

sub insertProteinSimilarities {
  my ($name, $subjectTable, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $file = "$mgr->{pipelineDir}/similarity/$name/blastSimilarity.out.gz";
  $file .= ".correctPK" if ($name eq "finalDots-cdd");
  my $prop = "iPSRestart_$name";
  my $restart = $propertySet->getProp($prop);

  my $args = "--file $file --restartAlgInvs $restart --queryTable DoTS::Assembly --subjectTable $subjectTable --subjectsLimit 50 --hspsLimit 10";

  $mgr->runPlugin("loadSims_$name",
		  "GUS::Common::Plugin::LoadBlastSimFast", $args,
		  "Loading $name similarities");
}

sub substituteCDDPKs {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $smartDB = $propertySet->getProp('smart_db_rls_id');

  my $pfamDB = $propertySet->getProp('pfam_db_rls_id');

  my $loadDB = $propertySet->getProp('load_db_rls_id');

  my $cogDB = $propertySet->getProp('cog_db_rls_id');

  my $cdDB = $propertySet->getProp('cd_db_rls_id');

  my $kogDB = $propertySet->getProp('kog_db_rls_id');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $signal = "substituteCDDPKs";

  return if $mgr->startStep("Substitituting primary keys into CDD similarities file",
			    $signal);

  my $file = "$mgr->{pipelineDir}/similarity/finalDots-cdd/blastSimilarity.out.gz";
  $mgr->runCmd("gunzip $file");
  my $infile = "$mgr->{pipelineDir}/similarity/finalDots-cdd/blastSimilarity.out";
  my $outfile = "$mgr->{pipelineDir}/similarity/finalDots-cdd/blastSimilarity.out.gz.correctPK";
  my $logfile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $sql = "select source_id,aa_sequence_id from dots.motifaasequence where external_database_release_id in ($smartDB,$pfamDB,$loadDB,$cogDB,$cdDB,$kogDB)";

  $mgr->runCmd("cat $infile | substitutePrimKeysInSimilarity --subjectSQL \"$sql\" --verbose --gusConfigFile $gusConfigFile > $outfile 2>> $logfile");

  $mgr->endStep($signal);
}

sub assignSequenceDescription {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $nrdbReleaseId = $propertySet->getProp('nrdb_db_rls_id');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $assignDescriptionRestart = $propertySet->getProp('assignDescriptionRestart');


  my $restart = $assignDescriptionRestart > 1 ? "--restart \"select na_sequence_id from dots.assembly where row_alg_invocation_id in ($assignDescriptionRestart)\"" : "";

  my $dotsMGIfile = "$mgr->{pipelineDir}/misc/${prefix}_bestNRDBHits.dat";
  my $sql = "select na_sequence_id,na_sequence_id from dots.Assembly where taxon_id = $taxonId";

  my $args = "--update_rna_descriptions --copy_manual_descriptions --taxon_id $taxonId --table DoTS::Assembly --query_table DoTS::Assembly $restart --doNotVersion --nrdb_ext_db_rls_id $nrdbReleaseId --dots_mgi_file $dotsMGIfile  --idSQL \"$sql\"";

  $mgr->runPlugin("assignSeqDescrip",
		  "DoTS::DotsBuild::Plugin::AssignSequenceDescription",
		  $args,
		  "Assigning sequence descriptions from MGI");
}

sub deleteIndexWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteIndexWords";

  return if $mgr->startStep("Deleting index_word_link_ids for assemblies from GUS", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $sql = "select index_word_link_id from dots.indexwordlink l,dots.assembly a where a.taxon_id = $taxonId and a.na_sequence_id = l.target_id and l.target_table_id = 56";

  my $cmd = "deleteEntries.pl --table DoTS::IndexWordLink --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub makeIndexWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $sql = "select na_sequence_id,description from dots.Assembly where taxon_id = $taxonId";

  my $args = "--attribute description --restart --table DoTS::Assembly --idSQL \"$sql\"";

  $mgr->runPlugin("makeIndexWords",
		  "GUS::Common::Plugin::MakeIndexWordLink", $args,
		  "Index assembly description words & add entries to IndexWordLink");
}


sub deleteNRDBIndexWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteNRDBIndexWords";

  return if $mgr->startStep("Deleting nrdb index word links from GUS",
			    $signal, 'downloadNRDB');

  my $buildDate = $propertySet->getProp('buildDate');
  my $nrdbReleaseId = $propertySet->getProp('nrdb_db_rls_id');

  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $sql = "select index_word_link_id from dots.indexwordlink l,dots.externalaasequence e where e.aa_sequence_id = l.target_id and l.target_table_id = 83 and e.external_database_release_id = $nrdbReleaseId  and e.modification_date > '$buildDate'";

  my $cmd = "deleteEntries.pl --table DoTS::IndexWordLink --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub deleteMotifIndexWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteMotifIndexWords";

  return if $mgr->startStep("Deleting motif index word links from GUS",
			    $signal, 'downloadCDD');

  my $buildDate = $propertySet->getProp('buildDate');

  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $sql = "select index_word_link_id from dots.indexwordlink l,dots.motifaasequence m where m.aa_sequence_id = l.target_id and l.target_table_id = 277 and m.modification_date > '$buildDate'";

  my $cmd = "deleteEntries.pl --table DoTS::IndexWordLink --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub makeNRDBIndexWordLink {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $nrdbReleaseId = $propertySet->getProp('nrdb_db_rls_id');

  my $sql = "select aa_sequence_id, description from dots.ExternalAASequence where external_database_release_id=$nrdbReleaseId and aa_sequence_id not in (select target_id from  dots.IndexWordLink where target_table_id=83)";

  my $args = "--attribute description --table DoTS::ExternalAASequence --idSQL \"$sql\"";

  $mgr->runPlugin("makeNRDBIndexWordLink",
		  "GUS::Common::Plugin::MakeIndexWordLink", $args,
		  "Making nrdb index word links", 'downloadNRDB');
}

sub makeMotifIndexWordLink {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $prodomReleaseId = $propertySet->getProp('prodom_db_rls_id');
  my $smartReleaseId = $propertySet->getProp('smart_db_rls_id');
  my $loadReleaseId = $propertySet->getProp('load_db_rls_id');
  my $pfamReleaseId = $propertySet->getProp('pfam_db_rls_id');
  my $cogReleaseId = $propertySet->getProp('cog_db_rls_id');
  my $cdReleaseId = $propertySet->getProp('cd_db_rls_id');
  my $kogReleaseId = $propertySet->getProp('kog_db_rls_id');

  my $sql = "select aa_sequence_id, description from dots.MotifAASequence where external_database_release_id in ($prodomReleaseId,$smartReleaseId,$loadReleaseId,$pfamReleaseId,$cogReleaseId,$cdReleaseId,$kogReleaseId) and aa_sequence_id not in (select target_id from  dots.IndexWordLink where target_table_id=277)";

  my $args = "--attribute description --table DoTS::MotifAASequence --idSQL \"$sql\"";

  $mgr->runPlugin("makeMotifIndexWordLink",
		  "GUS::Common::Plugin::MakeIndexWordLink", $args,
		  "Making motif index word links", 'downloadCDD');
}

sub indexSimilarityWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $sql = "select a.na_sequence_id from dots.assembly a where a.taxon_id = $taxonId";

  my $args = "--similarity_table DoTS::MotifAASequence --target_table DoTS::Assembly --idSQL \"$sql\"";

  $args .= " --restart"
    if ($propertySet->getProp('indexSimilarityWordsRestart') eq "yes");

  $mgr->runPlugin("indexSimilarityWords",
		  "GUS::Common::Plugin::MakeIndexWordSimLink", $args,
		  "Indexing words in similarity descriptions");
}

sub indexNRDBWords {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $sql = "select a.na_sequence_id from dots.assembly a where a.taxon_id = $taxonId";

  my $args = "--similarity_table DoTS::ExternalAASequence --target_table DoTS::Assembly --idSQL \"$sql\"";

  $args .= " --restart"
    if ($propertySet->getProp('indexNRDBWordsRestart') eq "yes");

  $mgr->runPlugin("indexNRDBWords",
		  "GUS::Common::Plugin::MakeIndexWordSimLink", $args,
		  "Indexing words in NRDB descriptions");
}

sub makePredictedProteinFile {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "makePredictedProteinFile";
  
  return if $mgr->startStep("Prepare predicted protein file", $signal);
  
  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";
  my $predictedProteinsFile = "${prefix}_predictedProteins.fasta";

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $logFile = "$mgr->{pipelineDir}/logs/${predictedProteinsFile}Extract.log";

  my $cmd = "dumpPreferredAssemblyTranslationSequences.pl --taxon $taxonId --outputFile $mgr->{pipelineDir}/seqfiles/$predictedProteinsFile --gusConfigFile $gusConfigFile  2>>  $logFile";
  
  $mgr->runCmd($cmd);
  
  $mgr->endStep($signal);

}

sub makeProteinChunks {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "makeProteinChunks";

  return if $mgr->startStep("dividing up the protein sequences", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $predictedProteinsFile = "$mgr->{pipelineDir}/seqfiles/${prefix}_predictedProteins.fasta";

  my $proteinChunkDir = "$mgr->{pipelineDir}/misc/proteinSequenceChunks";

  my $cmd = "mkdir $proteinChunkDir";

  $mgr->runCmd($cmd);

  my $cmd ="fasplit --InputFile $predictedProteinsFile --ChunkSize 1000 --OutputFileFormat ${proteinChunkDir}/%04.04d.fa";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub predictTmAndSignalP {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "predictTmAndSignalP";

  return if $mgr->startStep("making TM and signalP predictions", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $predictedProteinsFile = "$mgr->{pipelineDir}/seqfiles/${prefix}_predictedProteins.fasta";

  my $outputfile = "$mgr->{pipelineDir}/misc/${prefix}_predictedProteins";

  my $proteinChunkDir = "$mgr->{pipelineDir}/misc/proteinSequenceChunks";

  my $cmd = "all-predictions $predictedProteinsFile $outputfile $proteinChunkDir 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}


sub parseTMFile {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet}; 

  my $signal = "parseTMFile";

  return if $mgr->startStep("parsing the TM prediction file", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $predictionPath = "$mgr->{pipelineDir}/misc/*.tmhmm.*";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $outputFile = "$mgr->{pipelineDir}/misc/${prefix}_predictedProteins.tmhmm.parsed";

  my $predictedProteinsFile = "$mgr->{pipelineDir}/seqfiles/${prefix}_predictedProteins.fasta";

  my $type = 'tmhmm';
  my $cmd = " parse-predictions --PredictionPath \"$predictionPath\" --AaSequenceFile $predictedProteinsFile --Type $type > $outputFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub parseSGPSignalP {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseSGPSignalPFile";

  return if $mgr->startStep("parsing the signalP SGP prediction files", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $predictionPath = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/*.fa.sgp.Z";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $outputFile = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/${prefix}_predictedProteins.sgp.parsed";

  my $predictedProteinsFile = "$mgr->{pipelineDir}/seqfiles/${prefix}_predictedProteins.fasta";

  my $type = 'sgp';

  my $cmd = " parse-predictions --PredictionPath \"$predictionPath\" --AaSequenceFile $predictedProteinsFile --Type $type > $outputFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}
 
sub parseSGPHMMSignalP {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseSGPHMMSignalPFile";

  return if $mgr->startStep("parsing the signalP SGPHMM prediction files", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $predictionPath = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/*.fa.sgphmm.Z";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $outputFile = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/${prefix}_predictedProteins.sgphmm.parsed";

  my $predictedProteinsFile = "$mgr->{pipelineDir}/seqfiles/${prefix}_predictedProteins.fasta";

  my $type = 'sgphmm';

  my $cmd = " parse-predictions --PredictionPath \"$predictionPath\" --AaSequenceFile $predictedProteinsFile --Type $type > $outputFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub deletePredictedAAFeatures{
  my ($predictionTable, $mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "delete${predictionTable}";

  return if $mgr->startStep("Deleting predicted $predictionTable features from GUS", $signal);

  my $logFile = "$mgr->{pipelineDir}/logs/$signal.log";

  my $taxonId = $propertySet->getProp('taxonId');

  my $sql = "select p.aa_feature_id from dots.${predictionTable} p, dots.translatedaafeature f, dots.rnafeature r, dots.assembly a where a.taxon_id = $taxonId and a.na_sequence_id = r.na_sequence_id and r.na_feature_id = f.na_feature_id and f.aa_sequence_id = p.aa_sequence_id";

  my $cmd = "deleteAAFeatures.pl --table DoTS::${predictionTable} --idSQL \"$sql\" 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub loadTMHMM {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $inputFile = "$mgr->{pipelineDir}/misc/${prefix}_predictedProteins.tmhmm.parsed";

  my $project_id = $propertySet->getProp('project_id');

  my $args = "--filename $inputFile --project_id $project_id";

  $mgr->runPlugin("loadTMHMM", "DoTS::DotsBuild::Plugin::LoadPredictedAAFeatures", $args, "loading parsed TMHMM predictions into GUS");

}


sub loadSGPSignalP {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $inputFile = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/${prefix}_predictedProteins.sgp.parsed";

  my $project_id = $propertySet->getProp('project_id');

  my $args = "--filename $inputFile --project_id $project_id";

  $mgr->runPlugin("loadSGPSignalP", "DoTS::DotsBuild::Plugin::LoadPredictedAAFeatures", $args, "loading parsed SignalP predictions into GUS");
}

sub loadSGPHMMSignalP {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $inputFile = "$mgr->{pipelineDir}/misc/proteinSequenceChunks/${prefix}_predictedProteins.sgphmm.parsed";

  my $project_id = $propertySet->getProp('project_id');

  my $args = "--filename $inputFile --project_id $project_id";
  
  $mgr->runPlugin("loadSGPHMMSignalP", "DoTS::DotsBuild::Plugin::LoadPredictedAAFeatures", $args, "loading parsed SignalPHMM predictions into GUS");
}

sub assemblyProteinIntegration {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $args = "--taxon_id $taxonId";

  $mgr->runPlugin("integrateAssemblyProtein", "DoTS::DotsBuild::Plugin::AssemblyProteinInstance", $args, "integrating assemblies and proteins"); 
}

sub RNAProteinIntegration {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $sourceDB = $propertySet->getProp('sourceDB');

  my @DBs = split(/,/, $sourceDB);

  my @extRel;

  foreach my $db (@DBs) {
    my ($db_rel_id, $abrev) = split(/:/, $db);
    die "--sourceDB db_rel_id:abreviation pairs must be provided\n" if (!$db_rel_id || !$abrev);
    push (@extRel, $db_rel_id);
  }

  my $relList = join(',', @extRel);

  my $args = "--taxon_id $taxonId --ext_db_rel $relList";

  $mgr->runPlugin("integrateRNAProtein", "DoTS::DotsBuild::Plugin::RNAProteinIntegration", $args, "integrating RNA and proteins");
}

sub setPreferredProtein {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');
  my $refseq_rel_id = $propertySet->getProp('refseq_rel_id');
  my $swissprot_rel_id = $propertySet->getProp('swissprot_rel_id');

  my $args = "--taxon_id $taxonId --refseq_db_rel_id $refseq_rel_id --swissprot_db_rel_id $swissprot_rel_id";

  $mgr->runPlugin("setPreferrredProtein", "DoTS::DotsBuild::Plugin::FindPreferredProtein", $args, "determining preferred protein translation for assemblies and updating dots.proteininstance.is_referenced");
}

sub makeGeneForRNA  {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $args = "--taxon $taxonId";

  $mgr->runPlugin("makeGeneForRNA", "DoTS::DotsBuild::Plugin::MakesGenesForRNA", $args, "making genes for rna with null gene_id");
} 
 
sub getIdsPerAssembly {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "getIdsPerAssembly";

  return if $mgr->startStep("Getting ids per assembly", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $output  = "$mgr->{pipelineDir}/misc/${prefix}_accessionsPerAssembly.dat";

  my $cmd = "getIdsPerAssembly --taxon_id $taxonId --gusConfigFile $gusConfigFile  > $output";

  $mgr->runCmd($cmd);

  my $cmd = "gzip $output";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub getAssembliesPerGene {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "getAssembliesPerGene";

  return if $mgr->startStep("Getting assemblies per gene", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $output  = "$mgr->{pipelineDir}/misc/${prefix}_DTperDG.dat";

  my $cmd = "getAssPerGene --taxon_id $taxonId --gusConfigFile $gusConfigFile  > $output";

  $mgr->runCmd($cmd);

  my $cmd = "gzip $output";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub getmRNAPerAssembly {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "getmRNAPerAssembly";

  return if $mgr->startStep("Getting mRNA per Assembly", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $output  = "$mgr->{pipelineDir}/misc/${prefix}_mRNAaccessionsPerAssembly.dat";

  my $cmd = "getmRNAperAssembly --taxon_id $taxonId --gusConfigFile $gusConfigFile  > $output";

  $mgr->runCmd($cmd);

  my $cmd = "gzip $output";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub makeEpconFastaFiles {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "makeEpconFastaFiles";

  return if $mgr->startStep("Making EpConDB fasta files", $signal);

  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $logFile ="$mgr->{pipelineDir}/logs/${signal}.log";

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $outputFile  = "$mgr->{pipelineDir}/misc/EPConDB_rel${dotsRelease}_${speciesNickname}DoTS";

  my $taxonId = $propertySet->getProp('taxonId');

  my $species = $propertySet->getProp('speciesFullname');

  my $epconDB_anatomy_ids = $propertySet->getProp('epconDB_anatomy_ids');

  my $epconDB_array = $propertySet->getProp('epconDB_array');
  
  my $epconDB_chip = $propertySet->getProp('epconDB_chip');

  my $epconDB_makefile = $propertySet->getProp('epconDB_makefile');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $sql = "select 'DT.' ||a.na_sequence_id, '[ $species ]', a.description ,'('||number_of_contained_sequences || ' sequences)', 'length=' || a.length ,sequence from dots.assembly a, dots.assemblyanatomypercent p where p.anatomy_id in ($epconDB_anatomy_ids) and p.percent>0 and a.na_sequence_id=p.na_sequence_id and p.taxon_id=a.taxon_id and a.taxon_id=$taxonId";

  my $cmd = "dumpSequencesFromTable.pl --verbose --gusConfigFile $gusConfigFile --outputFile $outputFile --idSQL \"$sql\" 2>>  $logFile";

  $mgr->runCmd($cmd);

  my $outputFile  = "$mgr->{pipelineDir}/misc/${epconDB_chip}_rel${dotsRelease}_${speciesNickname}DoTS";

  my $sql = "select 'DT.' ||a.na_sequence_id, '[ $species ]', a.description ,'(' ||number_of_contained_sequences || ' sequences)' , 'length=' ||a.length, sequence from rad3.elementassembly_mv ea , rad3.spot s, dots.assembly a where s.array_id = $epconDB_array and s.element_id = ea.element_id and ea.assembly_na_sequence_id = a.na_sequence_id";

  my $cmd = "dumpSequencesFromTable.pl --verbose --gusConfigFile $gusConfigFile --outputFile $outputFile --idSQL \"$sql\" 2>>  $logFile";

  $mgr->runCmd($cmd) if $epconDB_makefile eq 'yes';

  $mgr->endStep($signal);

}

sub prepareEPConBlastSiteFiles {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "prepareEpConcdBlastSiteFiles";
  return if $mgr->startStep("Preparing EPConDB files for blast site", $signal);
  my $epconDB_makefile = $propertySet->getProp('epconDB_makefile');

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $blastBinDir = $propertySet->getProp('wuBlastBinPath');
  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $outputFile1  = "$mgr->{pipelineDir}/misc/EPConDB_rel${dotsRelease}_${speciesNickname}DoTS";
  my $fastalink1 = "$mgr->{pipelineDir}/blastSite/EPConDB_${speciesNickname}DoTS";

  $mgr->runCmd("ln -s $outputFile1 $fastalink1");
  $mgr->runCmd("$blastBinDir/xdformat -n $fastalink1");
  $mgr->runCmd("rm -rf $fastalink1");
  $mgr->runCmd("gzip $outputFile1");
  

  $mgr->endStep($signal) if $epconDB_makefile eq 'no';

  my $epconDB_chip = $propertySet->getProp('epconDB_chip');
  my $outputFile2  = "$mgr->{pipelineDir}/misc/${epconDB_chip}_rel${dotsRelease}_${speciesNickname}DoTS";
  my $fastalink2 = "$mgr->{pipelineDir}/blastSite/${epconDB_chip}_${speciesNickname}DoTS";

  $mgr->runCmd("ln -s $outputFile2 $fastalink2");
  $mgr->runCmd("$blastBinDir/xdformat -n $fastalink2");
  $mgr->runCmd("rm -rf $fastalink2");
  $mgr->runCmd("gzip $outputFile2");

  $mgr->endStep($signal);

}


sub prepareBlastSiteFiles {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "prepareBlastSiteFiles";

  return if $mgr->startStep("Preparing files for blast site", $signal);

  # make blast files for this species
  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $blastBinDir = $propertySet->getProp('wuBlastBinPath');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";
  my $fastafile = "$mgr->{pipelineDir}/downloadSite/${prefix}.fasta.gz";
  my $cmd = "gunzip $mgr->{pipelineDir}/downloadSite/${prefix}.fasta.gz";
  $mgr->runCmd($cmd);
  $fastafile = "$mgr->{pipelineDir}/downloadSite/${prefix}.fasta";
  my $fastalink = "$mgr->{pipelineDir}/blastSite/${speciesNickname}DoTS";
  $mgr->runCmd("ln -s $fastafile $fastalink");
  $mgr->runCmd("$blastBinDir/xdformat -n $fastalink");

  # make blast files for all species.  combine this build w/ latest of other species
  my $dotsBuildDir = $propertySet->getProp('dotsBuildDir');
  my $otherSpeciesNickname = $speciesNickname eq 'mus'? 'hum' : 'mus';
  my $otherSpeciesRelease = $propertySet->getProp('otherSpeciesRelease');
  my $otherSpeciesBuildName = &makeBuildNameRls($otherSpeciesNickname, $otherSpeciesRelease);
  my $otherPipelineDir = "$dotsBuildDir/$otherSpeciesBuildName";
  my $prefix = "${otherSpeciesNickname}DoTS_rel${otherSpeciesRelease}";

  my $cmd = "gunzip $otherPipelineDir/downloadSite/${prefix}.fasta.gz";
  $mgr->runCmd($cmd);
  my $fastafile2 = "$otherPipelineDir/downloadSite/${prefix}.fasta";
  my $fastalink2 = "$mgr->{pipelineDir}/blastSite/${otherSpeciesNickname}DoTS";
  $mgr->runCmd("ln -s $fastafile2 $fastalink2");
  $mgr->runCmd("cat $fastalink $fastalink2 > $mgr->{pipelineDir}/blastSite/allDoTS");
  $mgr->runCmd("$blastBinDir/xdformat -n $mgr->{pipelineDir}/blastSite/allDoTS");
  $mgr->runCmd("gzip $fastafile2");
  $mgr->runCmd("rm -rf $fastalink2");
  $mgr->runCmd("rm -rf $mgr->{pipelineDir}/blastSite/allDoTS");


  $mgr->runCmd("rm -rf $fastalink");

  $mgr->runCmd("gzip $fastafile");

  $mgr->endStep($signal);

}

sub markBadSeqs {
  my ($mgr) = @_;

  my $file = "$mgr->{pipelineDir}/repeatmask/assemSeqs/master/mainresult/blocked.err";

  my $regex = "\\>(\\d+)";

  my $args = "--filename $file --processed_category repeat --regex_id \"$regex\"";

  $mgr->runPlugin("markBadSeqs",
		  "DoTS::DotsBuild::Plugin::MarkAssemblySequencesBad",
		  $args,
		  "Marking bad assembly table sequences");
}

sub downloadHInvitationalFile {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "downloadHinvitational";

  return if $mgr->startStep("Download H-inviational file", $signal, 'downloadHinvitational');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/h-invitational/$date";

  my $logfile = "$mgr->{pipelineDir}/logs/${signal}.log";

  $mgr->runCmd("mkdir -p $downloadSubDir") unless (-e $downloadSubDir);

  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH  -A \"acc2hinv_id.txt.gz\" -P $downloadSubDir ftp://hinv.ddbj.nig.ac.jp/";

  $mgr->runCmd($cmd);
  
  $mgr->runCmd("gunzip $downloadSubDir/acc2hinv_id.txt.gz");

  $mgr->endStep($signal);
}

sub parseHinvitationalFile {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseHinvitational";

  return if $mgr->startStep("Parse H-inviational file", $signal, 'downloadHinvitational');

  my $taxonId = $propertySet->getProp('taxonId');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');

  my $logfile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $inputFile = "$externalDbDir/h-invitational/$date/acc2hinv_id.txt";

  my $outputFile = "$externalDbDir/h-invitational/$date/acc2hinv_id.txt.out";

  my $cmd = "makeHinvToNaSeqId --taxon_id $taxonId --inputFile  $inputFile > $outputFile 2>> $logfile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub deleteHinvitational2NaSeqId {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteHinvitational2NaSeqId";

  return if $mgr->startStep("Delete H-inviational to na_seq_id", $signal, 'downloadHinvitational');

  my $db_id = $propertySet->getProp('h-inv_db_id');

  my $db_rel_id = $propertySet->getProp('h-inv_db_rls_id');

  my $taxonId = $propertySet->getProp('taxonId');

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n, dots.assembly a where d.external_database_release_id = $db_rel_id and d.db_ref_id = n.db_ref_id and n.na_sequence_id = a.na_sequence_id and a.taxon_id = $taxonId";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub loadHinvitational2NaSeqId {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $taxonId = $propertySet->getProp('taxonId');

  my $db_id = $propertySet->getProp('h-inv_db_id');

  my $db_rel_id = $propertySet->getProp('h-inv_db_rls_id');

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');

  my $inputFile = "$externalDbDir/h-invitational/$date/acc2hinv_id.txt.out";

  my $args = "--mappingfiles $inputFile --pattern1 '(HIT\S+)' --pattern2 'HIT\S+\t(\d+)' --db_id $db_id --db_rel_id $db_rel_id";

  $mgr->runPlugin("loadHinvitational2NaSeqId", "GUS::Common::Plugin::InsertDbRefAndDbRefNaSequenceGeneral", $args, "Load H-inviational file to na_seq_id", 'downloadHinvitational');

}

sub prepareDownloadSiteFiles {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "prepareDownloadSiteFiles";
  
  return if $mgr->startStep("Preparing files for download site", $signal);

  my @files;
  my $htaccessString;
  my $descrip;

  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $speciesFullname = $propertySet->getProp('speciesFullname');
  my $gusConfigFile = $propertySet->getProp('gusConfigFile');
  my $taxonId = $propertySet->getProp('taxonId');
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  # readme file header
  &writeReadmeFileHeader($mgr);

  # fasta file of sequences
  my $fastaFile = "${prefix}.fasta";
  my $logFile = "$mgr->{pipelineDir}/logs/${fastaFile}Extract.log";
  my $sql = "select 'DT.'||na_sequence_id,'[$speciesNickname]',description,'('||number_of_contained_sequences||' sequences)','length='||length,sequence from dots.Assembly where taxon_id = $taxonId";
  my $cmd = "dumpSequencesFromTable.pl --outputFile $mgr->{pipelineDir}/seqfiles/$fastaFile --verbose --gusConfigFile $gusConfigFile  --idSQL \"$sql\" 2>>  $logFile";
  $mgr->runCmd($cmd);
  my $cmd = "gzip $mgr->{pipelineDir}/seqfiles/$fastaFile";
  $mgr->runCmd($cmd) unless -e "$mgr->{pipelineDir}/seqfiles/${fastaFile}.gz";
  push(@files, "$mgr->{pipelineDir}/seqfiles/${fastaFile}.gz");
  $descrip = "The sequence for the consensus transcripts";
  $htaccessString .= "AddDescription \"$descrip\" *fasta*\n";
  addFileToReadme($mgr, "${fastaFile}.gz", $descrip);

  # predicted proteins
  my $predictedProteinsFile = "${prefix}_predictedProteins.fasta";
  my $logFile = "$mgr->{pipelineDir}/logs/${predictedProteinsFile}Extract.log";
  my $sql = "select 'DT.'||a.na_sequence_id,'length of predicted protein sequence ='||ts.length,ts.sequence from dots.translatedaasequence ts,dots.translatedaafeature tf,dots.rnafeature rf,dots.assembly a where a.taxon_id = $taxonId and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_sequence_id = ts.aa_sequence_id";
  my $cmd = "dumpSequencesFromTable.pl --outputFile $mgr->{pipelineDir}/seqfiles/$predictedProteinsFile --verbose --gusConfigFile $gusConfigFile  --idSQL \"$sql\" 2>>  $logFile";
  $mgr->runCmd($cmd) unless -e "$mgr->{pipelineDir}/seqfiles/$predictedProteinsFile" || "$mgr->{pipelineDir}/seqfiles/${predictedProteinsFile}.gz";
  my $cmd = "gzip $mgr->{pipelineDir}/seqfiles/$predictedProteinsFile";
  $mgr->runCmd($cmd) unless -e "$mgr->{pipelineDir}/seqfiles/${predictedProteinsFile}.gz";
  push(@files, "$mgr->{pipelineDir}/seqfiles/${predictedProteinsFile}.gz");
  $descrip = "The predicted protein translation of each assembled transcript";
  $htaccessString .= "AddDescription \"$descrip\" *Proteins*\n";
  addFileToReadme($mgr, "${predictedProteinsFile}.gz", $descrip);

  #DTAnatomy
  my $DTAnatomyFile = "${prefix}_DTAnatomy";
  my $cmd = "makeDTAnatomyFile --taxonId $taxonId --outFile $mgr->{pipelineDir}/misc/$DTAnatomyFile";
  $mgr->runCmd($cmd);
  $mgr->runCmd("gzip $mgr->{pipelineDir}/misc/$DTAnatomyFile");
  push(@files, "$mgr->{pipelineDir}/misc/${DTAnatomyFile}.gz");
  $descrip = "The anatomy percent for assembled transcripts";
  $htaccessString .= "AddDescription \"$descrip\" *Anatomy*\n";
  addFileToReadme($mgr, "${DTAnatomyFile}.gz", $descrip);

  # NRDB Hits 
  my $nrdbHitsFile = "${prefix}_bestNRDBHits.dat";
  my $cmd = "gzip $mgr->{pipelineDir}/misc/$nrdbHitsFile" unless -e "$mgr->{pipelineDir}/misc/${nrdbHitsFile}.gz";
  $mgr->runCmd($cmd);
  push(@files, "$mgr->{pipelineDir}/misc/${nrdbHitsFile}.gz");
  $descrip = "The best hit in NRDB for each consensus transcript";
  $htaccessString .= "AddDescription \"$descrip\" *NRDB*\n";
  addFileToReadme($mgr, "${nrdbHitsFile}.gz", $descrip);

  # Accs Per Assembly
  push(@files, "$mgr->{pipelineDir}/misc/${prefix}_accessionsPerAssembly.dat.gz");
  $descrip = "The Genbank accessions of ESTs and mRNAs contained in each assembled transcript";
  $htaccessString .= "AddDescription \"$descrip\" *_acc*\n";
  addFileToReadme($mgr, "${prefix}_accessionsPerAssembly.dat.gz", $descrip);

  # DTs per DG
  push(@files, "$mgr->{pipelineDir}/misc/${prefix}_DTperDG.dat.gz");
  $descrip = "Assembled transcripts belonging to each gene";
  $htaccessString .= "AddDescription \"$descrip\" *DTperDG*\n";
  addFileToReadme($mgr, "${prefix}_DTperDG.dat.gz", $descrip);

  # mRNAs per DT
  push(@files, "$mgr->{pipelineDir}/misc/${prefix}_mRNAaccessionsPerAssembly.dat.gz");
  $descrip = "The Genbank accessions of mRNAs contained in each assembled transcript";
  $htaccessString .= "AddDescription \"$descrip\" *mRNA*\n";
  addFileToReadme($mgr, "${prefix}_mRNAaccessionsPerAssembly.dat.gz", $descrip);

  # Brain Anatomy terms
  my $brainFile = "${prefix}_brainTerms.dat";
  my $logFile = "$mgr->{pipelineDir}/logs/brainTerms.log";
  my $cmd = "makeAnatomyCountFile --taxonId $taxonId --root brain --rootLevel level_4 --estCount 2 --percent 10 --outputfile $mgr->{pipelineDir}/misc/$brainFile 2>> $logFile";
  $mgr->runCmd($cmd);
  push(@files, "$mgr->{pipelineDir}/misc/$brainFile");
  $descrip = "Brain anatomy terms for which there are DoTS Transcripts.  (Tab delimited: term, term ID, count of DTs)";
  $htaccessString .= "AddDescription \"$descrip\" *brainTerms*\n";
  addFileToReadme($mgr, $brainFile, $descrip);

  # move files to download directory
  foreach my $file (@files) {
    my $outputfile = "$mgr->{pipelineDir}/downloadSite/" . &basename($file);
    my $cmd = "cp $file $outputfile";
    $mgr->runCmd($cmd);
  }

  &updateHtaccessFile($mgr,$htaccessString);

  $mgr->runCmd("chmod g+w $mgr->{pipelineDir}/downloadSite");

  $mgr->endStep($signal);
}

sub writeReadmeFileHeader {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $speciesFullname = $propertySet->getProp('speciesFullname');
 
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $genbankRel = $propertySet->getProp('genbankRel');
  my $nrdbDate = $propertySet->getProp('nrdbDate');
  my $cddFileDates = $propertySet->getProp('cddFileDates');
  my $dbestDate = $propertySet->getProp('dbestDate');
  my $prodomVersion = $propertySet->getProp('prodomRelease');
  my $genomeVersion = $propertySet->getProp('genomeVersion');

  my $readme = "$mgr->{pipelineDir}/downloadSite/${speciesNickname}_README.txt";
  open(F, ">$readme") || $mgr->error("Can't open $readme for writing");

  my $date = `date`;
  chomp $date;
  
  print F "
Release ${dotsRelease} of the Database of Transcribed Sequences (DoTS) for $speciesFullname 
was completed on $date.

The data sources include:
    * GenBank (Release $genbankRel)
    * NRDB ($nrdbDate)
    * dbEST ($dbestDate)
    * Pfam  ($cddFileDates)
    * ProDom ($prodomVersion)
    * CDD ($cddFileDates)
    * Gene Ontology (GO) consortium ontologies and assignments
    * NCBI gene trap tag records from 8 original sources: GGTC, Baygenomics,
      SIGTR, MFGC, CMHD, Lexicon Genetics, and the H.E.Ruley and P.Soriano labs
    * UCSC Genome Bioinformatics Group (genome version $genomeVersion)


The files available are:

";
  close(F);
}

sub addFileToReadme {
  my ($mgr, $filename, $descrip) = @_;
  my $propertySet = $mgr->{propertySet};

  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $readme = "$mgr->{pipelineDir}/downloadSite/${speciesNickname}_README.txt";
  open(F, ">>$readme") || $mgr->error("Can't open $readme for writing");

  print F "
$filename
    $descrip
";
}

sub updateHtaccessFile {
  my ($mgr,$htaccessString) = @_;
  my $propertySet = $mgr->{propertySet};

  my $htaccess = "$mgr->{pipelineDir}/downloadSite/.htaccess";
  
  open(F, ">>$htaccess") || $mgr->error("Can't open $htaccess for writing");
  my $cmd = print F ("$htaccessString");
  close(F);
}

sub downloadGeneId {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "downloadGeneId";
  return if $mgr->startStep("Downloading GeneID", $signal,'downloadGeneId');
  my $pipelineDir = $mgr->{pipelineDir};
  my $logfile = "$pipelineDir/logs/downloadGeneId.log";
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');
  my $downloadSubDir = "$externalDbDir/geneId/$date";
  $mgr->runCmd("mkdir -p $downloadSubDir");
  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=2 -A \"gene2accession.gz,gene_info.gz\"  -P $downloadSubDir  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/";
  $mgr->runCmd($cmd);
  $mgr->runCmd("gunzip $downloadSubDir/gene2accession.gz");
  $mgr->runCmd("gunzip $downloadSubDir/gene_info.gz");
  $mgr->endStep($signal);
}

sub deleteGeneIdToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "deleteGeneIdToNaSeq";
  return if $mgr->startStep("Deleting GeneId to NASeq entries from DbRefNASequence", $signal,'downloadGeneId');
  my $pipelineDir = $mgr->{pipelineDir};
  my $taxonId = $propertySet->getProp('taxonId');
  my $gusConfigFile = $propertySet->getProp('gusConfigFile');
  my $geneIdDbRlsId = $propertySet->getProp('gene_db_rls_id');
  my $logFile = "$pipelineDir/logs/${signal}.log";
  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n, dots.externalnasequence a where d.external_database_release_id = $geneIdDbRlsId and d.db_ref_id = n.db_ref_id and n.na_sequence_id = a.na_sequence_id and a.taxon_id = $taxonId";
  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";
  $mgr->runCmd($cmd);
  $mgr->endStep($signal);
}

sub parseGeneId {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "parseGeneId";
  return if $mgr->startStep("Parsing Gene to na_seq_id via Genbank accesion from file", $signal,'downloadGeneId');
  my $pipelineDir = $mgr->{pipelineDir};
  my $logFile = "$pipelineDir/logs/${signal}.log";
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');
  my $inputFile = "$externalDbDir/geneId/$date/gene2accession";
  my $outputFile = "$pipelineDir/misc/gene2naseq";
  my $tax_id = $propertySet->getProp('ncbiTaxId');
  my $taxonId = $propertySet->getProp('taxonId'); 
  my $cmd = "makeGeneIdToNaSeqFile --inputFile $inputFile --taxon_id $taxonId --tax_id $tax_id > $outputFile 2>>$logFile";
  $mgr->runCmd($cmd);
  $mgr->endStep($signal);
}

sub loadGeneIdToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $pipelineDir = $mgr->{pipelineDir};
  my $file = "$pipelineDir/misc/gene2acc";
  my $geneIdDbRlsId = $propertySet->getProp('gene_db_rls_id');
  my $geneIdDbId = $propertySet->getProp('gene_db_id');
  my $geneIdDeleteDbRef = $propertySet->getProp('deleteGeneId');
  my $delete = $geneIdDeleteDbRef eq "yes" ? "--delete" : "";
  my $args = "--mappingfiles $file $delete --pattern1 '^(\\d+)\\t' --pattern2 '\\t(\\d+)' --db_id $geneIdDbId --db_rel_id $geneIdDbRlsId";
  $mgr->runPlugin("loadGeneIdToNaSeq", "GUS::Common::Plugin::InsertDbRefAndDbRefNASequence", $args, "loading GeneId to NaSeq mapping");
}

sub loadGeneIdInfo {
  my ($mgr)= @_;
  my $propertySet = $mgr->{propertySet};
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');
  my $infoFile = "$externalDbDir/geneId/$date/gene_info";
  my $externalDbRel = $propertySet->getProp('gene_db_rls_id');
  my $tax_id = $propertySet->getProp('ncbiTaxId');
  my $args = "--infoFile $infoFile --ncbiTaxId $tax_id --externalDbRel $externalDbRel";
  $mgr->runPlugin("loadGeneIdInfo", "DoTS::DotsBuild::Plugin::LoadGeneIdInfo", $args, "loading GeneId information");
}

sub downloadMGC {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "downloadMGC";
  return if $mgr->startStep("Downloading MGC", $signal);
  my $pipelineDir = $mgr->{pipelineDir};
  my $logfile = "$pipelineDir/logs/downloadMGC.log";
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  my $date = $propertySet->getProp('buildDate');
  my $downloadSubDir = "$externalDbDir/mgc/$date";
  $mgr->runCmd("mkdir -p $downloadSubDir");
  my $cmd = "wget -t5 -o mgc.log -m -np -nH --cut-dirs=1 -P  $downloadSubDir \"http://mgc.nci.nih.gov/Reagents/StaticCloneList?PAGE=0&STATUS=Confirmed&ORG=Mm\"";
  $mgr->runCmd($cmd);
  $mgr->runCmd("mv \"$downloadSubDir/StaticCloneList?PAGE=0&STATUS=Confirmed&ORG=Mm\" $downloadSubDir/StaticCloneList");
  $mgr->endStep($signal);
}

sub deleteMGC {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteMGCToNaSeq";

  return if $mgr->startStep("Deleting MGC to na_sequence_id entries from DbRefNASequence", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $MgcDbRlsId = $propertySet->getProp('mgc_db_rel_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n, dots.externalnasequence x where d.external_database_release_id = $MgcDbRlsId and d.db_ref_id = n.db_ref_id and n.na_sequence_id = x.na_sequence_id and x.taxon_id = $taxonId";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub parseMGC {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseMGC";

  return if $mgr->startStep("Parsing MGC to DoTS file", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $outFile = "$pipelineDir/misc/mgc2naseq";

  my $taxonId = $propertySet->getProp('taxonId');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/mgc/$date";

  my $inFile = "${downloadSubDir}/StaticCloneList";

  my $cmd = "makeMGC2NaSeqId --taxon_id $taxonId --inputFile $inFile --verbose > $outFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub loadMGCToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $pipelineDir = $mgr->{pipelineDir};

  my $file = "$pipelineDir/misc/mgc2naseq";

  my $db_rel_id = $propertySet->getProp('mgc_db_rel_id');

  my $db_id = $propertySet->getProp('mgc_db_id');

  my $args = "--mappingfiles $file --delete --pattern1 '(MGC:\\d+)' --pattern2 '\\t(\\d+)' --db_id $db_id --db_rel_id $db_rel_id";
  
  $mgr->runPlugin("loadMGCMapping", "GUS::Common::Plugin::InsertDbRefAndDbRefNaSequenceGeneral", $args, "loading MGC to DoTS mapping");
}

sub loadMGCInfo {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $pipelineDir = $mgr->{pipelineDir};
  
  my $logfile = "$pipelineDir/logs/loadMGCInfo.log";

  my $db_rel_id = $propertySet->getProp('mgc_db_rel_id');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/mgc/$date";

  my $file = "$downloadSubDir/StaticCloneList";

  my $args = "--infoFile $file --externalDbRel $db_rel_id";

  $mgr->runPlugin("loadMGCInfo", "DoTS::DotsBuild::Plugin::LoadMGCInfo", $args, "loading MGC Info from file");
}

sub deleteFantom {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteFantomToNaSeq";

  return if $mgr->startStep("Deleting Fantom to na_sequence_id entries from DbRefNASequence", $signal,'loadFantom');

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $FantomDbRlsId = $propertySet->getProp('fantom_db_rel_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n, dots.externalnasequence x where d.external_database_release_id = $FantomDbRlsId and d.db_ref_id = n.db_ref_id and n.na_sequence_id = x.na_sequence_id and x.taxon_id = $taxonId";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub parseFantom {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseFantom";

  return if $mgr->startStep("Parsing Fantom to na_sequence_id file", $signal,'loadFantom');

  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $outFile = "$pipelineDir/misc/fantom2naseq";

  my $taxonId = $propertySet->getProp('taxonId');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/fantom2/$date";

  my $inFile = "${downloadSubDir}/FANTOM2setInfo.txt";

  my $cmd = "makeFantom2ToNaSeqId --taxon_id $taxonId --inputFile $inFile --verbose > $outFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub loadFantomToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $pipelineDir = $mgr->{pipelineDir};

  my $file = "$pipelineDir/misc/fantom2naseq";

  my $db_rel_id = $propertySet->getProp('fantom_db_rel_id');

  my $db_id = $propertySet->getProp('fantom_db_id');

  my $args = "--mappingfiles $file --delete --pattern1 '^(\\S+)' --pattern2 '\\t(\\d+)' --db_id $db_id --db_rel_id $db_rel_id";

  $mgr->runPlugin("loadFantomMapping", "GUS::Common::Plugin::InsertDbRefAndDbRefNaSequenceGeneral", $args, "loading Fantom to na_sequence_id mapping", 'loadFantom');
}


sub downloadMGIInfo {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "downloadMGIInfo";

  return if $mgr->startStep("Downloading MGI Info", $signal,'loadMGI');

  my $pipelineDir = $mgr->{pipelineDir};

  my $logfile = "$pipelineDir/logs/downloadMGIInfo.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/mgi/$date";

  $mgr->runCmd("mkdir -p $downloadSubDir");
    
  my $cmd = "wget -t5 -o $logfile  -m -np -nd -nH --cut-dirs=2 -A \"MRK_Dump2.rpt,MRK_Sequence.rpt,MGI_EntrezGene.rpt\"  -P $downloadSubDir  ftp://ftp.informatics.jax.org/pub/reports/";

  $mgr->runCmd($cmd);

  my $cmd = "wget -t5 -o $logfile -m -np -nd -nH --cut-dirs=3 -A \"MGI_DT_via_GB_one_*\"  -P $downloadSubDir  ftp://ftp.informatics.jax.org/pub/reports/dotstigr/";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub deleteMGIToDots {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "deleteMGIToDots";

  return if $mgr->startStep("Deleting MGI to DoTS entries from DbRefNASequence", $signal,'loadMGI');
  my $pipelineDir = $mgr->{pipelineDir};

  my $mgiDbRlsId = $propertySet->getProp('mgi_db_rls_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n where d.external_database_release_id = $mgiDbRlsId and d.db_ref_id = n.db_ref_id";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub parseMgiToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseMGIToNaSeq";

  return if $mgr->startStep("Parsing MGI to na_sequence_ids from file", $signal, 'loadMGI');

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $taxonId = $propertySet->getProp('taxonId');

  my $inputFile = "$externalDbDir/mgi/$date/MRK_Sequence.rpt";
  my $pipelineDir = $mgr->{pipelineDir};
  my $logFile = "$pipelineDir/logs/${signal}.log";
  my $outputFile = "$pipelineDir/misc/mgi2naseq";

  my $cmd = "makeMgiToNaSeqFile --taxon_id $taxonId --fileFile $inputFile > $outputFile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub loadMgiToNaSeq {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $pipelineDir = $mgr->{pipelineDir};

  my $file = "$pipelineDir/misc/mgi2naseq";
  my $mgiDbRlsId = $propertySet->getProp('mgi_db_rls_id');
  my $mgiDbId = $propertySet->getProp('mgi_db_id');


  my $args = "--mappingfiles $file --pattern1 '(MGI:\\d+)\\t' --pattern2 '\\t(\\d+)' --db_id $mgiDbId --db_rel_id $mgiDbRlsId";
  $mgr->runPlugin("loadMgiToNaSeq", "GUS::Common::Plugin::InsertDbRefAndDbRefNASequence", $args, "loading MGI to NaSeq mapping", 'loadMGI');
}



sub loadMGIToDoTS {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $externalDbDir = $propertySet->getProp('externalDbDir');
  
  my $date = $propertySet->getProp('buildDate');
  
  my $downloadSubDir = "$externalDbDir/mgi/$date";

  my @file = split (/,/, $propertySet->getProp('mgiFiles'));
  
  my $files;
  
  foreach my $f (@file) {
    $files .= "$downloadSubDir/$f,";
  }

  chop($files);

  my $db_rel_id = $propertySet->getProp('mgi_db_rls_id');
  my $mgiDbId = $propertySet->getProp('mgi_db_id');

  my $args = "--mappingfiles $files --delete --pattern1 '(MGI:\\d+)' --pattern2 'DT\.(\\d+)' --db_id $db_rel_id  --db_rel_id $db_rel_id";
  
  $mgr->runPlugin("loadMGIToDoTS", "GUS::Common::Plugin::InsertDbRefAndDbRefNASequence", $args, "loading MGI to DoTS mapping",'loadMGI');
}

sub loadMGIInfo {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/mgi/$date";

  my $infoFile = "$downloadSubDir/MRK_Dump2.rpt";

  my $geneFile =  "$downloadSubDir/MGI_EntrezGene.rpt";

  my $db_rel_id = $propertySet->getProp('mgi_db_rls_id');
  
  my $args = "--infoFile $infoFile --geneFile $geneFile --external_db_release_id $db_rel_id";
  
  $mgr->runPlugin("loadMGIInfo", "DoTS::DotsBuild::Plugin::LoadMGIInfo", $args, "Loading MGI Info",'loadMGI' );
  
}

sub deleteGeneCardsToDots {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "deleteGeneCardsToDots";

  return if $mgr->startStep("Deleting GeneCards to DoTS entries from DbRefNASequence", $signal,'loadGeneCards');
  my $pipelineDir = $mgr->{pipelineDir};

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $geneCardsDbRlsId = $propertySet->getProp('genecards_db_rls_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n where d.external_database_release_id =$geneCardsDbRlsId and d.db_ref_id = n.db_ref_id";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub parseGeneCardsToDoTS {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "parseGeneCardsToDots";

  return if $mgr->startStep("Parsing GeneCards to DoTS entries from file", $signal,'loadGeneCards');
  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $date = $propertySet->getProp('buildDate');

  my $downloadSubDir = "$externalDbDir/genecards/$date";

  my $inputfile = "$downloadSubDir/dumpForDots.txt";

  my $outputfile = "$downloadSubDir/geneCards2DoTS.txt";

  my $cmd = "DoTS2GeneCardsParse --inputFile $inputfile --outputFile $outputfile 2>>$logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}
    
sub loadGeneCardsToDoTS {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet}; 

  my $externalDbDir = $propertySet->getProp('externalDbDir');
  
  my $date = $propertySet->getProp('buildDate');
  
  my $downloadSubDir = "$externalDbDir/genecards/$date";

  my $file = "$downloadSubDir/geneCards2DoTS.txt";

  my $geneCardsDbRlsId = $propertySet->getProp('genecards_db_rls_id');

  my $args = "--mappingfiles $file --delete --pattern1 '(\\S+)\\t' --pattern2 '\\tDT\.(\\d+)' --db_id 195 --db_rel_id $geneCardsDbRlsId";
  
  $mgr->runPlugin("loadGeneCardsMapping", "GUS::Common::Plugin::InsertDbRefAndDbRefNASequence", $args, "loading GeneCards to DoTS mapping",'loadGeneCards');
}

sub deleteGEAToDoTS {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "deleteGEAToDoTS";

  return if $mgr->startStep("Deleting GEA to DoTS entries from DbRefNASequence", $signal);

  my $gusConfigFile = $propertySet->getProp('gusConfigFile');

  my $pipelineDir = $mgr->{pipelineDir}; 

  my $geaDbRlsId = $propertySet->getProp('gea_db_rls_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $sql = "select n.db_ref_na_sequence_id from sres.dbref d, dots.dbrefnasequence n where d.external_database_release_id in ($geaDbRlsId) and d.db_ref_id = n.db_ref_id";

  my $cmd = "deleteEntries.pl --table DoTS::DbRefNASequence --idSQL \"$sql\" --gusConfigFile $gusConfigFile --verbose 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub parseGEA {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "parseGEAToDots";

  return if $mgr->startStep("Parsing GEA to DoTS entries from file", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $downloadSubDir = "$externalDbDir/gea";

  my $taxonId = $propertySet->getProp('taxonId');

  my @geaFiles = split(/,/, $propertySet->getProp('geaFiles'));
  
  foreach my $file (@geaFiles) {

    my ($geaFile, $db_id, $db_rel_id, $regex) = split(/:/, $file);

    my $inputfile = "$downloadSubDir/$geaFile";
  
    my $outputfile = "$pipelineDir/misc/${geaFile}2DoTS.txt";

    my $cmd = "makeGEA2DoTSfile --inputFile $inputfile --taxon_id $taxonId --regex \"$regex\" > $outputfile 2>>$logFile";

    $mgr->runCmd($cmd);

  }


  $mgr->endStep($signal);
}
    
sub loadGEA {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "loadGEA";

  return if $mgr->startStep("Loading GEA to na_sequence_id entries from file", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $externalDbDir = $propertySet->getProp('externalDbDir');

  my $downloadSubDir = "$externalDbDir/gea";

  my $species = $propertySet->getProp('speciesNickname');

  my @geaFiles = split(/,/, $propertySet->getProp('geaFiles'));

  foreach my $file (@geaFiles) {

    my ($geaFile, $db_id , $relId, $regex) = split(/:/, $file);

    my $file = "$pipelineDir/misc/${geaFile}2DoTS.txt";

    my $args = "--mappingfiles $file --pattern1 '(\\S+)\\t' --pattern2 '\\tDT\.(\\d+)' --db_id $db_id --db_rel_id $relId";

    $mgr->runPlugin("loadGEA${geaFile}Mapping", "GUS::Common::Plugin::InsertDbRefAndDbRefNaSequenceGeneral", $args, "loading GEA to DoTS mapping");
  }
  $mgr->endStep($signal);
}



sub makeProjectLink {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "makeProjectLink";

  return if $mgr->startStep("Insert links between projectinfo and nasequence into projectlink table", $signal);

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesVer = $propertySet->getProp('allgenesVersion');

  my $project_id = $propertySet->getProp('project_id');

  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";

  my $args;

  my $project = $project_id > 0 ? "--project_id ".$project_id : "";
  $args = "--commit --verbose --allgenes_num $allgenesVer --taxon $taxonId $project";

  $args .= " --restart"
    if ($propertySet->getProp('projectLinkRestart') eq "yes");

  my $cmd = "makeProjectLink.pl $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);
}

sub refreshMaterializedViews {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $signal = "refreshMaterializedViews";
  my $logFile = "$mgr->{pipelineDir}/logs/${signal}.log";
  return if $mgr->startStep("Refresh MaterializedViews", $signal);
  my $materializedViews = $propertySet->getProp('materializedViews');
  my @views = split (/\,/,$materializedViews);
  foreach my $schemaView (@views) {
    my ($schema,$view) = split (/\:/,$schemaView);
    my $args = "--materializedView $view --schema $schema --verbose";
    my $cmd = "refreshMaterializedView $args 2>> $logFile";
    $mgr->runCmd($cmd);
  }
  $mgr->endStep($signal);
}

sub updateProteinAssemblyTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateProteinAssemblyTable";

  return if $mgr->startStep("Insert entries in allgenes ProteinAssembly table ", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updateProteinAssemblyTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);


}

sub updateCentralDogmaTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateCentralDogmaTable";

  return if $mgr->startStep("Insert entries in allgenes CentralDogma table ", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updateCentralDogmaTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);


}

sub updateDTOrfPValueTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateDTOrfPValueTable";

  return if $mgr->startStep("Insert entries in allgenes makeDTOrfPValue table", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updateDTOrfPvalueTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub updatePromoterRegionTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updatePromoterRegionTable";

  return if $mgr->startStep("Insert entries into Allgenes.PromoterRegion table", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $genome_rls_id = $propertySet->getProp('genome_db_rls_id');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $args = "--taxon $taxonId --genomeExtRelDbId $genome_rls_id --allgenesSchema $allgenesSchema";

  my $cmd = "updatePromoterRegionTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub updatePancreasAssembliesTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updatePancreasAssembliesTable";

  return if $mgr->startStep("Insert entries in Allgenes.PancreasAssemblies table", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updatePancreasAssembliesTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub updateAssemblySignalPSummaryTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateAssemblySignalPSummaryTable";

  return if $mgr->startStep("Insert entries in allgenes makeAssemblySignalPSummary table", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updateAssemblySignalPSummaryTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub updateAssemblyTMDomainSummaryTable {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "updateAssemblyTMDomainSummaryTable";

  return if $mgr->startStep("Insert entries in allgenes makeAssemblyTMDomainSummary table", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $taxonId = $propertySet->getProp('taxonId');

  my $allgenesSchema = $propertySet->getProp('allgenesSchema');

  my $logFile = "$pipelineDir/logs/${signal}.log";

  my $args = "--verbose --taxon $taxonId --allgenesSchema $allgenesSchema";

  my $cmd = "updateAssemblyTMDomainSummaryTable $args 2>> $logFile";

  $mgr->runCmd($cmd);

  $mgr->endStep($signal);

}

sub createPredTranslDetailsFile {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "createPredTranslDetailsFile";

  return if $mgr->startStep("Preparing predicted translation details file", $signal);

  my $pipelineDir = $mgr->{pipelineDir};

  my $dotsRelease = $propertySet->getProp('dotsRelease');

  my $speciesNickname = $propertySet->getProp('speciesNickname');

  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";

  my $predTranslDetailsFile = "$pipelineDir/downloadSite/${prefix}_predictedProteinDetails.txt";

  my $taxonId = $propertySet->getProp('taxonId');
  
  my $cmd = "getAssTranslATGStartStop --taxon_id $taxonId >$predTranslDetailsFile";  
  $mgr->runCmd($cmd);

  my $cmd = "gzip $predTranslDetailsFile";
  $mgr->runCmd($cmd);

  my $descrip = "A tab delimited file of the translation details for proteins predicted using FrameFinder";
  &updateHtaccessFile($mgr, "AddDescription \"$descrip\" *Details*\n");
  &addFileToReadme($mgr, "${prefix}_predictedProteinDetails.txt.gz", $descrip);

  $mgr->endStep($signal);
}

sub createManuallyReviewedDoTSFile { 
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "createManuallyReviewedDoTSFile";

  return if $mgr->startStep("Preparing manually reviewed DoTs file", $signal);

  my $pipelineDir = $mgr->{pipelineDir};
  my $dotsRelease = $propertySet->getProp('dotsRelease');
  my $speciesNickname = $propertySet->getProp('speciesNickname');
  my $tempLogin = $propertySet->getProp('tempLogin');  
  my $tempPassword = $propertySet->getProp('tempPassword'); 
  my $taxonId = $propertySet->getProp('taxonId');
  my $allgenesLogin  = $propertySet->getProp('allgenesLogin');
  
  my $prefix = "${speciesNickname}DoTS_rel${dotsRelease}";
  
  my $manRevwDoTSReportFile = "$pipelineDir/downloadSite/${prefix}_manuallyReviewedTranscriptsReport.txt";
    
  # create report of manually reviewed dots
  my $cmd = "gusreport --configModule GUS::ReportMaker::SampleTranscriptReportConfig --tempTableName tempResult --requestedColumns 'Description, GeneSymbol, Length, GOid, DoTSGene, SeqsInAssem, MGI, LocusLink, ContainsMRNA, Motifs, Organism, mRNASeq, proteinSeq' --sql 'select distinct pa.na_sequence_id, 9999 from DoTS.rna r, allgenes.proteinassembly pa where r.review_status_id = 1 and r.rna_id = pa.rna_id and pa.taxon_id = $taxonId' > $manRevwDoTSReportFile";
  $mgr->runCmd($cmd);

  my $cmd = "gzip $manRevwDoTSReportFile";

  $mgr->runCmd($cmd);

  my $descrip = "A tab delimited report for all manually reviewed DoTS Transcripts";
  &updateHtaccessFile($mgr, "AddDescription \"$descrip\" *manuallyReviewed*\n");
  &addFileToReadme($mgr, "${prefix}_manuallyReviewedTranscriptsReport.txt.gz", $descrip);

  $mgr->endStep($signal);
}

sub makeStatisticsPage {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};

  my $signal = "makeStatisticsPage";

  return if $mgr->startStep("Making statistics html file", $signal,'makeStatPage');

  my $pipelineDir = $mgr->{pipelineDir};

  my $output = "$pipelineDir/misc/statistics.html.in";

  my $logfile = "$pipelineDir/logs/${signal}.log";

  my $cmd = "dotsStatistics --templateFile $ENV{GUS_HOME}/data/DoTS/DotsBuild/statsTemplate.html.in --goRowTemplateFile $ENV{GUS_HOME}/data/DoTS/DotsBuild/statsGoRowTemplate.txt --goBarGraphWidth 200 > $output 2>> $logfile"; 

  $mgr->runCmd($cmd);

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
  my $prog = `basename $0`;
  chomp $prog;
  print STDERR "usage: $prog propertiesfile\n";
  exit 1;
}

############################# Utility Subroutines #############################

# get list of taxonIds for the current species.  If the property
# "includeSubspecies" equals "yes", then that's a comma-separated list
# including the species taxonId (from the property "taxon_id") together with
# any descendent nodes in the SRes::Taxon tree.  If includeSubspecies isn't
# set, it's just the taxon_id (again, from the property).

sub getTaxonIdList {
  my ($mgr) = @_;
  my $propertySet = $mgr->{propertySet};
  my $returnValue;

  my $taxonId = $propertySet->getProp('taxonId');
  if ($propertySet->getProp('includeSubspecies') eq "yes") {
    $returnValue = $mgr->runCmd("getSubTaxa --taxon_id $taxonId");
    chomp $returnValue;
  } else {
    $returnValue = $taxonId;
  }

  return $returnValue;
}
##################################### main ####################################

1;
