package DoTS::DotsBuild::Plugin::SetAssSeqQualStartStop;

@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use DBI;

############################################################
# Add any specific objects (Objects::GUSdev::Objectname) here
############################################################
use GUS::Model::DoTS::AssemblySequence;
use CBIL::Bio::SequenceUtils;

sub new {
  my ($class) = @_;

  my $self = {};
  bless($self,$class); 

  my $usage = 'Set the quality start and stop of existing assembly sequences with ids obtained with user supplied sql';

  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'idSQL',
      t => 'string',
      h => 'SQL statement:  must return list of primary identifiers from AssemblySequence',
     },
     {o => 'setSeqStartEnd',
      t => 'boolean',
      h => 'also set the sequence_start and sequence_end...default is to not do this',
     },
     {o => 'repeatFile',
      t => 'string',
      h => 'repeat file with full path',
     },
     {o => 'phrapDir',
      t => 'string',
      h => 'directory with phrap files',
     },
    ];

  $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$', # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     revisionNotes => 'make consistent with GUS 3.0',
		     easyCspOptions => $easycsp,
		     usage => $usage
                    });
  return $self;
}


my $countProcessed = 0;
my $countBad = 0;
my $library;
my $debug = 0;
$| = 1;

sub run {
  my $self   = shift;

  my $tmpFile = "tmpFile.$$";

  my $dbiDb = $self->getDb();
  $dbiDb->setMaximumNumberOfObjects(100000);
  $dbiDb->setGlobalNoVersion(1);


  print $self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on ". $self->getArgs()->{'testnumber'}."\n" if $self->getArgs()->{'testnumber'};

  die "you must provide --idSQL \n" unless $self->getArgs()->{'idSQL'};
  die "you must provide a repeat file\n" unless $self->getArgs()->{'repeatFile'};
  die "you must provide the phrap dir\n" unless $self->getArgs()->{'phrapDir'};
  $library = $self->getArgs()->{'repeatFile'};

  ##DBI handle to be used for queries outside the objects...
  my $dbh = $self->getQueryHandle();

  my $stmt = $dbh->prepareAndExecute($self->getArgs()->{idSQL});
  my @ids;
  my $ct = 0;
  while(my($id) = $stmt->fetchrow_array()){
    push(@ids,$id);
    $ct++;
    print STDERR "Retrieving $ct\n" if $ct % 10000 == 0;
    last if $self->getArgs()->{testnumber} && $ct >= $self->getArgs()->{testnumber};
  }
  print "making ",scalar(@ids)," assemblysequences quality\n";

  my $miniLib = "";
  my $count = 0;
  foreach my $id (@ids){
    my $aseq = GUS::Model::DoTS::AssemblySequence->new({'assembly_sequence_id'=> $id});
    if(!$aseq->retrieveFromDB()){ print STDERR "unable to retrieve AssemblySequence.$id\n"; next;}
    my $ex = $aseq->getParent('GUS::Model::DoTS::ExternalNASequence',1);
    ##there could be more than one lenssequence child...want newest one...
    my @ls_sort = sort{$b->getId() <=> $a->getId()} $ex->getChildren('GUS::Model::DoTS::EST',1);
    my $ls;
    $ls = $ls_sort[0] if scalar(@ls_sort) > 0;
    my $qStop;
    $qStop = $ls->getQualityStop() if $ls;
    if($ls && defined $qStop){  ##exists and has qualityStop...assume quality_start is not relevant if qualitystop is null..
      my $qStart = $ls->getQualityStart() > 0 ? $ls->getQualityStart() - 1 : 0;
      ##turns out that $qStart may be > length of the sequence...set to the length...
      if($qStart > $ex->getLength() - 1){$qStart =  $ex->getLength() - 1;}
      ##set length to 20 if less than thatt so cross_match will work...will be removed from AssembySequence
      ##since is less than 50 bp.
      my $qLength = $qStop == 0 ? $ex->getLength() - $qStart : ($qStop - $qStart < 20 ? 20 : $qStop - $qStart);
      $miniLib .= ">".$id."\n".CBIL::Bio::SequenceUtils::breakSequence(substr($ex->getSequence(),$qStart,$qLength));
    }else{
      $miniLib .= ">".$id."\n".CBIL::Bio::SequenceUtils::breakSequence($ex->getSequence());
    }
    $count++;
    if ($count >= 10000) {
      $self->processSet($miniLib,$tmpFile);
      $miniLib = "";						##reset for next set of seqs
      $count = 0;
      $ctx->{'self_inv'}->undefPointerCache();
    }
  }
  &processSet($miniLib,$tmpFile);				##processes last set
  $dbiDb->undefPointerCache();

  ##delete thetmp files.
  unlink("$tmpFile");
  unlink("$tmpFile.screen");
  unlink("$tmpFile.log");

  my $ret = "Processed $countProcessed AssemblySequences, marked $countBad as low_quality";
  print "$ret\n";
  return $ret;
  
}


sub processSet {
  my($self,$miniLib,$tmpFile) = @_;

  print STDERR "Starting next set following $countProcessed ",`date`;
  open(S, ">$tmpFile");
  print S $miniLib;
  close S;
  
  ##NOTE that need to get this installed correctly...
  #  system("/usr/local/src/bio/PHRAP_etal/phrap.SUNWspro/ultra.bin/cross_match tmpLib /usr/local/db/others/repeat/vector -screen > cross.test");
  ##NOTE: need to use species specific library if human or mouse as am removing ribosomal and mitochondrial sequences for these as well...
  my $cmd = "$phrap_dir/cross_match $tmpFile $library -screen > cross.test 2> /dev/null";
  print STDERR "$cmd\n" if $debug;
  system($cmd);

  ##generate better sequence....
  open(S,"$tmpFile.screen");
  my $seq;
  my $na_seq_id;
  my @sub;
  while (<S>) {
    if (/^\>(\d+)/) {           ##$1 = na_sequence_id
      my $as = $self->setQualStartStop($na_seq_id,$seq) if $na_seq_id;
      if($as){
        push(@sub,$as);
      }elsif($na_seq_id){ 
        print STDERR "ERROR: Unable to set quality for $na_seq_id\n";
      }
      $na_seq_id = $1;
      $seq = "";
      if($countProcessed % 100 == 0){
        &submitAssSeqs(@sub);
        undef @sub;
        print STDERR "Processed: $countProcessed, low_quality: $countBad\n";
      }
    } else {
      $seq .= $_;
    }
  }
  my $as = $self->setQualStartStop($na_seq_id,$seq) if $na_seq_id;
  if($as){
    push(@sub,$as);
  }else{ 
    print STDERR "ERROR: Unable to set quality for $na_seq_id\n";
  }
  $self->submitAssSeqs(@sub);
  print STDERR "Processed: $countProcessed, low_quality: $countBad\n";
  close S;
}

sub submitAssSeqs {
  my ($self,@seqs) = @_;
  my $dbiDb = $self->getDb();
  $dbiDb->manageTransaction(undef,'begin');
  foreach my $a (@seqs){
    $a->submit(0,1);
  }
  $dbiDb->manageTransaction(undef,'commit');
}

sub setQualStartStop {
  my($self,$id,$seq) = @_;
  my $dbiDb = $self->getDb();

  my $ass = $dbiDb->getFromDbCache('GUS::Model::DoTS::AssemblySequence',$id);

  if(!$ass){ print STDERR "Unable to retrieve AssSeq.$id from dbCache\n"; return undef; }

  my $qseq = $self->returnQuality($seq);

  my($start,$stop) = $self->getStartStop($ass->getParent('GUS::Model::DoTS::ExternalNASequence')->getSequence(),$qseq);
  return undef if ! defined $start;
  $ass->setQualityStart($start);
  $ass->setQualityEnd($stop);
  if($stop - $start + 1 < 50){
    $ass->setProcessedCategory('low_quality') unless $ass->getProcessedCategory() eq 'low_quality';
  }elsif($ass->getProcessedCategory()){
    $ass->setProcessedCategory('NULL');
  }
  print $ass->toXML(0,1) if $debug;
  $countProcessed++;
  return $ass;
}

sub getStartStop {
  my($self,$nas,$sequence) = @_;
  $nas =~ tr/a-z/A-Z/;          #3upper case it as this is what do with assemblySequence..
  $nas =~ s/\s//g;
  $sequence =~ tr/a-z/A-Z/;
  $sequence =~ s/\s//g;
  ##there could be residual N's on the ends from blocking and trimming.....remove just in case....
  $sequence =~ s/^N*(\w*?)N*$/$1/;
  $countBad++ if length($sequence) < 50;
  my $index = index($nas,$sequence);
  if($index == -1){ return undef; } ##didn't match
  return ($index + 1, $index + length($sequence));
}

##trims leadint TTTT and trailing AAAAA and scans for stretch of good sequence
##by looking at the %Ns in a 20 bp window...must be less than 20%
sub returnQuality{
  my($self,$seq) = @_;
  my $newSeq;
  $seq =~ s/\s+//g;             ##gets rid of spaces and newlines
  $seq =~ s/X/N/g;              ##subsititutes the X's put in by cross_match with N's
  my $newstart = 0;
  for (my $i=0; $i<(length($seq) - 20); $i++) {
    my $tmp = substr($seq, $i, 20);
    $tmp =~ s/[ACGTacgt]//g;
    if (length($tmp) >= 5) {
      if (($i - $newstart) < 40) { ##looking for a 40 bp region of good sequence...=50 at end..
        $i = $i + 9;
        $newstart = $i;
				#	print STDERR "Trimming Poor qual beginning: newStart set to $newstart\n";
      } else {
        $newSeq = substr($seq, $newstart, ($i + 10 - $newstart));
        return $self->trimAT($newSeq); 
      }
    }
  }
  $newSeq = substr($seq, $newstart, length($seq));
  return $self->trimAT($newSeq);
}

##steps through sequence with 30 bp window...if >=36 As or Ts then truncates
##either beginning or end (whichever is closer)
sub trimAT {
  my($self,$seq) = @_;
  $seq =~ s/\s//g;
  print STDERR "\ntrimAT input: \n", CBIL::Bio::SequenceUtils::breakSequence($seq) if $debug == 1;
  my @seq = split('', $seq);
  my %nuc;
  my($startBase,$endBase);
  foreach my $n (@seq[0..29]) {
    $nuc{$n}++;
  }
  print STDERR "First 30: ", join('', @seq[0..29]), "\n" if $debug == 1;
  my $inRun = 0;
  for (my $i=1;$i<(length($seq)-30);$i++) {
    if ($inRun == 0 && ($nuc{A} >=27 || $nuc{T} >= 27)) {
      $startBase = $i;
      print STDERR "A=$nuc{A},T=$nuc{T} Sequence Length: ", length($seq), " StartBase: $startBase " if $debug == 1;
      $inRun = 1;
      next;
    }
    if ($inRun == 1 && ($nuc{A} < 27 && $nuc{T} < 27)) {
      ##process the sucker!!
      $endBase = $i + 30;
      print STDERR "EndBase: $endBase\n" if $debug == 1;
      ##test if closer to beginning or end and tructate accordingly
      if ($startBase <= (length($seq) - $endBase)) {
        print STDERR "substr\(\$seq,\($endBase-5\),\(length\(\$seq\)\)\)\n" if $debug == 1;
        return substr($seq,($endBase-5),(length($seq)));
      } else {
        print STDERR "substr\(\$seq,0,\($startBase+5\)\)\n" if $debug == 1;
        return substr($seq,0,($startBase+5));
      }
    }
    ##update %nuc
    $nuc{$seq[$i]}--;           ##removes the first character
    $nuc{$seq[$i+29]}++;        ##adds new nucleotide
  }
  if ($inRun == 1) {            ##stretch goes to end of sequence
    ##truncate from startBase
    print STDERR "substr\(\$seq,0,\($startBase+5\)\)\n" if $debug == 1;
    return substr($seq,0,$startBase);
  }
  return $seq;
}

1;

