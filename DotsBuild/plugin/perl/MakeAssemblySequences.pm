############################################################
## Change Package name....
############################################################
package MakeAssemblySequences;

use strict;

use Objects::GUSdev::ExternalNASequence;
use Objects::GUSdev::AssemblySequence;
use FileHandle;
use SequenceUtils;

##for query handle on database....
use DBI;

sub new {
  my $Class = shift;

  return bless {}, $Class;
}

sub Usage {
  my $M   = shift;
  return 'Retrieves ExternalNASequences that are not assembly sequences and inserts AssemblySequence entry';
}

############################################################
# put the options in this method....
############################################################
sub EasyCspOptions {
  my $M   = shift;
  {

    #		test_opt1 => {
    #									o => 'opt1=s',
    #									h => 'option 1 for test application',
    #									d => 4,
    #									l => 1,	ld => ':',
    #									e => [ qw( 1 2 3 4 ) ],
    #								 },

 testnumber        => {
                       o => 'testnumber=i',
                       h => 'number of iterations for testing',
                      },
                        
 date              => {
                       o => 'date=s',
                       h => 'earliest date for sequences to include: default = all',
                      },

 taxon_id          => {
                       o => 'taxon_id=i',
                       h => 'taxon_id for sequences to process: 8=Hum, 14=Mus.',
                      },
 idSQL             => {
                       o => 'idSQL',
                       t => 'string',
                       h => 'SQL query that returns na_sequence_ids from ExternalNASequence to be processed',
                      },
 export            => {
                       o => 'export',
                       t => 'string',
                       h => 'filename to export sequences to...default does not export',
                      },
                        
                    }
}

my $debug = 0;
my $ctx;
my $countProcessed = 0;
my $countBad = 0;
my $dbh;
my $phrap_dir = '/usr/local/src/bio/phrap/latest';
my %finished;
my $library;

$| = 1;

sub Run {
  my $M   = shift;
  $ctx = shift;

  print $ctx->{'cla'}->{'commit'} ? "COMMIT ON\n" : "COMMIT TURNED OFF\n";
  print "Testing on $ctx->{'cla'}->{'testnumber'}\n" if $ctx->{'cla'}->{'testnumber'};

  ##set the taxon_id...
  die "You must enter either the --taxon_id and optionally --idSQL on the command line\n" unless $ctx->{cla}->{taxon_id};
    ##set up the library:
  if($ctx->{cla}->{taxon_id} == 8){
    $library = '/usr/local/db/others/repeat/vector_humMitoRibo.lib';
  }elsif($ctx->{cla}->{taxon_id} == 14){
    $library = '/usr/local/db/others/repeat/vector_musMitoRibo.lib';
  }else{
    $library = '/usr/local/db/others/repeat/vector';
    print STDERR "No taxon_id specified....";
  }
  print STDERR "Running cross_match with $library\n";

  $ctx->{'self_inv'}->setMaximumNumberOfObjects(100000);

  $dbh = $ctx->{self_inv}->getQueryHandle();

  if ($ctx->{cla}->{export}) {
    if (-e $ctx->{cla}->{export}) {
      open(F,"$ctx->{cla}->{export}");
      while (<F>) {
        if (/^\>(\S+)/) {
          $finished{$1} = 1;
        }
      }
      close F;
      print STDERR "Already processed ",scalar(keys%finished)," ids\n";
    }
    open(EX,">>$ctx->{cla}->{export}");
  }

  if ($ctx->{cla}->{idSQL}) {
    &processQuery($ctx->{cla}->{idSQL});
  } else {

    print "Generating new AssemblySequences for taxon(s) $ctx->{cla}->{taxon_id}\n";

    ##first get the ESTs and mRNAs..
    my $sql = "select na_sequence_id from ExternalNASequence where taxon_id in ($ctx->{cla}->{taxon_id}) ".
      "and sequence_type_id in (7,8) " .
        "and na_sequence_id not in (select na_sequence_id from AssemblySequence) ";
    
    if ($ctx->{'cla'}->{'date'}) {
      $sql .= "and modification_date >= '$ctx->{'cla'}->{'date'}'";
    } 
    
    &processQuery($sql);
    
    ##next get the things from embl that are RNAs longer than 500 bp...
    ##need to check this for things that are not human or mouse...may need to use less sophisticated query!
    my $mRNASql = "select na_sequence_id from externalnasequence o where o.na_sequence_id in (
select s.na_sequence_id from externalnasequence s, transcript t, nalocation l
where s.taxon_id = $ctx->{cla}->{taxon_id} 
and s.sequence_type_id = 2 
and t.na_sequence_id = s.na_sequence_id
and t.name = 'CDS'
and l.na_feature_id = t.na_feature_id
group by s.na_sequence_id having count(*) = 1 )
and o.length > 400
and o.na_sequence_id not in (select na_sequence_id from AssemblySequence)";

    #more general query that may not be as good...
    #    my $mRNASql = "select na_sequence_id from ExternalNASequence
    # where sequence_type_id = 2
    # and external_db_id = 78
    # and taxon_id in ( $ctx->{'cla'}->{'taxon_id'} )
    # and length > 500
    # and na_sequence_id not in (select na_sequence_id from AssemblySequence)";
    
    if ($ctx->{'cla'}->{'date'}) {
      $mRNASql .= "and modification_date >= '$ctx->{'cla'}->{'date'}'";
    } 
    &processQuery($mRNASql);

  }
    
  # unlink "tmpLib";
  # unlink "tmpLib.mask";

  close EX;
  my $results = "Processed $countProcessed AssemblySequences and marked $countBad as 'low_quality'";

  print "\n$results\n";

  return $results;
}

sub processQuery {
  my($sql) = @_;

  print STDERR "\n$sql\n"; ## if $debug;

  my $stmt = $dbh->prepare($sql);
  $stmt->execute() || die "SQL ERROR: $stmt->errstr()";
	
  my $miniLib;
  my $count = 0;
  ##following segment for pre fetching theids...is less efficient than row...
  my @ids;
  while (my ($id) = $stmt->fetchrow_array()) {
    #    next unless $id > 1000000;
    next if exists $finished{$id};
    $count++;
    last if ($ctx->{'cla'}->{'testnumber'} && $count >= $ctx->{'cla'}->{'testnumber'}); ##testing..
    print STDERR "fetching $count ids to process\n" if $count % 10000 == 0;
    push(@ids,$id);
  }
  print STDERR "\nMaking $count AssemblySequences\n\n";
  return if $count == 0;
  $count = 0;
  foreach my $id (@ids) {
    my $ex = ExternalNASequence->new({'na_sequence_id' => $id});
    next unless $ex->retrieveFromDB();
    ##want to just use the quality_sequence if it exists...
    ##there could be more than one lenssequence child...want newest one...
    my @ls_sort = sort{$b->getId() <=> $a->getId()} $ex->getChildren('LENSSequence',1);
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
      $miniLib .= ">".$ex->getId()."\n".SequenceUtils::breakSequence(substr($ex->getSequence(),$qStart,$qLength));
    }else{
      $miniLib .= $ex->toFasta();
    }
    $count++;
    if ($count >= 10000) {
      &processSet($miniLib);
      $miniLib = "";            ##reset for next set of seqs
      $count = 0;
      $ctx->{'self_inv'}->undefPointerCache();
    }
  }
  &processSet($miniLib) if $miniLib;        ##processes last set
  $ctx->{'self_inv'}->undefPointerCache();
  $stmt->finish();              ##cancels when testing so can do second query....
}

sub processSet {
  my($miniLib) = @_;
  open(S, ">tmpLib");
  print S $miniLib;
  close S;
  
  ##NOTE that need to get this installed correctly...
  #  system("/usr/local/src/bio/PHRAP_etal/phrap.SUNWspro/ultra.bin/cross_match tmpLib /usr/local/db/others/repeat/vector -screen > cross.test");
  system("$phrap_dir/cross_match tmpLib $library -screen > cross.test 2> /dev/null");

  ##generate better sequence....
  open(S,"tmpLib.screen");
  my $seq;
  my $na_seq_id;
  while (<S>) {
    if (/^\>(\d+)/) {           ##$1 = na_sequence_id
      &makeAndInsertAssSeq($na_seq_id,$seq) if $na_seq_id;
      $na_seq_id = $1;
      $seq = "";
      print STDERR "Processed: $countProcessed, low_quality: $countBad ",($countProcessed % 1000 == 0 ? `date` : "\n") if $countProcessed % 100 == 0;
    } else {
      $seq .= $_;
    }
  }
  &makeAndInsertAssSeq($na_seq_id,$seq) if $na_seq_id;
  close S;
}

sub makeAndInsertAssSeq {
  my($na_seq_id,$seq) = @_;
  my $qseq = &returnQuality($seq);

  my $ass = AssemblySequence->new( { 'na_sequence_id' => $na_seq_id,
                                     'sequence_version' => 1,
                                     'assembly_offset' => 0,
                                     'assembly_strand' => 1 } );
  $ass->setSequence($qseq);
  $ass->setQualityStart($ass->getSequenceStart());
  $ass->setQualityEnd($ass->getSequenceEnd());
  if (length($qseq) < 50) {     ##poor quality sequence
    $ass->set('have_processed',1); ##mark processed so do not cluster
    $ass->set('processed_category','low_quality');
    $countBad++;
  } else {
    $ass->set('have_processed',0);
  }
  print $ass->toString(0,1) if $debug;
  $ass->submit();
  print EX $ass->toFasta(1) if ($ctx->{cla}->{commit} && $ctx->{cla}->{export} && $ass->getHaveProcessed() == 0);
  $countProcessed++;
}

##trims leadint TTTT and trailing AAAAA and scans for stretch of good sequence
##by looking at the %Ns in a 20 bp window...must be less than 20%
sub returnQuality{
  my($seq) = @_;
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
        return &trimAT($newSeq); 
      }
    }
  }
  $newSeq = substr($seq, $newstart, length($seq));
  return &trimAT($newSeq);
}

##steps through sequence with 30 bp window...if >=36 As or Ts then truncates
##either beginning or end (whichever is closer)
sub trimAT {
  my($seq) = @_;
  $seq =~ s/\s//g;
  print STDERR "\ntrimAT input: \n", SequenceUtils::breakSequence($seq) if $debug == 1;
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

__END__

=pod
=head1 Description
B<Template> - a template plug-in for C<ga> (GUS application) package.

=head1 Purpose
B<Template> is a minimal 'plug-in' GUS application.

=cut
