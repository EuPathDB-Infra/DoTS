############################################################
## Change Package name....
############################################################
package ExtractAndBlockAssemblySequences;

use strict;

############################################################
# Add any specific objects (GUSdev::) here
############################################################

use Objects::GUSdev::AssemblySequence;
##for query handle on database....
use DBI;

sub new {
  my $Class = shift;

  return bless {}, $Class;
}

sub Usage {
  my $M   = shift;
  return 'Extracts unprocessed AssembySequences, blocks them and  writes to a file for clustering';
}

############################################################
# put the options in this method....
############################################################
sub EasyCspOptions {
  my $M   = shift;
  {

  testnumber        => {
                        o => 'testnumber=i',
                        h => 'number of iterations for testing',
                       },
                         
  taxon_id          => {
                        o => 'taxon_id=i',
                        h => 'taxon_id for sequences to process: 8=Hum, 14=Mus.',
                       },

  outputfile        => {
                        o => 'outputfile=s',
                        h => 'Name of file for output sequences',
                       },
                         
  rm_options        => {
                        o => 'rm_options',
                        t => 'string',
                        h => 'RepeatMasker options',
                       },
                         
  idSQL             => {
                        o => 'idSQL',
                        t => 'string',
                        h => 'SQL query that returns assembly_sequence_ids to be processed',
                       },
  extractonly       => {
                        o => 'extractonly',
                        t => 'boolean',
                        h => 'if true then does not Block extracted sequences',
                       },
                         
                     }
}

my $ctx;
my $countProcessed = 0;
my $countBad = 0;
my $repLib;
my $debug = 0;
my $repMaskDir = '/usr/local/src/bio/RepeatMasker/04-04-1999';
my $tmpLib = "tmpLib.$$";
my $RepMaskCmd;

sub Run {
  my $M   = shift;
  $ctx = shift;

  die "You must enter repeat masker options on the command line to specify minimally the organism to be blocked\n" unless $ctx->{cla}->{rm_options} || $ctx->{cla}->{extractonly};
  die "You must enter either the --taxon_id or --idSQL on the command line\n" unless $ctx->{cla}->{taxon_id} || $ctx->{cla}->{idSQL};
  die "You must enter an outputfile name on the command line\n" unless $ctx->{cla}->{'outputfile'};

  print STDERR $ctx->{cla}->{'commit'} ? "***COMMIT ON***\n" : "COMMIT TURNED OFF\n";
  print STDERR "Testing on $ctx->{cla}->{'testnumber'}\n" if $ctx->{cla}->{'testnumber'};

  my $dbh = $ctx->{self_inv}->getQueryHandle();

  if ($ctx->{cla}->{extractonly}) {
    print STDERR "Extracting sequences without blocking\n";
  } else {
    $RepMaskCmd = "RepeatMasker $ctx->{cla}->{rm_options}";
    print STDERR "RepeatMasker command:\n  $RepMaskCmd\n";
  }

  ##implement restart here....
  my %finished;
  ##restart 
  if ( -e "$ctx->{cla}->{'outputfile'}") {
    open(F,"$ctx->{cla}->{'outputfile'}");
    while (<F>) {
      if (/^\>(\d+)/) {
        $finished{$1} = 1;
      }
    }
    close F;
    open(OUT,">>$ctx->{cla}->{'outputfile'}");
    print STDERR "outputFile $ctx->{'outputfile'} exists...Restarting: already have ".scalar(keys%finished)." sequences\n";
  } else {
    open(OUT,">$ctx->{cla}->{'outputfile'}");
  }

  my $getSeqs;
  if ($ctx->{cla}->{idSQL}) {
    $getSeqs = $ctx->{cla}->{idSQL};
  } else {
    $getSeqs = "select a.assembly_sequence_id from AssemblySequence a, ExternalNASequence e
  where a.have_processed = 0 and a.na_sequence_id = e.na_sequence_id
  and e.taxon_id = $ctx->{'taxon_id'} and a.quality_end - a.quality_start >= 50";
  }

  print STDERR "$getSeqs\n" if $debug;

  my $stmt = $dbh->prepare($getSeqs);
  $stmt->execute();
  my $count = 0;
  my $miniLib = "";
  my @todo;

  ##run it into an array so does not block!!
  while (my($id) = $stmt->fetchrow_array()) {
    next if exists $finished{$id};
    last if ($ctx->{'testnumber'} && $count >= $ctx->{'testnumber'}); ##breaks 
    print STDERR "Retrieving $count\n" if $count % 10000 == 0;
    push(@todo,$id);
    $count++;
  }
  print STDERR "Extracting",($ctx->{cla}->{extractonly} ? " " : " and blocking "),"$count sequences from taxon_id $ctx->{'taxon_id'}\n";

  $count = 0;
  my $countProc = 0;
  my $reset = 0;
  foreach my $id (@todo) {
    print STDERR "Processing $id\n" if $debug;
    my $ass = AssemblySequence->new( { 'assembly_sequence_id' => $id } );
    $ass->retrieveFromDB();

    ##want to set the sequence_start = quality_start etc here....has not been assembled...
    $ass->resetAssemblySequence();
#    $reset += $ass->submit() if $ass->hasChangedAttributes();

    $miniLib .= $ass->toFasta(1);
    $count++;
    $countProc++;
    print STDERR "Processing $countProc\n" if $countProc % 1000 == 0;
    if ($count >= 1000) {
      &processSet($miniLib);
      $miniLib = "";            ##reset for next set of seqs
      $count = 0;
      $ctx->{'self_inv'}->undefPointerCache();
    }
  }
  &processSet($miniLib);        ##processes last set

  ##clean up after self...
  unlink "$tmpLib";
  unlink "$tmpLib.masked";

  ############################################################
  ###  put an informative summary in the results variable
  ############################################################
  my $results = "Extracted and blocked $countProcessed AssemblySequences, marked $countBad as repeat and reset $reset to qality_start/end";
  $results = "Extracted $countProc AssemblySequences" if $ctx->{cla}->{extractonly};
  print STDERR "\n$results\n";
  return $results;
}

sub processSet {
  my($miniLib) = @_;

  if ($ctx->{cla}->{extractonly}) {
    print OUT $miniLib;
    return;
  }

  open(S, ">$tmpLib");
  print S $miniLib;
  close S;
  
  ##RepeatMasker
  system("$repMaskDir/$RepMaskCmd $tmpLib");

  ##generate better sequence....
  open(S,"$tmpLib.masked");
  my $seq;
  my $na_seq_id;
  while (<S>) {
    if (/^\>(\d+)/) {           ##$1 = na_sequence_id
      &processBlockedSequence($na_seq_id,$seq) if $na_seq_id;
      $na_seq_id = $1;
      $seq = "";
      print STDERR "Processed: $countProcessed, repeats: $countBad\n" if $countProcessed % 100 == 0;
    } else {
      $seq .= $_;
    }
  }
  &processBlockedSequence($na_seq_id,$seq) if $na_seq_id;
  close S;
}

sub processBlockedSequence{
  my($ass_seq_id,$seq) = @_;

  $countProcessed++;

  $seq =~ s/\s+//g;
  $seq =~ s/X/N/g;
  ##trim dangling NNNN.s
  my $sequence = &trimDanglingNNN($seq);

  ##check for lenth..
  my $tmpSeq = $sequence;
  $tmpSeq =~ s/N//g;
	
  ##if too short then update AssemblySquence else print to file...
  if (length($tmpSeq) < 50) {
    print STDERR "Sequence $ass_seq_id too short (".length($tmpSeq).") following blocking\n" if $debug;
    ##update AssSeq..
    my $ass = $ctx->{'self_inv'}->getFromDbCache('AssemblySequence',$ass_seq_id);
    if (!$ass) {
      print STDERR "ERROR: $ass_seq_id not in cache...retrieving from Database\n";
      $ass = AssemblySequence->new( { 'assembly_sequence_id' => $ass_seq_id });
      $ass->retrieveFromDB();
      if (!$ass->get('assembly_strand')) {
				##this is invalid sequence.....is reverse strand..
        print STDERR "ERROR:  AssemblySequence $ass_seq_id is invalid\n";
        return undef;
      }
    }
    $ass->set('have_processed',1);
    $ass->set('processed_category','repeat');
    $ass->submit();
    $countBad++;
  } else {
    print OUT "\>$ass_seq_id\n".CBIL::Bio::SequenceUtils::breakSequence($sequence);
  }
}

sub trimDanglingNNN {
  my($seq) = @_;
  if ($seq =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
    $seq = $2 if length($1) < 20; ##don't want to leave 20 bp at end...
  }

  if ($seq =~ /N/) {            ##still has at least one N so..
    my $rev = CBIL::Bio::SequenceUtils::reverseComplementSequence($seq);
    if ($rev =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
      #			print "matched ending NNNN length\$1=".length($1)." length\$2=".length($2)."\nSEQ:$seq\n";
      $rev = $2 if length($1) < 20; ##don't want to leave 20 bp at end...
    }
    if (length($rev) == length($seq)) {
      return $seq;
    } else {
      return CBIL::Bio::SequenceUtils::reverseComplementSequence($rev);
    }
  } else {
    return $seq;
  }
}


1;
