package DoTS::DotsBuild::Plugin::ExtractAndBlockAssemblySequences;


@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use GUS::Model::DoTS::AssemblySequence;
use CBIL::Bio::SequenceUtils;

sub new {
  my $Class = shift;

  my $self = {};
  bless($self,$class);

  my $usage = 'Extract unprocessed AssembySequences, block them and write to a file for clustering';

  my $easycsp =
    [
     {o => 'testnumber=i',
      h => 'number of iterations for testing',
     },
     {o => 'taxon_id=i',
      h => 'taxon_id for sequences to process: 8=Hum, 14=Mus.',
     },
     {o => 'outputfile=s',
      h => 'Name of file for output sequences',
     },
     {o => 'rm_options',
      t => 'string',
      h => 'RepeatMasker options',
     },
     {o => 'idSQL',
      t => 'string',
      h => 'SQL query that returns assembly_sequence_ids to be processed',
     },
     {o => 'extractonly',
      t => 'boolean',
      h => 'if true then does not Block extracted sequences',
     }
    ];

  $self->initialize({requiredDbVersion => {},
		  cvsRevision => '$Revision$', # cvs fills this in!
		  cvsTag => '$Name$', # cvs fills this in!
		  name => ref($m),
		  revisionNotes => 'make consistent with GUS 3.0',
		  easyCspOptions => $easycsp,
		  usage => $usage
		 });

  return $self;
}


my $countProcessed = 0;
my $countBad = 0;
my $repLib;
my $debug = 0;
my $repMaskDir = '/usr/local/src/bio/RepeatMasker/04-04-1999';
my $tmpLib = "tmpLib.$$";
my $RepMaskCmd;

sub Run {
    my $M   = shift;
    
    die "You must enter repeat masker options on the command line to specify minimally the organism to be blocked\n" unless $M->getCla{rm_options} || $M->getCla{extractonly};
    die "You must enter either the --taxon_id or --idSQL on the command line\n" unless $M->getCla{taxon_id} || $M->getCla{idSQL};
    die "You must enter an outputfile name on the command line\n" unless $M->getCla{'outputfile'};
    
    $M->log ($M->getCla{'commit'} ? "***COMMIT ON***\n" : "COMMIT TURNED OFF\n");
    $M->log ("Testing on $M->getCla{'testnumber'}\n") if $M->getCla{'testnumber'};
    
    my $dbh = $M->getQueryHandle();
    
    if ($M->getCla{extractonly}) {
	$M->log("Extracting sequences without blocking\n");
    } else {
	$RepMaskCmd = "RepeatMasker $M->getCla{rm_options}";
	$M->log ("RepeatMasker command:\n  $RepMaskCmd\n");
    }
    
    ##implement restart here....
    my %finished;
    ##restart 
    if ( -e "$M->getCla{'outputfile'}") {
	open(F,"$M->getCla{'outputfile'}");
	while (<F>) {
	    if (/^\>(\d+)/) {
		$finished{$1} = 1;
	    }
	}
	close F;
	open(OUT,">>$M->getCla{'outputfile'}");
	$M->log ("outputFile $M->getCla{'outputfile'} exists...Restarting: already have ".scalar(keys%finished)." sequences\n");
    } else {
	open(OUT,">$M->getCla{'outputfile'}");
    }
    
    my $getSeqs;
    if ($M->getCla{idSQL}) {
	$getSeqs = $M->getCla{idSQL};
    } else {
	$getSeqs = "select a.assembly_sequence_id 
  from dots.AssemblySequence a, dots.ExternalNASequence e
  where a.have_processed = 0 
  and a.na_sequence_id = e.na_sequence_id
  and e.taxon_id = $M->getCla{'taxon_id'} 
  and a.quality_end - a.quality_start >= 50";
    }
    
    $M->log ("$getSeqs\n") if $debug;
    
    my $stmt = $dbh->prepare($getSeqs);
    $stmt->execute();
    my $count = 0;
    my $miniLib = "";
    my @todo;
    
    ##run it into an array so does not block!!
    while (my($id) = $stmt->fetchrow_array()) {
	next if exists $finished{$id};
	last if ($M->getCla{'testnumber'} && $count >= $M->getCla{'testnumber'}); ##breaks 
	$M->log ("Retrieving $count\n") if $count % 10000 == 0;
	push(@todo,$id);
	$count++;
    }
    $M->log ("Extracting",($M->getCla{extractonly} ? " " : " and blocking "),"$count sequences from taxon_id $M->getCla{'taxon_id'}\n");
    
    $count = 0;
    my $countProc = 0;
    my $reset = 0;
    foreach my $id (@todo) {
	$M->log ("Processing $id\n") if $debug;
	my $ass = GUS::Model::DoTS::AssemblySequence->
	    new( { 'assembly_sequence_id' => $id } );
	$ass->retrieveFromDB();
	
	##want to set the sequence_start = quality_start etc here....has not been assembled...
	$ass->resetAssemblySequence();
#    $reset += $ass->submit() if $ass->hasChangedAttributes();
	
	$miniLib .= $ass->toFasta(1);
	$count++;
	$countProc++;
	$M->log ("Processing $countProc\n") if $countProc % 1000 == 0;
	if ($count >= 1000) {
	    &processSet($miniLib);
	    $miniLib = "";            ##reset for next set of seqs
	    $count = 0;
	    $M->undefPointerCache();
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
    $results = "Extracted $countProc AssemblySequences" if $M->getCla{extractonly};
    $M->log ("\n$results\n");
    return $results;
}

sub processSet {

    my $M   = shift;
    my($miniLib) = @_;
    
    if ($M->getCla{extractonly}) {
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
	    $M->log ("Processed: $countProcessed, repeats: $countBad\n") if $countProcessed % 100 == 0;
	} else {
	    $seq .= $_;
	}
    }
    &processBlockedSequence($na_seq_id,$seq) if $na_seq_id;
    close S;
}

sub processBlockedSequence{

    my $M   = shift; 
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
	$M->log ("Sequence $ass_seq_id too short (".length($tmpSeq).") following blocking\n") if $debug;
	##update AssSeq..
	my $ass = $ctx->{'self_inv'}->getFromDbCache('AssemblySequence',$ass_seq_id);
	if (!$ass) {
	    $M->log ("ERROR: $ass_seq_id not in cache...retrieving from Database\n");
	    $ass = GUS::Model::DoTS::AssemblySequence->
		new( { 'assembly_sequence_id' => $ass_seq_id });
	    $ass->retrieveFromDB();
	    if (!$ass->get('assembly_strand')) {
		##this is invalid sequence.....is reverse strand..
		$M->log ("ERROR:  AssemblySequence $ass_seq_id is invalid\n");
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

    my $M   = shift;
    my($seq) = @_;
    if ($seq =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
	$seq = $2 if length($1) < 20; ##don't want to leave 20 bp at end...
    }
    
    if ($seq =~ /N/) {            ##still has at least one N so..
	my $rev = CBIL::Bio::SequenceUtils::reverseComplementSequence($seq);
	if ($rev =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
	    #	 $M->log ("matched ending NNNN length\$1=".length($1)." length\$2=".length($2)."\nSEQ:$seq\n");
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
