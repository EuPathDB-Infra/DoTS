package DoTS::DotsBuild::Plugin::ExtractAndBlockAssemblySequences;


@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use GUS::Model::DoTS::AssemblySequence;
use CBIL::Bio::SequenceUtils;

sub new {

  my $class = shift;

  my $self = {};
  bless($self,$class);

  my $usage = 'Extract unprocessed AssembySequences, block them and write to a file for clustering';

  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'taxon_id_list',
      t => 'string',
      h => 'comma delimited list of taxon_ids for sequences to process: 8=Hum, 14=Mus.',
     },
     {o => 'outputfile',
      t => 'string',
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
		  name => ref($self),
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

sub run {
    my $self   = shift;

    my $cla =$self->getCla; 
    
    die "You must enter repeat masker options on the command line to specify minimally the organism to be blocked\n" unless $cla->{rm_options} || $cla->{extractonly};
    die "You must enter either the --taxon_id_list or --idSQL on the command line\n" unless $cla->{taxon_id_list} || $cla->{idSQL};
    die "You must enter an outputfile name on the command line\n" unless $cla->{'outputfile'};
    
    $self->logAlert ($cla->{'commit'} ? "***COMMIT ON***\n" : "COMMIT TURNED OFF\n");
    $self->logAlert ("Testing on $cla->{'testnumber'}\n") if $cla->{'testnumber'};
    
    my $dbh = $self->getQueryHandle();
    
    if ($cla->{extractonly}) {
	$self->logAlert("Extracting sequences without blocking\n");
    } else {
	$RepMaskCmd = "RepeatMasker $cla->{rm_options}";
	$self->logAlert ("RepeatMasker command:\n  $RepMaskCmd\n");
    }
    
    ##implement restart here....
    my %finished;
    ##restart 
    if ( -e "$cla->{'outputfile'}") {
	open(F,"$cla->{'outputfile'}");
	while (<F>) {
	    if (/^\>(\d+)/) {
		$finished{$1} = 1;
	    }
	}
	close F;
	open(OUT,">>$cla->{'outputfile'}");
	$self->logAlert ("outputFile $cla->{'outputfile'} exists...Restarting: already have ".scalar(keys%finished)." sequences\n");
    } else {
	open(OUT,">$cla->{'outputfile'}");
    }
    
    my $getSeqs;
    if ($cla->{idSQL}) {
	$getSeqs = $cla->{idSQL};
    } else {
    my $taxon = $cla->{'taxon_id_list'};
	$getSeqs = "select a.assembly_sequence_id 
  from dots.AssemblySequence a, dots.ExternalNASequence e
  where a.have_processed = 0 
  and a.na_sequence_id = e.na_sequence_id
  and e.taxon_id in (". $taxon .") and a.quality_end - a.quality_start >= 50";
    }
    
    $self->logVerbose ("$getSeqs\n");
    
    my $stmt = $dbh->prepare($getSeqs);
    $stmt->execute();
    my $count = 0;
    my $miniLib = "";
    my @todo;
    
    ##run it into an array so does not block!!
    while (my($id) = $stmt->fetchrow_array()) {
	next if exists $finished{$id};
	last if ($cla->{'testnumber'} && $count >= $cla->{'testnumber'}); ##breaks 
	$self->logAlert ("Retrieving $count\n") if $count % 10000 == 0;
	push(@todo,$id);
	$count++;
    }
    $self->logAlert ("Extracting",($cla->{extractonly} ? " " : " and blocking "),"$count sequences from taxon_id(s) $cla->{'taxon_idist'}\n");
    
    $count = 0;
    my $countProc = 0;
    my $reset = 0;
    foreach my $id (@todo) {
	$self->logAlert ("Processing $id\n") if $debug;
	my $ass = GUS::Model::DoTS::AssemblySequence->
	    new( { 'assembly_sequence_id' => $id } );
	$ass->retrieveFromDB();
	
	##want to set the sequence_start = quality_start etc here....has not been assembled...
	$ass->resetAssemblySequence();
#    $reset += $ass->submit() if $ass->hasChangedAttributes();
	
	$miniLib .= $ass->toFasta(1);
	$count++;
	$countProc++;
	$self->logAlert ("Processing $countProc\n") if $countProc % 1000 == 0;
	if ($count >= 1000) {
	    $self->processSet($miniLib);
	    $miniLib = "";            ##reset for next set of seqs
	    $count = 0;
	    $self->undefPointerCache();
	}
    }
    $self->processSet($miniLib);        ##processes last set
    
    ##clean up after self...
    unlink "$tmpLib";
    unlink "$tmpLib.masked";
    
    ############################################################
    ###  put an informative summary in the results variable
    ############################################################
    my $results = "Extracted and blocked $countProcessed AssemblySequences, marked $countBad as repeat and reset $reset to qality_start/end";
    $results = "Extracted $countProc AssemblySequences" if $cla->{extractonly};
    $self->logAlert ("\n$results\n");
    return $results;
}

sub processSet {

    my $self   = shift;

    my($miniLib) = @_;

    my $cla = $self->getCla;
    if ($cla->{extractonly}) {
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
	    $self->processBlockedSequence($na_seq_id,$seq) if $na_seq_id;
	    $na_seq_id = $1;
	    $seq = "";
	    $self->logAlert ("Processed: $countProcessed, repeats: $countBad\n") if $countProcessed % 100 == 0;
	} else {
	    $seq .= $_;
	}
    }
    $self->processBlockedSequence($na_seq_id,$seq) if $na_seq_id;
    close S;
}

sub processBlockedSequence{


    my $self   = shift; 
    my $ctx = shift;
    my($ass_seq_id,$seq) = @_;
    
    $countProcessed++;
    
    $seq =~ s/\s+//g;
    $seq =~ s/X/N/g;
    ##trim dangling NNNN.s
    my $sequence = $self->trimDanglingNNN($seq);
    
    ##check for lenth..
    my $tmpSeq = $sequence;
    $tmpSeq =~ s/N//g;
    
    ##if too short then update AssemblySquence else print to file...
    if (length($tmpSeq) < 50) {
	$self->logAlert ("Sequence $ass_seq_id too short (".length($tmpSeq).") following blocking\n") if $debug;
	##update AssSeq..
	my $ass = $ctx->{'self_inv'}->getFromDbCache('AssemblySequence',$ass_seq_id);
	if (!$ass) {
	    $self->logAlert ("ERROR: $ass_seq_id not in cache...retrieving from Database\n");
	    $ass = GUS::Model::DoTS::AssemblySequence->
		new( { 'assembly_sequence_id' => $ass_seq_id });
	    $ass->retrieveFromDB();
	    if (!$ass->get('assembly_strand')) {
		##this is invalid sequence.....is reverse strand..
		$self->logAlert ("ERROR:  AssemblySequence $ass_seq_id is invalid\n");
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

    my $self   = shift;
    my($seq) = @_;
    if ($seq =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
	$seq = $2 if length($1) < 20; ##don't want to leave 20 bp at end...
    }
    
    if ($seq =~ /N/) {            ##still has at least one N so..
	my $rev = CBIL::Bio::SequenceUtils::reverseComplementSequence($seq);
	if ($rev =~ /^(.*?)NNNNNNNNNN+(.*?)$/) {
	    #	 $self->logAlert ("matched ending NNNN length\$1=".length($1)." length\$2=".length($2)."\nSEQ:$seq\n");
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
