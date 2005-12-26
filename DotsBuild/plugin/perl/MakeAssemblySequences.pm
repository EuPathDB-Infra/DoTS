package DoTS::DotsBuild::Plugin::MakeAssemblySequences;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;

use GUS::Model::DoTS::ExternalNASequence;
use GUS::Model::DoTS::AssemblySequence;
use CBIL::Bio::SequenceUtils;
use GUS::PluginMgr::Plugin;

my $argsDeclaration =
[
 integerArg({name => 'testnumber1',
	     descr => 'number of iterations for testing',
	     constraintFunc => undef,
	     reqd => 0,
	     isList => 0
	     }),

 stringArg({name => 'date',
	    descr => 'earliest date for sequences to include: default = all',
	    constraintFunc => undef,
	    reqd => 0,
	    isList => 0
	    }),

 stringArg({name => 'taxon_id_list',
	    descr => 'comma delimited taxon_id list for sequences to process: 8=Hum, 14=Mus.',
	    constraintFunc => undef,
	    reqd => 1,
	    isList => 0
	    }),

 stringArg({name => 'idSQL',
	    descr => 'SQL query that returns na_sequence_ids from ExternalNASequence to be processed',
	    constraintFunc => undef,
	    reqd => 0,
	    isList => 0
	    }),

 fileArg({name => 'export',
	    descr => 'filename to export sequences to...default does not export',
	    constraintFunc => undef,
	    reqd => 0,
	    mustExist => 0,
	    isList => 0,
    	format => 'Text'
	    }),

 fileArg({name => 'repeatFile',
	  descr => 'full path of file of repeats',
	  constraintFunc => undef,
	  reqd => 1,
	  isList => 0,
	  mustExist => 1,
	  format => 'Text'
        }),

 fileArg({name => 'phrapDir',
	  descr => 'full path of the directory containing phrap\'s cross_match program',
	  constraintFunc => undef,
	  reqd => 1,
	  isList => 0,
	  mustExist => 1,
	  format => 'Directory'
        }),


 ];


my $purposeBrief = <<PURPOSEBRIEF;
Retrieves ExternalNASequences that are not assembly sequences and inserts AssemblySequence entry
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
Retrieves ExternalNASequences that are not assembly sequences and inserts AssemblySequence entry
PLUGIN_PURPOSE

#check the documentation for this
my $tablesAffected = [];

my $tablesDependedOn = [
    ['DoTS::ExternalNASequence', ''],
    ['DoTS::AssemblySequence', '']
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

my $debug = 0;

my $countProcessed = 0;
my $countBad = 0;

my %finished;
my $library;

$| = 1;

sub run {
  my $self = shift;
  my $ctx = shift;  
  print STDERR ($self->getCla->{'commit'} ? "COMMIT ON\n" : "COMMIT TURNED OFF\n");
  print STDERR ("Testing on". $self->getCla->{'testnumber'}."\n") if $self->getCla->{'testnumber'};
  
  ##set the taxon_id_list...
  die "You must enter either the --taxon_id_list and optionally --idSQL on the command line\n" unless $self->getCla->{taxon_id_list};

  die "You must provide the subdirectory and repeat file, e.g. unknown_release/vector_humMitoRibo.lib\n" unless $self->getCla->{repeatFile};

  die "You must provide a --phrapDir in which an executable cross_match program resides" unless -x $self->getCla->{phrapDir} . "/cross_match";

  $library = $self->getCla->{repeatFile};
  print STDERR ("Running cross_match with $library\n");
  
  $ctx->{'self_inv'}->setMaximumNumberOfObjects(100000);
  
  my $dbh = $self->getQueryHandle();
  
  if ($self->getCla->{export}) {
    my $export = $self->getCla->{export}; 
    if (-e $self->getCla->{export}) {
      open(F,"$export");
      while (<F>) {
	if (/^\>(\S+)/) {
	  $finished{$1} = 1;
	}
      }
      close F;
      print STDERR ("Already processed ",scalar(keys%finished)," ids\n");
    }
    open(EX,">>$export");
  }
  
  if ($self->getCla->{idSQL}) {
    $self->processQuery($self->getCla->{idSQL});
  } else {
    
    print STDERR ("Generating new AssemblySequences for taxon(s)". $self->getCla->{taxon_id_list}."\n");
    
    my $taxonIdList = $self->getCla->{taxon_id_list};
    
    ##first get the ESTs and mRNAs..
    my $sql =
        "select e.na_sequence_id 
         from dots.ExternalNASequence e, dots.sequencetype st 
         where e.taxon_id in($taxonIdList)
         and st.name in ('mRNA', 'EST')
         and e.sequence_type_id = st.sequence_type_id    
         and e.na_sequence_id not in 
         (select a.na_sequence_id from dots.AssemblySequence a) ";
    
    if ($self->getCla->{'date'}) {
      $sql .= "and modification_date >= '".$self->getCla->{'date'}."'";
    } 
    
    $self->processQuery($sql);
    
    ##next get the things from embl that are RNAs longer than 500 bp...
    ##need to check this for things that are not human or mouse...may need to use less sophisticated query!
    my $mRNASql = 
        "select o.na_sequence_id 
        from dots.externalnasequence o 
        where o.na_sequence_id 
        in ( select s.na_sequence_id 
             from dots.externalnasequence s, dots.transcript t, 
             dots.nalocation l, dots.sequencetype st
             where s.taxon_id in ($taxonIdList) 
             and st.name = 'RNA'
             and s.sequence_type_id = st.sequence_type_id    
             and t.na_sequence_id = s.na_sequence_id
             and t.name = 'CDS'
             and l.na_feature_id = t.na_feature_id
             group by s.na_sequence_id having count(*) = 1 )
             and o.length > 400
             and o.na_sequence_id 
             not in (select a.na_sequence_id from dots.AssemblySequence a)";

    
    if ($self->getCla->{'date'}) {
      $mRNASql .= "and modification_date >= '".$self->getCla->{'date'}."'";
    } 
    $self->processQuery($mRNASql);
    
  }
  
  # unlink "tmpLib";
  # unlink "tmpLib.mask";
  
  close EX;
  my $results = "Processed $countProcessed AssemblySequences and marked $countBad as 'low_quality'";
  
  print STDERR ("\n$results\n");
  
  return $results;
}

sub processQuery {
  
  my $self   = shift;
  my($sql) = @_;
  
  print STDERR ("\n$sql\n"); ## if $debug;
  my $dbh = $self->getQueryHandle();
  
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
    last if ($self->getCla->{'testnumber'} && $count >= $self->getCla->{'testnumber'}); ##testing..
    print STDERR ("fetching $count ids to process\n") if $count % 10000 == 0;
    push(@ids,$id);
  }
  print STDERR ("\nMaking $count AssemblySequences\n\n");
  return if $count == 0;
  $count = 0;
  foreach my $id (@ids) {
    my $ex = GUS::Model::DoTS::ExternalNASequence->new({'na_sequence_id' => $id});
    next unless $ex->retrieveFromDB();
    ##want to just use the quality_sequence if it exists...
    my ($ls) = $ex->getChildren('GUS::Model::DoTS::EST',1);
    my $qStop;
    $qStop = $ls->getQualityStop() if $ls;
    if($ls && defined $qStop){  ##exists and has qualityStop...assume quality_start is not relevant if qualitystop is null..
      my $qStart = $ls->getQualityStart() > 0 ? $ls->getQualityStart() - 1 : 0;
      ##turns out that $qStart may be > length of the sequence...set to the length...
      if($qStart > $ex->getLength() - 1){$qStart =  $ex->getLength() - 1;}
      ##set length to 20 if less than thatt so cross_match will work...will be removed from AssembySequence
      ##since is less than 50 bp.
      my $qLength = $qStop == 0 ? $ex->getLength() - $qStart : ($qStop - $qStart < 20 ? 20 : $qStop - $qStart);
      $miniLib .= ">".$ex->getId()."\n".CBIL::Bio::SequenceUtils::breakSequence(substr($ex->getSequence(),$qStart,$qLength));
    }else{
      $miniLib .= $ex->toFasta();
    }
    $count++;
    if ($count >= 10000) {
      $self->processSet($miniLib);
      $miniLib = "";            ##reset for next set of seqs
      $count = 0;
      $self->undefPointerCache();
    }
  }
  $self->processSet($miniLib) if $miniLib;        ##processes last set
  $self->undefPointerCache();
  $stmt->finish();              ##cancels when testing so can do second query....
}

sub processSet {
  
  my $self   = shift;
  my($miniLib) = @_;
  
  open(S, ">tmpLib") || die "Can't open tmpLib";
  print S $miniLib;
  close S;

  my $phrapDir = $self->getCla->{phrapDir};
    ##NOTE that need to get this installed correctly...
  #  system("/usr/local/src/bio/PHRAP_etal/phrap.SUNWspro/ultra.bin/cross_match tmpLib /usr/local/db/others/repeat/vector -screen > cross.test");
  die unless -x "$phrapDir/cross_match";
  my $retCode = system("$phrapDir/cross_match tmpLib $library -screen > cross.test 2> cross_match.err");  
  die "Failed with status $retCode running cross_match" if $retCode;

  ##generate better sequence....
  open(S,"tmpLib.screen") || die "Can't open tmpLib.screen";
  my $seq;
  my $na_seq_id;
  while (<S>) {
    if (/^\>(\d+)/) {           ##$1 = na_sequence_id
      $self->makeAndInsertAssSeq($na_seq_id,$seq) if $na_seq_id;
      $na_seq_id = $1;
      $seq = "";
      print STDERR ("Processed: $countProcessed, low_quality: $countBad ",($countProcessed % 1000 == 0 ? `date` : "\n")) if $countProcessed % 100 == 0;
    } else {
      $seq .= $_;
    }
  }
  $self->makeAndInsertAssSeq($na_seq_id,$seq) if $na_seq_id;
  close S;
}

sub makeAndInsertAssSeq {
  
  my $self   = shift;
  my($na_seq_id,$seq) = @_;
  my $qseq = $self->returnQuality($seq);
  
  my $ass = GUS::Model::DoTS::AssemblySequence->new( { 'na_sequence_id' => $na_seq_id,
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
  print STDERR ($ass->toString(0,1)) if $debug;
  $ass->submit();
  print EX $ass->toFasta(1) if ($self->getCla->{commit} && $self->getCla->{export} && $ass->getHaveProcessed() == 0);
  $countProcessed++;
}

##trims leadint TTTT and trailing AAAAA and scans for stretch of good sequence
##by looking at the %Ns in a 20 bp window...must be less than 20%
sub returnQuality{
  
  my $self   = shift;
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
	#	print STDERR ("Trimming Poor qual beginning: newStart set to $newstart\n");
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
  
  my $self   = shift;
  my($seq) = @_;
  $seq =~ s/\s//g;
  print STDERR ("\ntrimAT input: \n", CBIL::Bio::SequenceUtils::breakSequence($seq)) if $debug == 1;
  my @seq = split('', $seq);
  my %nuc;
  my($startBase,$endBase);
  foreach my $n (@seq[0..29]) {
    $nuc{$n}++;
  }
  print STDERR  ("First 30: ", join('', @seq[0..29]), "\n") if $debug == 1;
  my $inRun = 0;
  for (my $i=1;$i<(length($seq)-30);$i++) {
    if ($inRun == 0 && ($nuc{A} >=27 || $nuc{T} >= 27)) {
      $startBase = $i;
      print STDERR  ("A=$nuc{A},T=$nuc{T} Sequence Length: ", length($seq), " StartBase: $startBase ") if $debug == 1;
      $inRun = 1;
      next;
    }
    if ($inRun == 1 && ($nuc{A} < 27 && $nuc{T} < 27)) {
      ##process the sucker!!
      $endBase = $i + 30;
      print STDERR ("EndBase: $endBase\n") if $debug == 1;
      ##test if closer to beginning or end and tructate accordingly
      if ($startBase <= (length($seq) - $endBase)) {
	print STDERR ("substr\(\$seq,\($endBase-5\),\(length\(\$seq\)\)\)\n") if $debug == 1;
	return substr($seq,($endBase-5),(length($seq)));
      } else {
	print STDERR ("substr\(\$seq,0,\($startBase+5\)\)\n") if $debug == 1;
	return substr($seq,0,($startBase+5));
      }
    }
    ##update %nuc
    $nuc{$seq[$i]}--;           ##removes the first character
    $nuc{$seq[$i+29]}++;        ##adds new nucleotide
  }
  if ($inRun == 1) {            ##stretch goes to end of sequence
    ##truncate from startBase
    print STDERR ("substr\(\$seq,0,\($startBase+5\)\)\n") if $debug == 1;
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
