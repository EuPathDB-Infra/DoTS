package DoTS::DotsBuild::Plugin::FrameFinder;

# Run framefinder on Assemblies, storing the results in the database.
#
# Created:  05.06.2001
# Last edited: 10.22.2001
# 
# ----------------------------------------------------------

@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use GUS::Model::DoTS::Assembly;
use GUS::Model::Core::Algorithm;
use GUS::Model::DoTS::TranslatedAAFeatSeg;
use GUS::Common::Sequence;
use CBIL::Bio::SequenceUtils;

# ----------------------------------------------------------
# GUSApplication
# ----------------------------------------------------------

sub new {
  my ($class) = @_;

  my $self = {};
  bless($self, $class);

  my $usage = 'Plug-in for reconstruction of ORF by framefinder on assembly sequences and record the results in the TranslatedAAFeature and TranslatedAASequence tables';

  my $easycsp =
    [
     {o => 'restart',
      t => 'string',
      h => 'restarts from last entry in TranslatedAASequence....
                             takes list of row_alg_invocation_ids "234, 235"!',
     },
     {o => 'idSQL',
      t => 'string',
      h => 'SQL statement:  must return list of primary identifiers (na_sequence_id) from --table_name',
     },
     {o => 'wordfile',
      t => 'string',
      h => 'word probability file for framefinder',
      e => [ qw ( hum_GB123.wordprob mouse_GB123.wordprob ) ],
     },
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'ffdir',
      t => 'string',
      h => 'directory for framefinder_GUS location',
     },
     {o => 'dianadir',
      t => 'string',
      h => 'directory for diana program location',
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
my $debug = 0;

sub run {
  my $i = 0;
  my $M = shift;
  my $ctx = shift;
  print STDERR "framefinder: COMMIT ", $ctx->{cla}->{'commit'} ? "****ON****" : "OFF", "\n";
  print STDERR "Establishing dbi login\n";

  my $dbh = $ctx->{'self_inv'}->getDatabase()->getQueryHandle();                
  my $time1 = scalar localtime;
  my $framefinderdir=$ctx->{cla}->{ffdir} if ($ctx->{cla}->{ffdir});
  my $dianadir=$ctx->{cla}->{dianadir} if ($ctx->{cla}->{dianadir});
  my $Framefinder = $framefinderdir.'/bin/framefinder_GUS'; # location of FrameFinder
  my $Diana = $ctx->{cla}->{dianadir}.'/atg';
  if (!(-e $Diana))
	{die "Framefinder: No diana program at the site mentioned\n";}
  if (!(-e $Framefinder))
        {die "Framefinder: No Framefinder program at the site mentioned\n";}


  my %ignore;                   # skipping entries already processed
  ## want to be able to ignore entries already done!!
  # current key is non-zero tranlation score: if non-zero, then processed.
  if ($ctx->{cla}->{'restart'}) {
    my $query = 
"select rf.na_sequence_id 
from dots.rnafeature rf, dots.translatedaafeature tf  
where rf.na_feature_id = tf.na_feature_id 
and tf.row_alg_invocation_id in ($ctx->{cla}->{'restart'})";
#    my $query = "select distinct r.na_sequence_id from rnasequence r, translatedaafeatute tf, assembly a where r.na_feature_id = tf.na_feature_id and r.na_sequence_id = a.na_sequence_id and a.description != 'DELETED' and tf.translation_score is not null";
    print STDERR "Restarting: Querying for the ids to ignore\n$query\n";
    my $stmt = $dbh->prepare($query);
    $stmt->execute();
    while ( my($id) = $stmt->fetchrow_array()) {
      $ignore{$id} = 1;
    }
    print STDERR "Ignoring ".scalar(keys%ignore)." entries\n";
  }

  #         Starting the run
  my $verbose; # = $ctx->{cla}->{verbose};
  die "Error: idSQL query string parameter should not be empty\n" if (not defined($ctx->{cla}->{idSQL}));
  die "Error: cannot find wordfile $ctx->{cla}->{wordfile} for framefinder" unless (-e "$framefinderdir/wordProb/".$ctx->{cla}->{wordfile});
  print STDERR "$ctx->{cla}->{idSQL}\n"; 
  my $stmt = $dbh->prepare($ctx->{cla}->{idSQL});
  $stmt->execute();
  my @todo;
  my $cte = 0;
  while ((my($nas) = $stmt->fetchrow_array())) {
    $cte++;
    if ($ctx->{cla}->{testnumber} && $cte > $ctx->{cla}->{testnumber}) {
      $stmt->cancel(); last;
    }
    print STDERR "Retrieving entry $cte\n" if($cte % 10000 == 0);
    push(@todo,$nas) unless exists $ignore{$nas};
  }
  my $totalToDo = scalar(@todo);
  print STDERR "$totalToDo entries remaining to process.\n";
  undef %ignore;                ##free this memory...

######
# reading from the file
# open(F,"bad_ff_ids.acc");
# while (my $ss = <F>) {chomp($ss); push (@todo,$ss);}
# my $totalToDo = scalar(@todo);
# print STDERR "$totalToDo entries remaining to process.\n";


  # Select word probability file according to taxon id
  #
#  my $wordprob='embl59.wordprob'; #for human and mouse.
#  my $wordprob='hum_GB_123.wordprob'; #for human


  # creating single algorithm object;
  my $alg = GUS::Model::Core::Algorithm->new({'name'=>'FrameFinder'});
  $alg->retrieveFromDB();
  my $alg_id = $alg->get('algorithm_id');
  $alg =  GUS::Model::Core::Algorithm->new({'name'=>'TrivialTrans'});
  $alg->retrieveFromDB();
  my $triv_alg_id = $alg->get('algorithm_id');
  my $countEntries = 0;
  my $triv_count = 0;

  foreach my $naSeqId (@todo) {
    my $naSeq = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $naSeqId});
    if (!$naSeq->retrieveFromDB()) {
      print STDERR "ERROR: unable to retrieve assembly for $naSeqId\n";
      next;
    }
    my $seq = $naSeq->getSequence();

    if (not defined($seq)){
      print STDERR "Unable to retrieve sequence with na_sequence_id = $naSeqId\n";
      next;
    }

    # Write the sequence to a temporary file
    #
    my $tmpFile = "$$.fr.tmp";
    open(TFILE, ">$tmpFile");
    print TFILE &GUS::Common::Sequence::toFasta($seq, $naSeq->get('description'), 80), "\n";
    close(TFILE);
    print "\nProcessing the sequence ".$naSeqId."\n" if $verbose;

    #                        print STDERR "Created $tmpFile\n";

    # Summary data for Framefinder output
    #
    my $parameterset= "-O 14 -F -17 -w $framefinderdir"."wordProb/$ctx->{cla}->{wordfile}";
    my @F = `$Framefinder $parameterset $tmpFile`;
    #                        print STDERR "Running: @F\n";
    my $start_pos;
    my $i;
    my $end_pos;
    my $score;
    my $strand;
    my $trans_model;

    # new things for segments
    my $numofsegs;
    my @breakpoints;
    my @segments;
    my @numofnucs;
    my @shift_type;
    my @scores;
    my $flag = 0;
    my $contt = 0;
    my $cont2 = 0;
			
    # Rtun framefinder
    foreach my $line (@F) {
      if ($line =~/number of segments\s+(\d+)/) {
        $numofsegs = $1;
      }                        
      if ($line =~/(\S+)\s+after pos\s+(\d+)\s+for\s+(\d+)(\s+)nucleotide\(s\)(\s+)score\s+(\d+\.\d+)(\s+)/) 
	  {
	@shift_type[$contt] = $1;
        @breakpoints[$contt] = $2; 
        @numofnucs[$contt] = $3;
        @scores[$contt]=$6;

#        if ($line =~/deletion/) {
#          @shift_type[$contt]='deletion';
#        } else {
#          @shift_type[$contt]='insertion';
#        }
        $contt++;
      }

      if ($line =~/\[framefinder\s+\((\d+)\,(\d+)\)\s+score=(\d+\.\d+).*\{(\w+)\,(\w+)\}.*$/) {
        $start_pos = $1;
        $end_pos = $2;
        $score = $3;
        $strand = $4;
        $trans_model = $5;
        $flag= 1;
	if ($strand=~/revcomp/) {
	$strand=1;
      } elsif ($strand=~/forward/) {
        $strand=0;
      }
        next;
      }
      if ($flag ) {
        @segments[$cont2] = $line; $cont2++;
      }
    }
    ;
   
    #                        print STDERR "$start_pos $end_pos $score $strand \n";

    #                              DEBUG INFORMATION                
    #print STDERR " parameters $numofsegs, $start_pos, $end_pos $contt , $cont2 \n";

    #              foreach $flag(@breakpoints) {print "$flag\n";}
    #              foreach $flag(@scores) {print "$flag\n";}
    #               foreach $flag(@segments) {print "$flag\n";}
    #		foreach $flag(@numofnucs) {print "$flag\n";}
    #		foreach $flag(@shift_type) {print "$flag\n";}
    #                             end debug
    my $aa_seq = "";
    my $first_seg = "";  # first aa segment
    my $num_start_x=0;
    my $num_end_x=0;
    my $in_beg = 1;
    foreach my $lin (@segments) {
      if ($lin =~ /^(X*)\n/)    # only XXXX or empty string
        {# not now: should adjust for ins/del: 
#	if ($in_beg) {$num_start_x +=length($1)-1;}
        } else {
          chomp($lin); $aa_seq .= $lin; 
	$in_beg = 0;   # starting Xs ended
	if (length($first_seg)==0) {$first_seg=$lin; $first_seg =~ s/^(X+)//;
#	print " first segment is $first_seg  and the first Xs number are ". length($1)."\n";
				}
        }
    }
	if ($aa_seq =~ /^(X+)/) {$num_start_x+=length($1);}
	$aa_seq =~ s/^(X+)//;
        if ($aa_seq =~ /(X+)$/) {$num_end_x=length($1);}
	$aa_seq =~ s/(X+)$//;

#   print "Number of xs in front of the sequence is $num_start_x behind is $num_end_x\n";
 
    my $clev = CalcPvalue($seq,$aa_seq);
#################    Trying to check start position in na_sequence
 if ($debug)
 {
my $aaseg = substr($first_seg,0,length($first_seg) < 13 ? length($first_seg) - 3 : 10);
if ($end_pos>length($seq))
   {$end_pos = length($seq);}
    if (!$start_pos) {$start_pos = -2;}
my $stpos = $start_pos -2 + $num_start_x*3;    # here is -3 
if ($stpos<=0)
   {$stpos=1;}   
   my $o;
   for($o=0; $o<4; $o++) {
	my $seq1 = $seq;
	if ($strand) {$seq1 = CBIL::Bio::SequenceUtils::reverseComplementSequence($seq1)};
      my $naseg1 = substr($seq1,$stpos+$o-1,length($aaseg) * 3);
      my $aaseg1 = CBIL::Bio::SequenceUtils::translateSequence($naseg1);
	my $seq2 = substr($seq1, $o, length($seq1));
	my $aaseg2 = CBIL::Bio::SequenceUtils::translateSequence($seq2);
      $aaseg1 =~ s/X/\./g;
      $aaseg1 =~ s/\*/\./g;
      if (not ($aaseg =~ /$aaseg1/))
	{
   print STDERR "START POSITION $stpos failed in frame $o for $naSeqId\n na seq is $naseg1 and aa $aaseg1   aaseq  $aaseg\n" if $verbose;
#      print " in".CBIL::Bio::SequenceUtils::breakSequence($seq1)." translated is $aaseg2\n 
	}
	else
	{$start_pos = $stpos+$o; 
#	print "START POSITION IS OK: $start_pos END position is $end_pos in na_seq in $o frame with strand $strand\n na seq is $aaseg1   aaseq  $aaseg\n"; 
#	print " Here the sequences are ".CBIL::Bio::SequenceUtils::breakSequence($aaseg2)."\n".CBIL::Bio::SequenceUtils::breakSequence($aa_seq)."\n";
	last;}
 } # for
 } # if

    my $updating = 1;           #in the Run mode it should be update at least for now

    my $triv_flag = 0;
    my $triv_aa_length;
    my $triv_seq;
    my $triv_trans_start,
    my $triv_trans_revcomp;
    my $triv_threshold = 1e-16;

    if (($clev > 0.2))      #small proteins:  if trivial is 3 times longer then we take it (not sure if it's needed)
#  ||	(($numofsegs>1) && ($clev<$triv_threshold)))  	     #very long ORF: if trivial is longer - we take it in any case (most probably it comprises mRNA) - this is now abandoned due to Extremely small number of prots corrected;

	{
      #checking against trivial frame
      ($triv_aa_length, $triv_seq, $triv_trans_start, $triv_trans_revcomp) 
        = Trivial($seq);        #caculating the trivial length;
#	print "triv trans options start $triv_trans_start lenght $triv_aa_length  reverse $triv_trans_revcomp\n";
      #			die "Triv_trans made";


#	if ($triv_trans_revcomp) 
#	$triv_trans_start = $triv_trans_start-length($triv_seq)*3+1;

   ##     debug print
#	print "length of trivial sequence is $triv_aa_length framfinder is ".length($aa_seq)."\n"; 
#	print "trivial aa_sequence\n".CBIL::Bio::SequenceUtils::breakSequence($triv_seq);
#	print "source sequence\n".CBIL::Bio::SequenceUtils::breakSequence($seq);

	if (($triv_trans_start+$triv_aa_length*3-1) > length($seq)) {print STDERR "Something WRONG with triv_trans\n"; print STDERR "triv_trans_start".$triv_trans_start." na_sequence length ".length($seq)." aa_seq_length ".$triv_aa_length." trivial is reversed ".$triv_trans_revcomp."\n";}

      if (
	(($triv_aa_length > 3*length($aa_seq)) && ($clev>0.2)) ||
	(($triv_aa_length >= length($aa_seq)) && ($clev<$triv_threshold)))

		 {
        $triv_count++; 
        $clev = CalcPvalue($seq,$triv_seq);
        #	print "framefinder is NOT optimal here - triv trans is longer\n";
        $triv_flag = 1;
      }
    }
    # Set results into table
    #
    my @trSeqs = $naSeq->getTranslatedAASequences(1);
    my $tr_AASeq;  # TranslatedAASequence object
    my @transsegs; # TranslatedAAFeatureSegment objects;
    if (scalar(@trSeqs) > 1) {  # check if there is prediction_algo
      ##want to only use the one that has  correct prediction_algo...
      print STDERR "ERROR:  There is more than one translatedaasequence for $naSeqId\n";
    } elsif (scalar(@trSeqs) ==  1) {
      $tr_AASeq = @trSeqs[0];
    } else {
      ##now need to create these guys...
      print STDERR "ERROR:  There aren't any translatedaasequences for $naSeqId\n";
      next;
    }

    my $translate = $tr_AASeq->getChild('DoTS::TranslatedAAFeature'); #should be already there
    if ($updating)              #for now it is abundant
      {
        foreach my $seg ($translate->getChildren('DoTS::TranslatedAAFeatureSegment',1)) {
          $seg->markDeleted();  #deleting all old segments
        }
      }
    
    $translate->setIsPredicted(1) unless $translate->getIsPredicted() == 1;
    $translate->setReviewedStatusId(0) unless $translate->getReviewedStatusId() == 0;
    $translate->setPValue($clev) unless $translate->getPValue()==$clev;


    if ($triv_flag) {                 #trivial case


      $translate->setIsSimple(1) unless $translate->getIsSimple() == 1;
      $translate->setTranslationModel('') unless $translate->getTranslationModel() == '';
      $translate->setTranslationScore('0') unless $translate->getTranslationScore() == 0;
      $translate->setNumberOfSegments('1') unless $translate->getNumberOfSegments() == 1;
      $translate->setDianaAtgScore('0') unless $translate->getDianaAtgScore() == 0;
      $translate->setDianaAtgPosition('0') unless $translate->getDianaAtgPosition() ==0;
      $translate->setIsReversed($triv_trans_revcomp) unless $translate->getIsReversed() == $triv_trans_revcomp;
      $translate->setPredictionAlgorithmId('3289') unless $translate->getPredictionAlgorithmId() == 3289;
      my $triv_trans_stop = $triv_trans_start+3*length($triv_seq)-1;
     
     my $TransStart = &getPosition(length($seq),$translate->getIsReversed(),$triv_trans_start);
     my $TransStop  = &getPosition(length($seq),$translate->getIsReversed(),$triv_trans_stop ); #naSeq->getLength() sucks??
#	print "Trivial translation start: $TransStart  stop $TransStop\n";
#        print "Trivial translation start: $triv_trans_start length->getlength = ". $naSeq->getLength()." reversed:". $translate->getIsReversed()." length(seq) = ".length($seq)." \n";
        if ($TransStop>$TransStart)
                {
     $translate->setTranslationStart($TransStart);
     $translate->setTranslationStop($TransStop);
        }
        else
                {
     $translate->setTranslationStart($TransStop);
     $translate->setTranslationStop($TransStart);
        }

      print "setting sequence with the trivial translation\n" if $verbose;
#      print CBIL::Bio::SequenceUtils::breakSequence($triv_seq)."\n"; 
      $tr_AASeq->setSequence($triv_seq); # unless $tr_AASeq->getSequence() == $triv_seq;
#      print "set the trivial sequence".CBIL::Bio::SequenceUtils::breakSequence($tr_AASeq->getSequence())."\n";

###  Setting the segment for trivial translation;
     undef @transsegs;
     @transsegs[0]= GUS::Model::DoTS::TranslatedAAFeatSeg->new();
     @transsegs[0]->setTranslationScore(0);
     @transsegs[0]->setAaStartPos(1);
     @transsegs[0]->setAaEndPos(length($triv_seq));
     @transsegs[0]->setNucleotidesShifted(0);
     @transsegs[0]->setTypeOfShift('');
     @transsegs[0]->setStartPos($translate->getTranslationStart());
     @transsegs[0]->setEndPos($translate->getTranslationStop());

####################end non-trivial translation
    } else {                    #non-trivial translation
      my $shift_start = $start_pos<3 ? -1 : 2; #framefinder locates 3-d nucleotide
      $shift_start = $shift_start==0 ? -1: $shift_start;
#      my $shift_start = 2;
      $translate->setIsSimple(0) unless $translate->getIsSimple() == 0;
      $translate->setPredictionAlgorithmId('64') unless $translate->getPredictionAlgorithmId() == 64;
      $translate->setTranslationScore($score) unless $translate->getTranslationScore() == $score;
      $translate->setIsReversed($strand) unless $translate->getIsReversed() == $strand;
      my $TransStart = &getPosition(length($seq),$translate->getIsReversed(),$start_pos);
      $translate->setTranslationStart($TransStart)-$shift_start unless $translate->getTranslationStart() == $TransStart-$shift_start; #framefinder locates 3-d nucleotide
      my $TransStop = &getPosition(length($seq),$translate->getIsReversed(),$end_pos);
      $translate->setTranslationStop($TransStop) unless $translate->getTranslationStop() == $TransStop;      
      $translate->setParameterValues($parameterset) unless $translate->getParameterValues() == $parameterset;
      $translate->setTranslationModel($trans_model) unless $translate->getTranslationModel() == $trans_model;
      
      
      # performing DIANA (atg) in case of strict model
      my ($atg_pos, $atg)  = (0,0);
#      if ($trans_model =~ /strict/) {
##        my @FF = `$dianadir.est $tmpFile`;
        my @FF = `$Diana $tmpFile`;
#        			print STDERR "Running: @FF\n";
        foreach my $line (@FF) {
          if ($line =~ /pos:\s+(\d+)\s+atg_score:\s+(\d+\.\d+)/) {
            $atg_pos = $1; $atg = $2;
            #	print " atgpos $atg_pos atg $atg \n";
          } else {
            print "DIANA: did not find anything\n" if $verbose;
          }
        }
#      } translate in any case
      $translate->setDianaAtgScore($atg);
      $translate->setDianaAtgPosition($atg_pos);
      
      my $currpos =0;
#      my @transsegs;
      # Set children TranslatedAAFeatureSegment objects;
      
      # checking if there are pseudo segments consisting only of XXX or empty string (framefinder artifact)
      my $realsegs = $numofsegs;
      
      if ($numofsegs>1) { 
        foreach my $segm (@segments) {
          if ($segm =~ /^(X*)\n/) # only XXXX or empty string
            {
              $realsegs--;
            }
            
        }
      }
                     	
      $translate->setNumberOfSegments($realsegs) unless $translate->getNumberOfSegments() == $realsegs;

      print "setting sequence with the ff translation\n" if $verbose;
      $tr_AASeq->setSequence($aa_seq); #  unless $tr_AASeq->getSequence() == $aa_seq;

	
      if ($realsegs)          #more than 0 REAL segment - then chreate Tables
        {		
          @breakpoints[$numofsegs-1] = $end_pos-$num_end_x*3; #setting the last breakpoint for cycle
	  @shift_type[$numofsegs-1] = 0; # just in case there is some calculation;
          @scores[$numofsegs-1] = $score; #setting the last score for cycle
          @breakpoints = sort {$a <=> $b} @breakpoints;
 
          my $j = 0;            # this is for realsegs;
		
          for ($i=0; $i<$numofsegs; $i++) {
            if (not (@segments[$i] =~ /^(X*)\n/)) #the real segment is observed that's why j is introduced	
              {
                #			print "number of real segments $realsegs\n";
                chomp(@segments[$i]);
                @transsegs[$j]= GUS::Model::DoTS::TranslatedAAFeatSeg->new();

                @transsegs[$j]->setTranslationScore($scores[$i+1]-$scores[$i]);
                @transsegs[$j]->setAaStartPos($currpos+1);

#		print " This is for the parse: aa start pos is $currpos ";
                $currpos=$currpos+length(@segments[$i]);
		if ($j==0) {$currpos -= $num_start_x;} #shifting the number of X-s deleted in the beginning;
#		print "And aa end post is $currpos\n";
                @transsegs[$j]->setAaEndPos($currpos);
		
                if ($j>0) {
		my $grad;                      
		if (@shift_type[$i-1] =~ /inser/) {$grad = -1;} else {$grad = 1};
		print STDERR "Type of shift is ".@shift_type[$i-1]." so the grad is $grad\n" if $debug;
		# $grad = 0; # there is no need to do this correction;
		my $StartPos = @breakpoints[$i-1];
										# was -5
                if ($StartPos<=0) {$StartPos=1;}
                print STDERR "Initial Start is $StartPos   length is".$naSeq->getLength()."\n" if $debug;
		$StartPos = &getPosition(length($seq),$translate->getIsReversed(),$StartPos);
                  @transsegs[$j]->setTranslationScore($scores[$i]-$scores[$i-1]);
                  @transsegs[$j]->setStartPos($StartPos);
                } else {
                  @transsegs[$j]->setTranslationScore($scores[$i]);
		  if ($start_pos==0) {$start_pos = 1;}
                  @transsegs[$j]->setStartPos(&getPosition(length($seq),$translate->getIsReversed(),$start_pos));  #taking the initial start position
        	print "Start is $start_pos while setting the segs  length is".$naSeq->getLength()."\n" if $verbose;
                }               # starting from the 1-st position
		
		# print "EndPosition shift is ";
		my $grad = (@shift_type =~ /insert/ ? 1:-1);
		my $EndPos = &getPosition(length($seq),$translate->getIsReversed(),@breakpoints[$i]+$grad*@numofnucs[$i]);
		if ($EndPos<=0) {$EndPos = 1;}
                @transsegs[$j]->setEndPos($EndPos) unless @transsegs[$j]->getEndPos() == $EndPos;       
                @transsegs[$j]->setNucleotidesShifted(@numofnucs[$i]) 
		unless @transsegs[$j]->getNucleotidesShifted() == @numofnucs[$i];
                @transsegs[$j]->setTypeOfShift(@shift_type[$i]) unless @transsegs[$j]->getTypeOfShift() == @shift_type[$i];
#		print "Start Position is ". @transsegs[$j]->getStartPos()."  End Position is ".@transsegs[$j]->getEndPos()."Length of sequence".length($seq)."\n";
                $j++;
              }  #if not (@segments[$i...

          } #for
         if ($realsegs != $j) {
                 print "Error: too much real sequences $realsegs  $j \n";
             } 
} # if realsegs
else {
      print STDERR "ERROR: There aren't any segments for $naSeqId\n";
         #      next;
             }

    }                           #else non-trivial
	foreach my $seg(@transsegs)
	{
  print "\nCorrecting the segments\nStart Position is".$seg->getStartPos()."  End Position is ".$seg->getEndPos()."\n"if $debug;
                      ### getting the reverse complement in accordance to GUS denotations
	if ($seg->getStartPos()>$seg->getEndPos())
		{
		my $tmp = $seg->getStartPos();
		$seg->setStartPos($seg->getEndPos());
		$seg->setEndPos($tmp);
	}

		### Here is the module adjusting the positions
		my $naseg = substr($seq, $seg->getStartPos()-1, $seg->getEndPos()-$seg->getStartPos()+1);

#    		my $naseg = $seg->getNASequenceSegment();
		if ($translate->getIsReversed())
			{$naseg = CBIL::Bio::SequenceUtils::reverseComplementSequence($naseg);}
		my $aaseg = substr($tr_AASeq->getSequence(), $seg->getAaStartPos()-1, $seg->getAaEndPos() - $seg->getAaStartPos()+1);
		my $aaseg = substr($aaseg,0,length($aaseg) <13 ? length($aaseg) - 3 : 10);
		for (my $o =0; $o <7; $o++) {
			my $naseg1 = substr($naseg, $o, length($aaseg) * 3);
		my $aaseg1 = CBIL::Bio::SequenceUtils::translateSequence($naseg1);
        	$aaseg1 =~ s/X/\./g;
        	$aaseg1 =~ s/\*/\./g;
		if (not($aaseg =~ /$aaseg1/)) {
#        print STDERR "Correction failed in frame $o for $naSeqId\n na seq: $aaseg1 aaseg: $aaseg\n reversed is". $translate->getIsReversed()."\n";
		}
		else {
#			print "It is fine in frame $o\n";
			if ($translate->getIsReversed())
				{$seg->setEndPos($seg->getEndPos()-$o);}
			else
				{$seg->setStartPos($seg->getStartPos()+$o);}
#        print STDERR "Correction is OK for $naSeqId in frame $o\n";
#	print " na seq: $aaseg1  aaseq:  $aaseg\n reversed is". $translate->getIsReversed()."\n";
			last;
		}
		} #for	

	}
    
    $translate->addChildren(@transsegs);
    $tr_AASeq->submit();
#    print " The aa sequences for trivial after submission and initial\n ".CBIL::Bio::SequenceUtils::breakSequence($tr_AASeq->getSequence())."\n".CBIL::Bio::SequenceUtils::breakSequence($triv_seq)."\n";

    ##put in the sanity check here to make certain that the locations are correct
    my @failed = 0;
    my $i = 0;
    print "Number of segments is ".$translate->getChildren('DoTS::TranslatedAAFeatureSegment')."\n" if $verbose;
    foreach my $seg(@transsegs) {
#           ($translate->getChildren('TranslatedAAFeatureSegment')){
    
#   print " NA and AA positions ".$seg->getStartPos()."  ".$seg->getEndPos()." ".$seg->getAaStartPos()."  ".$seg->getAaEndPos()."\n";
    my $naseg = $seg->getNASequenceSegment(); 
#	print length($naseg)."\n";
#	print "NA sequence" .CBIL::Bio::SequenceUtils::breakSequence($naseg);

    my $aaseg0 = $seg->getAASequenceSegmentFromTranslatedAASequence(); # - have not submitted - other seq;
   
#        print "AA sequence" .CBIL::Bio::SequenceUtils::breakSequence($aaseg);
   
    my $aaseg = substr($aaseg0,0,length($aaseg0) < 13 ? length($aaseg0) - 3 : 10);
    if (length($aaseg) < 3) {print "Sanity check is abandoned due to short segment\n" if $verbose; next;}
    for(my $o=0; $o<7; $o++) {
        my $naseg1 = substr($naseg,$o,length($aaseg) * 3);
	my $aaseg1 = CBIL::Bio::SequenceUtils::translateSequence($naseg1);
	$aaseg1 =~ s/X/\./g;
        $aaseg1 =~ s/\*/\./g;	
	$aaseg =~ s/\*/\./g;
#        next if $naseg1 =~ /N/;
        if (length($naseg1)<length($aaseg)*3) {$aaseg = substr($aaseg,0,length($naseg1)/3);}
        if (not($aaseg =~ /$aaseg1/)){
            print STDERR "Sanity Check failed in frame $o for $naSeqId\n na seq: $aaseg1  aaseq:  $aaseg\n reversed is ". $translate->getIsReversed()."\n" if $debug;
            @failed[$i] = 1;
        } #if
	else
	{
#	print "It is fine - locations are correct for frame $o\n";
        print STDERR "Sanity Check is OK for $naSeqId\n na seq: ".CBIL::Bio::SequenceUtils::translateSequence($naseg1)." aaseq:  $aaseg\n ".substr($aaseg0, -4,4)." reversed is". $translate->getIsReversed()."\n" if $verbose;
        my $aafStartPos = $o;
        my $goodFirst = 1;
	@failed[$i]=0;
	last;
	} #else
    } #for
    $i++;
    }
#    if($failed)
{  ##failed sanity check...
        print" Overall test for this sequence: " if $verbose;
        for (my $i=0; $i<5;$i++) {print " @failed[$i]" if $verbose;}
        print "\n" if $verbose;
      ##what to do
}
      # Set children
      #

#    $translate->addChildren(@transsegs);
    $tr_AASeq->submit();
    # Remove temporary file(s)
    # 
    unlink $tmpFile;
    $countEntries++;
    print STDERR "$countEntries processed ",($totalToDo - $countEntries)," remaining\n" if $countEntries % 10 == 0;
    $ctx->{ self_inv }->undefPointerCache();
  } 
  #                       printf "run finished, processed $countEntries entries\n";
  print STDERR "Triv vs total number of entries: $triv_count vs $countEntries\n";
  my $time2 = scalar localtime;
  print STDERR "start time: $time1  End time: $time2\n";
  my $results = "run finished, processed $countEntries.";
  return $results;
}                               #Run

sub CalcPvalue {
  # calculation of the P-value based on Poisson approximation of geometric distribution;
  my ($source, $target) = @_;
  my $target_len = length($target);
  my $Npep = 2*length($source) - 6*$target_len+6; #it's scanning length in 6-frames;
  if ($Npep <=0)                #it happens due to the insertions/deletions by framefinder
    {
      $Npep = 1;
    }
  my $l_lambda = log(3/64)+$target_len*log(61/64)+log($Npep);
  my $lambda = exp($l_lambda);
  my $confid_level = 1-exp(-$lambda);
  # print "P-value is $confid_level\n";
  return $confid_level;
}


sub getPosition {
  my($len,$isrev,$pos) = @_;
  if($isrev){
    return $len - $pos + 1;
  }else{
    return $pos;
  }
}
1;

