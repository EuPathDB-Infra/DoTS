package DoTS::Gene::ESTOrientationPlot;

# -----------------------------------------------------------
# plotting density of 5' ends of 5' ESTs or 3' ends of 3' ESTs
# along assemblies or the genome in the case of aligned gene
#
# Yongchang Gan, Jan 24, 03
# ---------------------------------------------------------

use strict;

my $CN = 'ESTOrientationPlot';

# constructor
#
sub new {
    my($class, $args) = @_;

    my $TAG = $CN . "::new:";

    my $type = lc ($args->{record}->getType);
    die "$TAG type must of DT or DG" unless $type eq 'dt' or $type eq 'dg';
    my $bincount = $args->{bincount};
    $bincount = 10 unless $bincount;

    my $self = { type => $type,
		 record => $args->{record},
		 bincount => $bincount,
		 debug => $args->{debug}
	     };

    # cache the plot to avoid repeated processing
    $self->{plot} = undef;
    $self->{id} = $args->{id};
    $self->{id} = $self->{record}->getId unless $self->{id};
    $self->{strand} = $args->{strand};

    bless $self, $class;
    return $self;
}

sub getPlot {
    my $self = shift;

    my $TAG = $CN . "::getPlot:";
    my $dbg = $self->{debug};

    my $plot = $self->{plot};
    return $plot if $plot;

    my $plot = $self->_getRawPlot;
    my @bins;
    my $p53 = $plot->{p5} + $plot->{p3};
    my $raw_bins = $plot->{bins};
    foreach (@$raw_bins) {
	my ($bp5, $bp3) = ($_->{p5}, $_->{p3}); 
	print "$TAG (bp5, bp3)=($bp5, $bp3)\n" if $dbg;
	if ($p53 < 1) {
	    push @bins, [0,0];
	} else {
	    my $p5 = sprintf("%3.3f", $bp5/$p53);
	    my $p3 = sprintf("%3.3f", $bp3/$p53);
	    push @bins, [$p5, $p3];
	}
    }
    $plot->{bins} = \@bins;

    $self->{plot} = $plot;

    return $self->{plot};
}

sub getScore {
    my $self = shift;

    my $TAG = $CN ."::getScore:";
    my $dbg = $self->{debug};
    my $plot = $self->getPlot;
    my $strand = $self->{strand};
    my $bins = $plot->{bins};
    my ($p5, $p3, $px) = ($plot->{p5}, $plot->{p3}, $plot->{px});

    # compute a empirical score

    # score of -5 .. +5 for plot shape fitness
    my $fit_score = &_getPlotShapeFitnessScore($strand, $bins, $dbg);
    print "$TAG calculated the fitness score of the plot as $fit_score\n" if $dbg;

    # score of  0 ..  5 for count of informative ESTs 
    my $inf_score = &_getInformativeEstScore($p5+$p3);
    print "$TAG calculated the informative EST score of the plot as $inf_score\n" if $dbg;
    
    # compound score
    my $score;
    if ($fit_score == -5) {
	$score = -100;    # merge detected!
    } elsif ($fit_score == 5) {
	$score = 100;     # ideal shape detected!
    } else {
	$score = int 4 * $fit_score * $inf_score;
    }

    $score;
}

sub printPlot {
    my ($self, $opts) = @_;

    my $type = uc $self->{type};
    my $id = $self->{id};

    my $plot = $self->getPlot;
    my $score = $self->getScore;
    my $bins = $plot->{bins};
    my ($p5, $p3, $px) = ($plot->{p5}, $plot->{p3}, $plot->{px});

    my $dot = $opts->{dot}; $dot = '*' unless $dot;
    my $top = $opts->{top}; $top = 40 unless $top;
    my $graph = $opts->{graph};

    my $tot = scalar(@$bins);
    my $strand = $self->{strand}; # needs to call after getPlot so that "strand" is set
    print "# $type.$id; ", ($strand ? "($strand); " : ""),
          "ESTs(5,3,\?)=($p5,$p3,$px), score: $score\n";
    if ($graph) {
	printf("Index\t5\'%" . ($top-2) . "s\t3\'%" . ($top-2) . "s\n", '', '');
    } else {
	printf("Index\t 5'  \t 3'  \n");
    }
    for(my $i=0; $i<$tot; $i++) {
	my ($bpp5, $bpp3) = ($bins->[$i]->[0], $bins->[$i]->[1]);
	if ($graph) {
	    my $s5 = &_makeBar($bpp5, $top, $dot);
	    my $s3 = &_makeBar($bpp3, $top, $dot);
	    printf("%5d\t$s5\t$s3\n", $i);
	} else {
	    printf("%5d\t%1.3f\t%1.3f\n", $i, $bpp5, $bpp3);
	}
    }
}

sub _getRawPlot {
    my $self = shift;

    my $TAG = $CN . "::_getRawPlot:";
    my $dbg = $self->{debug};

    my $type = $self->{type};
    my $num_bins = $self->{bincount};
    my $record = $self->{record};

    my ($p5, $p3, $px) = (0, 0, 0);
    my @bins;
    for (my $i=0; $i<$num_bins; $i++) { 
	push @bins, { p5 => 0, p3 => 0, px => 0 };
    }

    # call getEntries first to populate record
    my $entries = $record->getEntries;
    my $id = $record->getId;
    my $strand = $record->getStrand;
    $self->{strand} = $strand unless $self->{strand};

    my $cs = $record->getChromosomeStart;
    my $ce = $record->getChromosomeEnd;
    foreach my $ent (@$entries) {
	my $p = $ent->getPEnd;
	my $gslen = $ent->getGappedSeqLength;
	my $o = $ent->getAssemblyOffset;
	my $alen = $ent->getAssemblyLength;
	my $gap_con = $ent->getGappedConsensus;

        my $sid = $ent->getSequenceId;
        my $asm_id = $ent->getAssemblyId;

        my $msg = ($type eq 'dt' ? 'DT' : 'AG') . ".$asm_id:$sid"; 
	my ($ugs, $uge) = &_estInUngappedAssembly($o, $gslen, $alen, $gap_con, $dbg, $msg);

	if ($p eq '5') { $p5++; } 
	elsif ($p eq '3') { $p3++; }
	else { $px++; }
       
	my $idx;
	if ($type eq 'dt') {
	    $idx = &_determineIndexDT($p, $alen, $ugs, $uge, $num_bins, $dbg);
            print "$TAG found bin index $idx (EST$sid, DT$asm_id)\n" if $dbg;
	} else {
            my $bs = $ent->getBlocksizes;
            my $qs = $ent->getQstarts;
            my $ts = $ent->getTstarts;
            my $blat_id = $ent->getBLATAlignmentId;
            print "$TAG bin index for (EST$sid, DT$asm_id, BLAT$blat_id)...\n" if $dbg;
	    $idx = &_determineIndexAG($p, $alen, $ugs, $uge, $strand, $cs, $ce,
				      $bs, $qs, $ts, $num_bins, $dbg);
	    print "$TAG found bin index: $idx!\n" if $dbg;
	}
	if ($idx >= 0 && $idx < $num_bins) { 
	    my $bin = $bins[$idx];

	    if ($p eq '5') { $bin->{p5}++; } 
	    elsif ($p eq '3') { $bin->{p3}++; }
	    else { $bin->{px}++; }
	}
    }

    { p5 => $p5, p3 => $p3, px => $px, bins => \@bins };
}

########## file scoped sub ##############

sub _makeBar {
    my ($fraction, $tot, $char) = @_;

    my $bar = '';
    my $c = int ($fraction * $tot);
    for(my $i=0; $i<$tot; $i++) {
	$bar .= ($i >= $c ? ' ' : $char);
    }
    return $bar;
}

# what is the start and end of gapped EST on gapped assembly consensus?
sub _estInUngappedAssembly {
    my ($o, $gslen, $alen, $gap_con, $dbg, $msg) = @_;

    my $TAG = $CN . "::_estInUngappedAssembly:";

    print "$TAG (asmOffset, gapSeqLen, asmLen, gapAsmLen) = ($o, $gslen, $alen, "
	. length($gap_con) . ")\n" if $dbg;

    my $g_s = $o;
    my $g_e = $o + $gslen;
    my $gcl = length($gap_con);
    if ($g_e > $gcl) {
        print STDERR "$TAG WARNING: $msg, $o (offset) + $gslen (gaps) > $gcl (gapped consensus length)\n ";
        $g_e = length($gap_con);
    }

    # correct for gaps
    my $gc1 = $g_s > 0 ? &_getGapCount($gap_con, 0, $g_s, '-') : 0;
    my $gc2 = &_getGapCount($gap_con, $g_s, $g_e, '-');
    my $ug_s = $g_s - $gc1;
    my $ug_e = $g_e - $gc1 - $gc2;

    print "$TAG gapped consensus=$gap_con\n" if $dbg > 2;
    print "$TAG (gapCount1, gapCount2, unGapS, unGapE)=($gc1, $gc2, $ug_s, $ug_e)\n" if $dbg;

    return ($ug_s, $ug_e);
}

sub _determineIndexDT {
    my ($p, $alen, $ugs, $uge, $num_bins) = @_;

    # determin what bin collects EST at $ugs - $uge
    my $bin_step = $alen / $num_bins;
    my $idx;
    if ($p eq '5') {
	$idx = int ($ugs/$bin_step);
    } elsif ($p eq '3') {
	$idx = int ($uge/$bin_step);
    } else {
	$idx = int (($ugs + ($uge-$ugs+1)/2)/$bin_step);
    }
    $idx = $num_bins-1 if $idx == $num_bins;

    return $idx;
}

sub _determineIndexAG {
    my ($p, $alen, $ugs, $uge, $strand, $cs, $ce,
	$blocksizes, $qstarts, $tstarts, $num_bins, $dbg) = @_;

    my $TAG = $CN . "::_determineIndexAG:";

    # assembly's genomic alignment info
    my $blocks = scalar(@$blocksizes);
    die "$TAG invalid genomic alignment blocks info!" 
	unless scalar(@$qstarts) == $blocks and scalar(@$tstarts) == $blocks;
    my $qs = $qstarts->[0];
    my $qe = $qstarts->[$blocks-1] + $blocksizes->[$blocks-1];

    # when on '-' strand, BLAT use query coords w.r.t. 3'
    if ($strand eq '-') {
	($ugs, $uge) = ($alen - $uge, $alen - $ugs);
	print "$TAG (-), (asmLen,newUnGapS,newUnGap)=($alen,$ugs,$uge)\n" if $dbg;
    }

    # project EST to aligned portion of query sequence
    my $prj_qs = ($ugs > $qs ? $ugs : $qs);
    my $prj_qe = ($uge < $qe ? $uge : $qe);
    print "$TAG (qs,ugs)=($qs,$ugs)=>prj_qs=$prj_qs;(qe,uge)=($qe,$uge)=>prj_qe=$prj_qe\n" if $dbg;

    # project further to the chromosome sequence
    my ($prj_cs, $prj_ce);
    for (my $i=0; $i<$blocks; $i++) {
	my $q_s = $qstarts->[$i];
	my $q_e = $qstarts->[$i] + $blocksizes->[$i];
	my $q_s_2 = ($i == $blocks ? -1 : $qstarts->[$i+1]);

	my $t_s = $tstarts->[$i];
	my $t_e = $tstarts->[$i] + $blocksizes->[$i];
	my $t_s_2 = ($i == $blocks ? -1 : $tstarts->[$i+1]);

	print "$TAG block $i, [$q_s, $q_e]:[$t_s, $t_e], ", $i+1, ", [$q_s_2]:[$t_s_2]\n" if $dbg;
	
	if ($prj_qs >= $q_s && $prj_qs <= $q_e) {
	    $prj_cs = $t_s + ($prj_qs - $q_s);
	    print "$TAG block $i [$q_s, $q_e] contains prj_qs=$prj_qs\n" if $dbg;
	}
	if ($prj_qs > $q_e && $prj_qs < $q_s_2) {
	    $prj_cs = $t_s_2;
	    print "$TAG gap between block ($i, ", $i+1, ") contains prj_qs=$prj_qs\n" if $dbg;
	}

	if ($prj_qe >= $q_s && $prj_qe <= $q_e) {
	    $prj_ce = $t_s + ($prj_qe - $q_s);
	    print "$TAG block $i [$q_s, $q_e] contains prj_qe=$prj_qe\n" if $dbg;
	}
	if ($prj_qe > $q_e && $prj_qe < $q_s_2) {
	    $prj_ce = $q_e;
	    print "$TAG gap between block ($i, ", $i+1, ") contains prj_qs=$prj_qs\n" if $dbg;
	}
    }

    # determin what bin will collect EST at $prj_cs - $prj_ce
    my $idx;
    my $bin_step = ($ce - $cs) / $num_bins;
    if ($p eq '5') {
	$idx = ($strand eq '+' ? ($prj_cs - $cs)/$bin_step : ($prj_ce - $cs)/$bin_step);
    } elsif ($p eq '3') {
	$idx = ($strand eq '+' ? ($prj_ce - $cs)/$bin_step : ($prj_cs - $cs)/$bin_step);
    } else {
        $idx = ($prj_cs + ($prj_ce - $prj_cs)/2 - $cs)/$bin_step;
    }

    $idx = int $idx; $idx = $num_bins - 1 if $idx == $num_bins;

    unless ($prj_cs && $prj_ce) {
	$idx = -1;
	print "$TAG could not project $p EST in ${alen}bp assembly ($ugs, $uge) to "
	    . "genome ($strand:$cs-$ce). Alignment info:\n"
	    . "blocksizes: " . join(', ', @$blocksizes) . "\n" 
	    . "qstarts: " . join(', ', @$qstarts) . "\n" 
	    . "tstarts: " . join(', ', @$tstarts) . "\n" if $dbg;
    }

    return $idx;
}

sub _getGapCount {
    my ($gapped_seq, $s, $e, $gap_char) = @_;

    my $sub_seq = substr($gapped_seq, $s, $e - $s);
    my $res = ($sub_seq =~ s/$gap_char{1}//g);

    return $res;
}

sub _getPlotShapeFitnessScore {
    my ($strand, $bins, $dbg) = @_;

    my $TAG = $CN ."::_getPlotShapeFitnessScore:";

    # special case: small number of bins
    my $num_bins = scalar(@$bins);
    if ($num_bins < 1) {
	return -1;
    } elsif ($num_bins < 3) {
	return 0;
    }

    my $score = 0;
    my ($ppp5s, $ppp3s) = &_findBinPeaks($bins);

    print $TAG . " p5 bins: ", join('; ', map { $_->[0] } @$bins),  "\n" if $dbg;
    print $TAG . " p5 peaks: ",
          join("; ", map { $_->[0] . ',' . $_->[1] } @$ppp5s), "\n" if $dbg;
    print $TAG . " p3 bins: ", join('; ', map { $_->[1] } @$bins),  "\n" if $dbg;
    print $TAG . " p3 peaks: ",
          join("; ", map { $_->[0] . ',' . $_->[1] } @$ppp3s), "\n" if $dbg;

    my $p5pks = scalar(@$ppp5s);
    my $p3pks = scalar(@$ppp3s);

    my $i0 = int $num_bins/3;
    my $i1 = int $num_bins * 2/3;
    my ($sum51, $sum52, $sum53, $sum31, $sum32, $sum33) = (0, 0, 0, 0, 0, 0);
    for (my $i = 0; $i < $num_bins; $i++) {
	my ($p5, $p3) = ($bins->[$i]->[0], $bins->[$i]->[1]);
	if ($i <= $i0) {
	    $sum51 += $p5; $sum31 += $p3;
	} elsif ($i >= $i1) {
	    $sum53 += $p5; $sum33 += $p3;
	} else {
	    $sum52 += $p5; $sum32 += $p3;
	}
    }
    my $sum5 = $sum51 + $sum52 + $sum53;
    my $sum3 = $sum31 + $sum32 + $sum33;
    my $hasMerge = &_detectMerges($ppp5s, $ppp3s, $sum5, 0.2, $sum3, 0.25, $dbg);

    if ($hasMerge)
    {
	print "$TAG detected merge\n" if $dbg;
	$score = -5;
    }
    elsif ($p5pks == 1 && $p3pks == 1 &&
	   (($strand eq '+' && $ppp5s->[0]->[0] < $i0 && $ppp3s->[0]->[0] >= $num_bins-2) || 
	    ($strand eq '-' && $ppp3s->[0]->[0] <= 1 && $ppp5s->[0]->[0] > $i1)))
    {
	print "$TAG detected ideal shape\n" if $dbg;
	$score = 5;
    } else {
	if (($strand eq '+' && $sum51 > 0.8*$sum5 && $sum33 > 0.9*$sum3) ||
	    ($strand eq '-' && $sum53 > 0.8*$sum5 && $sum31 > 0.9*$sum3))
	{
	    $score = 4;
	}
	elsif (($strand eq '-' && $sum51 > 0.8*$sum5 && $sum33 > 0.9*$sum3) ||
	       ($strand eq '+' && $sum53 > 0.8*$sum5 && $sum31 > 0.9*$sum3))
	{
	    $score = -4;
	}
	elsif (($strand eq '+' && $sum51 > 0.65*$sum5 && $sum33 > 0.8*$sum3) ||
	       ($strand eq '-' && $sum53 > 0.65*$sum5 && $sum31 > 0.8*$sum3) ||
	       ($p5pks > 0 && 
		(($strand eq '+' && $ppp5s->[0]->[0] < $i0 && $sum33 > 0.9*$sum3) ||
		 ($strand eq '-' && $ppp5s->[$p5pks-1]->[0] > $i1 && $sum31 > 0.9*$sum3))))
        {
            $score = 3;
        }
        elsif (($strand eq '-' && $sum51 > 0.65*$sum5 && $sum33 > 0.8*$sum3) ||
               ($strand eq '+' && $sum53 > 0.65*$sum5 && $sum31 > 0.8*$sum3) ||
	       ($p5pks > 0 &&
		(($strand eq '-' && $ppp5s->[0]->[0] < $i0 && $sum33 > 0.9*$sum3) ||
		 ($strand eq '+' && $ppp5s->[$p5pks-1]->[0] > $i1 && $sum31 > 0.9*$sum3))))
        {
            $score = -3;
        }
        elsif (($strand eq '+' && (($sum51 > 0.8*$sum5 && $sum3 == 0) ||
				   ($sum33 > 0.9*$sum3 && $sum5 == 0))) ||
               ($strand eq '-' && (($sum53 > 0.8*$sum5 && $sum3 == 0) ||
				   ($sum31 > 0.9*$sum3 && $sum5 == 0))))
        {
            $score = 2;
        }
        elsif (($strand eq '-' && (($sum51 > 0.8*$sum5 && $sum3 == 0) ||
				   ($sum33 > 0.9*$sum3 && $sum5 == 0))) ||
               ($strand eq '+' && (($sum53 > 0.8*$sum5 && $sum3 == 0) ||
				   ($sum31 > 0.9*$sum3 && $sum5 == 0))))
        {
            $score = -2;
        }
	elsif (($sum5 < 0.1*$sum3 &&
		$sum31 > 0.2*$sum3 && $sum33 > 0.2*$sum3 && $sum32 < 0.2 * $sum3) ||
	       ($sum3 < 0.1*$sum5 &&
		$sum51 > 0.2*$sum5 && $sum53 > 0.2*$sum5 && $sum52 < 0.2 * $sum5))
	{
	    $score = -1;
	}

    }
    
    $score;
}

sub _findBinPeaks {
    my ($bins) = @_;

    my $tot = scalar(@$bins);
    die "expecting at least 3 bins to find a peak\n" unless $tot >= 3;

    my @padded_bins = @$bins;
    unshift @padded_bins, [0,0]; push @padded_bins, [0,0];

    my (@ppp5s, @ppp3s);
    for (my $i=0; $i<$tot; $i++) {
	my $bin1 = $padded_bins[$i];
	my $bin2 = $padded_bins[$i+1];
	my $bin3 = $padded_bins[$i+2];
	push @ppp5s, [$i, $bin2->[0]]
	    if $bin1->[0] <= $bin2->[0] && $bin2->[0] > $bin3->[0];
	push @ppp3s, [$i, $bin2->[1]]
	    if $bin1->[1] <= $bin2->[1] && $bin2->[1] > $bin3->[1];
    }

    (\@ppp5s, \@ppp3s);
}

# detect existence of 5'(3') peak between two 3'(5') peak
# or a 5'/3' peak between a pair of 5'-3' peak pair
#
sub _detectMerges {
    my ($ppp5s, $ppp3s, $sum5, $threshold5, $sum3, $threshold3, $dbg) = @_;

    my $TAG = $CN ."::_detectMerges:";

    my $c5 = scalar(@$ppp5s);
    my $c3 = scalar(@$ppp3s);
    my $sum5t = $sum5 * $threshold5;
    my $sum3t = $sum3 * $threshold3;

    my @hi_peak5s = ();
    foreach (@$ppp5s) { push @hi_peak5s, $_ if $_->[1] >= $sum5t; }
    my @hi_peak3s = ();
    foreach (@$ppp3s) { push @hi_peak3s, $_ if $_->[1] >= $sum3t; }

    # 5'-5'-5' merge
    my $c5 = scalar(@hi_peak5s);
    if ($c5 >= 3) {
        print "$TAG detected 5'-5'-5' merge\n" if $dbg;
        return 1;
    }
    # 3'-3'-3' merge
    my $c3 = scalar(@hi_peak3s);
    if ($c3 >= 3) {
	print "$TAG detected 3'-3'-3' merge\n" if $dbg;
	return 1;
    }

    for(my $i=0; $i<$c5; $i++) {
	my $i51 = $hi_peak5s[$i]->[0];
	my $i52 = -1;
	$i52 = $hi_peak5s[$i+1]->[0] if $i < $c5-1;
	for (my $j=0; $j<$c3; $j++) {
	    my $p = $hi_peak3s[$j];
	    # 5'-3'-5'
	    if ($i52 > 0 && $p->[0] > $i51+1 &&
		$p->[0] < $i52-1 && $p->[1] >= $sum3t && $p->[1] > 0) 
	    {
		print "$TAG detected 5'-3'-5' merge\n" if $dbg;
		return 1;
	    }
	    # 5'-5'-3'
	    elsif ($i52 > 0 && $p->[0] > $i52 + 1 && $i52 - $i51 > 1 &&
		   $p->[1] >= $sum3t && $p->[1] > 0) 
	    {
		print "$TAG detected 5'-5'-3' merge\n" if $dbg;
		return 1;
	    } 
	    if ($j < $c3-1) {
		# 5'-3'-3'
		my $p1 = $hi_peak3s[$j+1];
		if ($p->[0] > $i51 + 1 && $p1->[0] - $p->[0] > 1 &&
		    $p->[1] >= $sum3t && $p->[1] > 0 &&
		    $p1->[1] >= $sum3t && $p1->[1] > 0) {
		    print "$TAG detected 5'-3'-3' merge\n" if $dbg;
		    return 1;
		}
	    }
	}
    }
    
    for(my $i=0; $i<$c3; $i++) {
	my $i31 = $hi_peak3s[$i]->[0];
	my $i32 = -1;
	$i32 = $hi_peak3s[$i+1]->[0] if $i < $c3-1;
	for (my $j=0; $j<$c5; $j++) {
	    my $p = $hi_peak5s[$j];
	    # 3'-5'-3'
	    if ($i32 > 0 && $p->[0] > $i31 + 1 && 
		$p->[0] < $i32 - 1 && $p->[1] >= $sum5t && $p->[1] > 0)
	    {
		print "$TAG detected 3'-5'-3' merge\n" if $dbg;
		return 1;
	    }
	    # 3'-3'-5'
	    elsif ($i32 > 0 && $p->[0] > $i32 + 1 && $i32 - $i31 > 1 &&
		   $p->[1] >= $sum5t && $p->[1] > 0)
	    {
		print "$TAG detected 3'-3'-5' merge\n" if $dbg;
		return 1;
	    }
	    if ($j < $c5-1) {
		# 3'-5'-5'
		my $p1 = $hi_peak5s[$j+1];
		if ($p->[0] > $i31 + 1 && $p1->[0] - $p->[0] > 1 &&
		    $p->[1] >= $sum5t && $p->[1] > 0 &&
		    $p1->[1] >= $sum5t && $p1->[1] > 0) {
		    print "$TAG detected 3'-5'-5' merge\n" if $dbg;
		    return 1;
		}
	    }
	}
    }

    print "$TAG detected no merge\n" if $dbg;
    return 0;
}

sub _getInformativeEstScore {
    my $tot = shift;
    
    my $score;
    if ($tot >= 100) {
	$score = 5;
    } elsif ($tot >= 50) {
	$score = 4;
    } elsif ($tot >= 25) {
	$score = 3;
    } elsif ($tot >= 2) {
	$score = 2;
    } elsif ($tot > 0) {
	$score = 1;
    } else {
	$score = 0;
    }

    $score;
}

1;
