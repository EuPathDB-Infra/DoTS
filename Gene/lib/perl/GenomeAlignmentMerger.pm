package DoTS::Gene::GenomeAlignmentMerger;

use strict;
use DoTS::Gene::Util;
use DoTS::Gene::GenomeFeature;
use DoTS::Gene::CompositeGenomeFeature;
use CBIL::Util::Disp;

# constructor
#
sub new {
    my($class, $sorted_alignments, $merge_criteria, $dbh) = @_;

    my $om = $merge_criteria->{overlap_merge};
    my $cm = $merge_criteria->{clonelink_merge};
    my $epc = $merge_criteria->{est_pair_cache};
    my $pm = $merge_criteria->{proximity_merge};
    die "dbh/EstPairCache not provided" if $cm && !($epc && $dbh);
    my $self = { sa => $sorted_alignments, dbh => $dbh,
		 om => $om, cm => $cm, epc => $epc, pm => $pm };

    $self->{cachePlus} = [];
    $self->{cacheMinus} = [];

    bless $self, $class;
    return $self;
}

########## subroutine ###################################################

sub setVerboseLevel { $_[0]->{debug} = $_[1]; }
sub getVerboseLevel { $_[0]->{debug}; }

# get merged alignments in this genomic region,
# return array ref of DoTS::Gene::CompositeGenomeFeature
#
sub getCompositeGenomeFeatures {
    my $self = shift;

    my $res = $self->{merged};
    return $res if defined $res;

    print "seed genome features with alignments\n";
    $self->_seed;

    if ($self->{om}) {
	print "merge transitively alignments with span overlap\n";
	$self->_doOverlapMerge;
    }

    if ($self->{cm} > 0) {
	print "merge transitively alignments within certain distance if linked by est pairs\n";
	$self->_doClonelinkMerge;
    }

    if ($self->{pm} > 0) {
	print "merge transitively alignments within certain distance\n";
	$self->_doProximityMerge;
    }

    print "combine merged groups into CompositeGenomeFeatures\n";
    $res = $self->_pack;

    $self->{merged} = $res;
    return $self->{merged};
}

#####################################################################################

sub _seed {
    my $self = shift;

    my $srt_aln = $self->{sa};
    foreach my $aln (@$srt_aln) {
	my $id = $aln->getBlatAlignmentId();
	my $q_seq_id = $aln->getQueryNaSequenceId();
	my @coords = $aln->getTargetCoordinates();

	my $seed = { id => $id, coords => \@coords, alns => [$aln], qseqs => { $q_seq_id => 1 } };
	if ($aln->getIsReversed()) {
	    push @{ $self->{cacheMinus} }, $seed;
	} else {
	    push @{ $self->{cachePlus} }, $seed;
	}
    }
}

sub _doOverlapMerge {
    my $self = shift;

    print "(+) overlap merge: before count " . scalar(@{ $self->{cachePlus} })  . "\n";
    $self->_doMerge('+', 'om', 0);
    print "(+) overlap merge: after count " . scalar(@{ $self->{cachePlus} }). "\n";

    print "(-) overlap merge: before count " . scalar(@{ $self->{cacheMinus} })  . "\n";
    $self->_doMerge('-', 'om', 0);
    print "(-) overlap merge: after count " . scalar(@{ $self->{cacheMinus} }). "\n";
}

sub _doProximityMerge {
    my $self = shift;
    my $pm = $self->{pm};

    print "(+) proximity ($pm bp) merge: before count " . scalar(@{ $self->{cachePlus} }) . "\n";
    $self->_doMerge('+', 'pm', $pm);
    print "(+) proximity ($pm bp) merge: after count " . scalar(@{ $self->{cachePlus} }) . "\n";

    print "(-) proximity ($pm bp) merge: before count " . scalar(@{ $self->{cacheMinus} }) . "\n";
    $self->_doMerge('-', 'pm', $pm);
    print "(-) proximity ($pm bp) merge: after count " . scalar(@{ $self->{cacheMinus} }) . "\n";
}

sub _doClonelinkMerge {
    my $self = shift;

    my $cm = $self->{cm};

    print "(+) clonelink merge ($cm bp): before count ", scalar(@{ $self->{cachePlus} }), "\n";
    $self->_doMerge('+', 'cm', $cm);
    print "(+) clonelink merge ($cm bp): after count ", scalar(@{ $self->{cachePlus} }), "\n";

    print "(-) clonelink merge ($cm bp): before count ", scalar(@{ $self->{cacheMinus} }), "\n";
    $self->_doMerge('-', 'cm', $cm);
    print "(-) clonelink merge ($cm bp): after count ", scalar(@{ $self->{cacheMinus} }), "\n";
}

# actual merge
# NOTE: in place merge in the cache
#
sub _doMerge {
    my ($self, $strand, $merge_mode, $merge_param) = @_;

    my $vl = $self->getVerboseLevel();

    my $old_seeds = ($strand eq '+' ? $self->{cachePlus} : $self->{cacheMinus});
    my $old_tot =  scalar(@$old_seeds);

    for (my $i=0; $i<$old_tot-1; $i++) {
	my $os1 = $old_seeds->[$i]; next unless $os1; my $id1 = $os1->{id};

	for (my $j=$i+1; $j<$old_tot; $j++) {
	    $os1 = $old_seeds->[$i]; next unless $os1; $id1 = $os1->{id};
	    my $os2 = $old_seeds->[$j]; next unless $os2; my $id2 = $os2->{id};

	    print "_doMerge::[$i, $j], seeds [$id1, $id2]\n" if $vl >= 3;

	    my ($merge, $stop) = $self->_detectMerge($os1, $os2, $merge_mode, $merge_param);
	    if($merge) {
		$os1->{coords} =
		    DoTS::Gene::CompositeGenomeFeature::mergeCoordinates($os1->{coords},
									 $os2->{coords});
		push @{ $os1->{alns} }, @{ $os2->{alns} };
		foreach (keys %{ $os2->{qseqs} }) { $os1->{qseqs}->{$_} = 1; }
		$os2->[$j] = undef;
		print "mered seed at index $j to $i\n" if $vl >= 2;
	    }
	    last if $stop;
	}
    }

    my $new_seeds = [];
    foreach (@$old_seeds) { push @$new_seeds, $_ if $_; }
    if ($strand eq '+') { $self->{cachePlus} = $new_seeds; } else { $self->{cacheMinus} = $new_seeds; }
}

sub _detectMerge {
    my ($self, $seed1, $seed2, $merge_mode, $merge_param) = @_;

    my $coords1 = $seed1->{coords};
    my $gf1 = DoTS::Gene::GenomeFeature->new({ coords => $coords1 });
    my ($c1s, $c1e) = $gf1->getGenomicBoundaries();
    my $coords2 = $seed2->{coords};
    my $gf2 = DoTS::Gene::GenomeFeature->new({ coords => $coords2 });
    my ($c2s, $c2e) = $gf2->getGenomicBoundaries();

    my $shouldMerge = 0;
    my $stopHere = ($c2s - $c1e > $merge_param);

    if ($merge_mode eq 'om' && !$stopHere) {
	my $overlap = DoTS::Gene::GenomeFeature::getSpanOverlap($coords1, $coords2);
	$shouldMerge = ($overlap >= $merge_param);
    } elsif ($merge_mode eq 'cm') {
	my @clonelinks = $self->_getCloneLinks($seed1, $seed2, $merge_param);
	$shouldMerge = (scalar(@clonelinks) >= 1 && ($c2e-$c1s) <= $merge_param);
    } else {
	# do not merge intertwined seeds (e.g. seed2start < seed1end - "negative" distance)
	$shouldMerge = ($c2s-$c1e>=0 && $c2s-$c1e<=$merge_param);
    }

    return ($shouldMerge, $stopHere);
}

sub _getCloneLinks {
    my ($self, $seed1, $seed2)  = @_;

    my $dbh = $self->{dbh};
    my $epc = $self->{epc};

    my $qseqs1 = $seed1->{qseqs};
    my $qseqs2 = $seed2->{qseqs};
    # TRICKY: if two seeds have the same query seqs aligning to two location, 
    # exclude them in search for common clone
    my %shared_qseqs = DoTS::Gene::Util::sharedHashKeys($qseqs1, $qseqs2);

    my $bigSeedQ;
    my @smallSeedQ;
    if (scalar(keys %$qseqs1) > scalar(keys %$qseqs2)) {
	$bigSeedQ = $qseqs1; @smallSeedQ = keys %$qseqs2; 
    } else {
	$bigSeedQ = $qseqs2; @smallSeedQ = keys %$qseqs1; 
    }

    my $inClause;
    # HACK: oracle in clause does not allow >255 values
    if (scalar(@smallSeedQ) > 254) {
	print "WARNING: >254 qseqs in smaller of seed pair, use 254 for clonelink search";
	$inClause = join(',', $smallSeedQ[0..253]);
    } else {
	$inClause = join(',', @smallSeedQ);
    }
    my $sql = "select seq_id_1, lib_clone from $epc"
	. " where seq_id_1 != seq_id_2 and seq_id_2 in ($inClause)"
	. " union"
	. " select seq_id_2, lib_clone from $epc"
	. " where seq_id_1 != seq_id_2 and seq_id_1 in ($inClause)";
    my $sth = $dbh->prepareAndExecute($sql);

    my @cls;
    while (my ($seqId, $lnk) = $sth->fetchrow_array) {
	push @cls, $lnk if !$shared_qseqs{$seqId} && $bigSeedQ->{$seqId};
    }
    @cls;
}

sub _pack {
    my $self = shift;

    print "(+) packing seeds into composite genome features\n";
    my $cgfPlus = $self->_packOneStrand('+');

    print "(-) packing seeds into composite genome features\n";
    my $cgfMinus = $self->_packOneStrand('-');

    print "combine (+) and (-) results in order\n";
    return _combineOrderedGenomeFeatures($cgfPlus, $cgfMinus);
}

sub _packOneStrand {
    my ($self, $strand) = @_;

    my $gene_seeds = ($strand eq '+' ? $self->{cachePlus} : $self->{cachMinus});

    my @res = ();
    foreach my $seed (@$gene_seeds) {
	my $id = $seed->{id};
	my $coords = $seed->{coords};
	my $alns = $seed->{alns};
	if (scalar(@$alns) > 0) {
	    my $genomeId = $alns->[0]->getTargetExternalDbReleaseId();
	    my $chr = $alns->[0]->getChromosome();
	    my $gf_args = { id=>$id, genome_id=>$genomeId, chr=>$chr, coords=>$coords };
	    my @parts;
	    foreach (@$alns) {
		my $pid = $_->getBlatAlignmentId();
		my $pSeqId = $_->getQueryNaSequenceId();
		my @coords = $_->getTargetCoordinates();
		push @parts, { id=>$pid, na_sequence_id=>$pSeqId, coords=>\@coords };
	    }

	    push @res, DoTS::Gene::CompositeGenomeFeature->new($gf_args, \@parts);
	}
    }

    if ($strand eq '+') { $self->{cachePlus} = undef; } else { $self->{cacheMinus} = undef; }

    return \@res;
}

################ file scoped subroutines ################################

# combine ordered genome features from both strand
sub _combineOrderedGenomeFeatures {
    my ($pgfs, $mgfs) = @_;

    my @res;
    my $c = scalar(@$mgfs);
    my $j = 0;
    foreach my $g1 (@$pgfs) {
	my $g1s = $g1->getChromStart;
	my $g1e = $g1->getChromEnd;
	for (; $j < $c; $j++) {
	    my $g2 = $mgfs->[$j];
	    my $g2s = $g2->getChromStart;
	    my $g2e = $g2->getChromEnd;
	    if ($g2s >= $g1s) {
		push @res, $g1;
		last;
	    }
	    push @res, $g2;
	}
	push @res, $g1 if $j == $c;
    }
    return \@res;
}

1;
