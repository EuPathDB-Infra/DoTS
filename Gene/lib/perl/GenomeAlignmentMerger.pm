package DoTS::Gene::GenomeAlignmentMerger;

use strict;
use DoTS::Gene::Util;
use DoTS::Gene::GenomeFeature;
use DoTS::Gene::CompositeGenomeFeature;
use CBIL::Util::Disp;
use Carp;

# constructor
#
sub new {
    my($class, $db, $sorted_alignments, $merge_criteria) = @_;

    my $om = $merge_criteria->{overlap_merge};
    my $cm = $merge_criteria->{clonelink_merge};
    my $epc = $merge_criteria->{est_pair_cache};
    my $pm = $merge_criteria->{proximity_merge};
    my $self = { sa => $sorted_alignments, db => $db,
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

    if (defined $self->{om}) {
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
    print "Total number of CompositeGenomeFeatures in the region: " . scalar(@$res) . "\n";

    $self->{merged} = $res;
    return $self->{merged};
}

#####################################################################################

sub _seed {
    my $self = shift;

    my $srt_aln = $self->{sa};
    my $c = scalar(@$srt_aln);
    if ($c > 3000) {
        print "increasing object cache size to 3 * $c + 1000 (>10000)\n";
        $self->{db}->setMaximumNumberOfObjects(3 * $c + 1000);
    }

    foreach my $aln (@$srt_aln) {
	my $id = $aln->getBlatAlignmentId();
	$aln->retrieveFromDB();
	my $q_seq_id = $aln->getQueryNaSequenceId();
	if (!$q_seq_id) { warn("query seq id not found for blat alignment $id, skip"); next; }
	my @coords = $aln->getTargetCoordinates();
	my $isRev = $aln->getIsReversed();

	my $seed = { id => $id, coords => \@coords, alns => [$aln], qseqs => { $q_seq_id => 1 },
		 omc => 0, cmc => 0, pmc => 0 };
	if ($isRev) {
	    push @{ $self->{cacheMinus} }, $seed;
	} else {
	    push @{ $self->{cachePlus} }, $seed;
	}
    }
}

sub _doOverlapMerge {
    my $self = shift;
    my $om = $self->{om};

    print "(+) overlap (of at least $om bp) merge: before count " . scalar(@{ $self->{cachePlus} })  . "\n";
    $self->_doMerge('+', 'om', $om);
    print "(+) overlap merge: after count " . scalar(@{ $self->{cachePlus} }). "\n";

    print "(-) overlap (of at least $om bp) merge: before count " . scalar(@{ $self->{cacheMinus} })  . "\n";
    $self->_doMerge('-', 'om', $om);
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
	my $localBoundR;
	my @localRecurseQ = ();
	for (my $j=$i+1; $j<$old_tot; $j++) {
	    # TRICKY: seed of outer loop needs to be refreshed since it can change!
	    my $os1 = $old_seeds->[$i]; last unless $os1;
	    my $id1 = $os1->{id}; my $crd1 = $os1->{coords};
	    $localBoundR = (reverse @$crd1)[0]->[1];
	    my $os2 = $old_seeds->[$j]; next unless $os2;
	    my $id2 = $os2->{id}; my $crd2 = $os2->{coords};

	    print "coord1: " . join(':', map { $_->[0] . '-' . $_->[1] } @$crd1) . "\n" if $vl >= 3;
	    print "coord2: " . join(':', map { $_->[0] . '-' . $_->[1] } @$crd2) . "\n" if $vl >= 3;
	    print "_doMerge::[$i, $j], seeds [$id1, $id2]\n" if $vl >= 3;

	    my ($merge, $bounds) = $self->_detectMerge($os1, $os2, $merge_mode, $merge_param);
	    last if $bounds->[2] - $localBoundR > $merge_param;
	    $localBoundR = ($bounds->[1] > $bounds->[3] ? $bounds->[1] : $bounds->[3]);
	    if($merge) {
		my $mc = DoTS::Gene::CompositeGenomeFeature::mergeCoordinateSets($crd1, $crd2);
		$os1->{coords} = $mc;

		print "merged: " . join(':', map { $_->[0] . '-' . $_->[1] } @$mc) . "\n" if $vl >= 3;

		push @{ $os1->{alns} }, @{ $os2->{alns} };
		foreach (keys %{ $os2->{qseqs} }) { $os1->{qseqs}->{$_} = 1; }
		$os1->{"${merge_mode}c"}++; # increment merge count
		$old_seeds->[$j] = undef;
		print "mered seed at index $j to $i\n" if $vl >= 2;
	    } else {
		push @localRecurseQ, $j;
		print "no merge detected\n" if $vl >= 3; 
	    }
	}

	# localRecurse
	my $changed = scalar(@localRecurseQ);
	while ($changed) {
	    $changed = 0;
	    my $os1 = $old_seeds->[$i];
	    my $id1 = $os1->{id}; my $crd1 = $os1->{coords};

	    foreach my $j (@localRecurseQ) {
		my $os2 = $old_seeds->[$j]; next unless $os2;
		my $id2 = $os2->{id}; my $crd2 = $os2->{coords};

		my ($merge) = $self->_detectMerge($os1, $os2, $merge_mode, $merge_param);
		if($merge) {
		    $changed = 1;
		    my $mc = DoTS::Gene::CompositeGenomeFeature::mergeCoordinateSets($crd1, $crd2);
		    $os1->{coords} = $mc;
		    print "merged: detected in local recurse over " . scalar(@localRecurseQ)
			. " indices not merged in first pass\n";

		    push @{ $os1->{alns} }, @{ $os2->{alns} };
		    foreach (keys %{ $os2->{qseqs} }) { $os1->{qseqs}->{$_} = 1; }
		    $os1->{"${merge_mode}c"}++; # increment merge count
		    $old_seeds->[$j] = undef;
		    print "mered seed at index $j to $i\n" if $vl >= 2;
		}
	    }
	}
    }

    my $new_seeds = [];
    foreach (@$old_seeds) { push @$new_seeds, $_ if $_; }
    if ($strand eq '+') { $self->{cachePlus} = $new_seeds; } else { $self->{cacheMinus} = $new_seeds; }
}

sub _detectMerge {
    my ($self, $seed1, $seed2, $merge_mode, $merge_param) = @_;

    &confess("merge param is required") unless defined $merge_param;

    my $coords1 = $seed1->{coords};
    my $gf1 = DoTS::Gene::GenomeFeature->new({ coords => $coords1 });
    my ($c1s, $c1e) = $gf1->getGenomicBoundaries();
    my $coords2 = $seed2->{coords};
    my $gf2 = DoTS::Gene::GenomeFeature->new({ coords => $coords2 });
    my ($c2s, $c2e) = $gf2->getGenomicBoundaries();

    my $shouldMerge = 0;
    if ($merge_mode eq 'om') {
	my $overlap = DoTS::Gene::GenomeFeature::getSpanOverlap($coords1, $coords2);
	$shouldMerge = ($overlap >= $merge_param);
    } else {
	if ($merge_mode eq 'cm') {
	    my @clonelinks = $self->_getCloneLinks($seed1, $seed2, $merge_param);
	    $shouldMerge = (scalar(@clonelinks) >= 1 && ($c2e-$c1s) <= $merge_param);
	} else {
	    # do not merge intertwined seeds (e.g. seed2start < seed1end - "negative" distance)
	    $shouldMerge = ($c2s-$c1e>=0 && $c2s-$c1e<=$merge_param);
	}
	if ($shouldMerge) {
	    my $hasFlc1 = $self->_hasFullLengthClone($seed1);
	    my $hasFlc2 = $self->_hasFullLengthClone($seed2);
	    if ($hasFlc1 && $hasFlc2) {
		$shouldMerge = 0;
		print "merge disqualified: both seeds contain full length clone(s)\n";
	    }
	}
    }

    return ($shouldMerge, [$c1s, $c1e, $c2s, $c2e]);
}

sub _hasFullLengthClone {
    my ($self, $seed) = @_;

    my $dbh = $self->{db}->getQueryHandle();
    my @qseqs = keys %{ $seed->{qseqs} };
    my $c = scalar(@qseqs);

    my @bin = ();
    for (my $i=0; $i<$c; $i++) {
	push @bin, $qseqs[$i];
	unless ($i % 253) {
	    my $sql = "select count(*) from DoTS.AssemblySequence s, DoTS.ExternalNaSequence e "
		. "where s.na_sequence_id = e.na_sequence_id "
		. "and s.assembly_na_sequence_id in (" . join(',', @bin) . ") "
		. "and e.external_database_release_id = 992";
	    my $sth = $dbh->prepareAndExecute($sql);
	    my $hasRefseq = $sth->fetchrow_array();
	    return 1 if $hasRefseq;

	    $sql = "select count(*) from DoTS.AssemblySequence s, "
		. "dots.dbrefnasequence dbrn, sres.dbref dbr "
		. "where s.na_sequence_id = dbrn.na_sequence_id "
		. "and dbrn.db_ref_id = dbr.db_ref_id "
		. "and s.assembly_na_sequence_id in (" . join(',', @bin) . ") "
		. "and dbr.external_database_release_id in (8494, 8495, 7414)";
	    $sth = $dbh->prepareAndExecute($sql);
	    my $hasMgcOrFantom = $sth->fetchrow_array();
	    return 1 if $hasMgcOrFantom;

	    @bin = ();
	}
    }
    return 0;
}

sub _getCloneLinks {
    my ($self, $seed1, $seed2)  = @_;

    my $dbh = $self->{db}->getQueryHandle();
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
    my $sql = "select dt_id_2, lib_clone from $epc"
	. " where dt_id_1 != dt_id_2 and dt_id_1 in ($inClause)";
    my $sth = $dbh->prepareAndExecute($sql);

    my @cls_all;
    while (my ($seqId, $lnk) = $sth->fetchrow_array) { push @cls_all, [$seqId, $lnk]; }
    my @cls;
    foreach (@cls_all) {
	my ($seqId, $lnk) = @$_;
	push @cls, $lnk if !$shared_qseqs{$seqId} && $bigSeedQ->{$seqId};
    }

    my $vl = $self->getVerboseLevel();
    if ($vl && scalar(@cls)) {
	print "*** qseq in small seed: " . join(', ', @smallSeedQ) . "\n" if $vl >= 2;
	print "*** qseq in big seed: " . join(', ', keys %$bigSeedQ) . "\n" if $vl >= 2;
	print "*** qseq in both seeds: " . join(', ', keys %shared_qseqs) . "\n" if $vl >= 2;
        print "*** links all: " . join(', ', map { $_->[1] } @cls_all) . "\n" if $vl >= 2;
	print "*** links informative: " . join(', ', @cls) . "\n";
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

    my $gene_seeds = ($strand eq '+' ? $self->{cachePlus} : $self->{cacheMinus});

    my @res = ();
    foreach my $seed (@$gene_seeds) {
	my $id = $seed->{id};
	my $coords = $seed->{coords};
	my $alns = $seed->{alns};
	if (scalar(@$alns) > 0) {
	    my $genomeId = $alns->[0]->getTargetExternalDbReleaseId();
	    my $chr = $alns->[0]->getChromosome();
	    my $str = $alns->[0]->getStrand();
	    my $gf_args = { id=>$id, genome_id=>$genomeId, chr=>$chr, strand=>$str, coords=>$coords };
	    my @parts;
	    foreach (@$alns) {
		my $pid = $_->getBlatAlignmentId();
		my $pSeqId = $_->getQueryNaSequenceId();
		my @coords = $_->getTargetCoordinates();
		push @parts, { id=>$pid, na_sequence_id=>$pSeqId, coords=>\@coords };
	    }

	    my $cgf = DoTS::Gene::CompositeGenomeFeature->new($gf_args, \@parts);
	    $cgf->setAnnotationProperty('omc', $seed->{omc});
	    $cgf->setAnnotationProperty('cmc', $seed->{cmc});
	    $cgf->setAnnotationProperty('pmc', $seed->{pmc});
	    push @res, $cgf;
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

    for (; $j < $c; $j++) { push @res, $mgfs->[$j]; }

    return \@res;
}

1;
