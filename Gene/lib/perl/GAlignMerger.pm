#!/usr/bin/perl

# -------------------------------------------------------
# create genes from genomic alignments of DoTS assemblies
#
# Yongchang Gan July 26, 2002
#
# -------------------------------------------------------

package DoTS::Gene::GAMerger;

use strict;

use DBI;
use Cwd;
use Util;
use GeneModel;

# constants for quality filter
#
my $MRNA_VERYGOOD = 1;
my $MRNA_VERYGOOD_GAP = 2;
my $MRNA_GOOD = 3;
my $MRNA_OK = 4;
my $SPLICED_VERYGOOD = 5;
my $SPLICED_VERYGOOD_GAP = 6;
my $SPLICED_GOOD = 7;
my $SPLICED_OK = 8;
my $VERYGOOD = 9;
my $VERYGOOD_GAP = 10;
my $GOOD = 11;
my $OK = 12;

my $MIN_BLAT_SCORE = 85;

# constants for merge mode selection
#
my $MERGE_BY_ALIGNMENT = 1;
my $MERGE_BY_CLONE = 2;
my $MERGE_BY_LOCATION = 3;

# constant for intron size cutoff
#
my $INTRON_SIZE_CUTOFF = 15;

$| = 1;

# constructor
#
sub new {
    my($class, $args) = @_;

    my $self = {
        taxon_id => $args->{taxon_id},
	chrom_id => $args->{chrom_id},
        start => $args->{start},
        end => $args->{end},
        quality_filter => $args->{quality_filter},
	exclude_singleton => $args->{exclude_singleton},
    };

    die "must specify taxon_id and chrom_id (DoTS.VirtualSequence.na_sequence_id)\n"
      unless $self->{taxon_id} && $self->{chrom_id};

    $self->{alignment_merge} = $args->{alignment_merge};
    $self->{clone_merge} = $args->{clone_merge};
    $self->{location_merge} = $args->{location_merge};
    $self->{include_desc} = $args->{show_desc};
    $self->{test} = $args->{test};
    $self->{debug} = $args->{debug};
    $self->{geneClusters} = { sortedKeys => [], clusters => undef, };
    $self->{genes} = [];
    $self->{plusSeeds} = [];
    $self->{minusSeeds} = [];

    my $ss;
    if ($self->{taxon_id} == 14) {
	$ss = 'm';
    } elsif ($self->{taxon_id} == 8) {
	$ss = 'h';
    } else {
	die "species $self->{taxon_id} not supported!\n";
    }

    my $fi = $args->{common_clone_image_id};
    $fi = &cwd . "\/common_clone_${ss}.image_id" unless $fi;
    my $fw = $args->{common_clone_washu_name};
    $fw = &cwd . "\/common_clone_${ss}.washu_name" unless $fw;
    $self->{commonCloneCacheFile} = { image_id => $fi, washu_name => $fw };
    $self->{commonCloneCache} = { image_id => undef, washu_name => undef };
    my $base =  "${ss}.$self->{chrom_id}.$self->{start}.$self->{end}"
             . ".$self->{quality_filter}" . ($self->{exclude_singleton} ? '.xs' : '');
    $self->{statusFile} = { dir => &cwd, suffix => ".sql_running", base => $base };

    # quality threshold for gene content
    my $qf = $self->{quality_filter};
    # mRNA-containment
    if ($qf==$MRNA_VERYGOOD || $qf==$MRNA_VERYGOOD_GAP || $qf==$MRNA_GOOD || $qf==$MRNA_OK) {
        $self->{mrna_clause} = " and a.contains_mrna = 1 ";
    }
    # splice
    if ($qf==$SPLICED_VERYGOOD || $qf==$SPLICED_VERYGOOD_GAP || $qf==$SPLICED_GOOD || $qf==$SPLICED_OK) {
        $self->{splice_clause} = " and b.max_target_gap >= $INTRON_SIZE_CUTOFF ";
    }
    # blatAlignmentQuality
    if ($qf == $MRNA_VERYGOOD || $qf == $SPLICED_VERYGOOD || $qf == $VERYGOOD) {
	$self->{quality_clause} = " and b.blat_alignment_quality_id = 1 ";
    } elsif ($qf == $MRNA_VERYGOOD_GAP || $qf == $SPLICED_VERYGOOD_GAP || $qf == $VERYGOOD_GAP) {
	$self->{quality_clause} = " and b.blat_alignment_quality_id <= 2 ";
    } elsif ($qf == $MRNA_GOOD || $qf == $SPLICED_GOOD || $qf == $GOOD) {
	$self->{quality_clause} = " and (b.blat_alignment_quality_id <= 3 "
	  . "or (b.is_best_alignment = 1 and score >= $MIN_BLAT_SCORE))";
    } elsif ($qf == $MRNA_OK || $qf == $SPLICED_OK || $qf == $OK) {
	$self->{quality_clause} = " and b.blat_alignment_quality_id <= 4 ";
    }

    # gene id base
    if ($self->{taxon_id} == 8) {
      $self->{gene_id_base_plus} = 1000000;
      $self->{gene_id_base_minus} = 2000000;
    } elsif ($self->{taxon_id} == 14) {
      $self->{gene_id_base_plus} = 3000000;
      $self->{gene_id_base_minus} = 4000000;
    } else {
      die "taxon $self->{taxon_id} not supported\n";
    }
    bless $self, $class;
    return $self;
}

########## subroutine ###################################################

# reture array of genes in this genomic region,
# create genes from genomic alignments if necessary
#
sub getGenes {
    my $self = shift;
    return $self->{genes} if scalar(@{ $self->{genes} }) > 0;

    my $SUB = 'getGenes';
    print "<PRE>\n" if $self->{debug};
    print "#debug ${SUB}::initially each gene cluster as a gene\n" if $self->{debug};
    $self->_initializeGeneSeeds;

    if ($self->{alignment_merge}) {
	print "#debug ${SUB}::merge clusters with alignment overlap\n" if $self->{debug};
	$self->_doAlignmentMerge;
    }

    if ($self->{clone_merge} > 0) {
	print "#debug ${SUB}::merge clusters within specified distance if supported by clone info\n"
	    if $self->{debug};
	$self->_doCloneMerge;
    }

    if ($self->{location_merge} > 0) {
	print "#debug ${SUB}::merge clusters within specified distance on genome\n" if $self->{debug};
	$self->_doLocationMerge;
    }

    print "#debug ${SUB}::flattern genomic alignments for each gene to create gene models\n"
	if $self->{debug};
    $self->_createGeneModels;
    print "#debug ${SUB}::total number of genes: ", scalar(@{ $self->{genes} }), "\n" if $self->{debug};
    print "</PRE>\n" if $self->{debug};

    return $self->{genes};
}

# return a hash of gene clusters
# and the genomic alignments by their component DoTS assemblies
# NOTE: a gene cluster is split to geneId- and geneId+
#       if all assemblies do not align to the same strand of genome
#       also split gene clusters where alignments of component assemblies do not overlap
#
sub getGeneClusterGenomicAlignments{
    my $self = shift;
    return $self->{geneClusters} if ($self->{geneClusters}->{clusters});

    my $SUB = 'getGeneClusterGenomicAlignments';

    my $chrom_id = $self->{chrom_id};
    my $taxon_id = $self->{taxon_id};
    my $start = $self->{start};
    my $end = $self->{end};
    my $desc = $self->{include_desc};
    my $exc_sgl = $self->{exclude_singleton};

    my $dbh = Util::getLogin();
    my ($sql, $sth, %res);
    my $sortedKeys = [];

    $sql = "SELECT r.gene_id, a.na_sequence_id,  a.contains_mrna, a.length, "
        . ($desc ? 'a.description,' : '') . " b.blat_alignment_id, b.number_of_spans, b.tstarts, "
        . "b.blocksizes, b.target_start, b.target_end, b.is_reversed, b.is_consistent "
        . "FROM DoTS.RNA r, DoTS.RNAInstance ri, DoTS.NAFeature naf, "
	. "  DoTS.Assembly a, DoTS.BlatAlignment b "
        . "WHERE r.rna_id = ri.rna_id and ri.na_feature_id = naf.na_feature_id "
	. "and naf.na_sequence_id = a.na_sequence_id "
	. "and b.query_na_sequence_id = a.na_sequence_id "
        . " " . $self->{mrna_clause} . " "
        . " " . $self->{splice_clause} . " "
        . " " . $self->{quality_clause} . " "
        . ($exc_sgl ? "and a.number_of_contained_sequences > 1 " : "")
        . "and b.target_na_sequence_id = $chrom_id "
        . ($end ? "and b.target_start <= $end " : "")
        . ($start ? "and b.target_end >= $start " : "")
	. "order by b.target_start asc, b.target_end asc ";

    my $statFile = $self->{statusFile};
    my ($dir, $base, $suffix) = ($statFile->{dir}, $statFile->{base}, $statFile->{suffix});

    my @running = ();
    do {
        my $str = `ls $dir\/\*$suffix` or print "";
	@running = split(/\n/, $str);
	if (scalar(@running) > 0) {
	    print "#info ${SUB}::", scalar(@running), " other chrom(s) currently running: ",
	          join(', ', @running), "\n";
	    print "#info ${SUB}::sleeping 5 minute...\n";
	    sleep 300;
	}
    } while scalar(@running > 0);

    my $sf = "$dir\/$base$suffix";
    open STATUS_FILE, ">$sf";
    print "#info\t${SUB}::writing to status file $sf\n";
    print STATUS_FILE $sql;
    close STATUS_FILE;

    $sth = $dbh->prepare($sql);
    print "#debug ${SUB}::running sql: $sql \n" if $self->{debug};
    $sth->execute();
    while (my $h = $sth->fetchrow_hashref('NAME_lc')) {
	my $gid = $h->{'gene_id'};
	my $bid = $h->{'blat_alignment_id'};
        my $rev = $h->{'is_reversed'};

	$rev = !$rev if &reversedTranscript($dbh, $bid);
        $h->{'is_reversed'} = $rev;

        delete $h->{gene_id};
        my $suffix = $rev ? '-' : '+';
	my $gckey = "$gid$suffix";
	# Tricky: handle the situation where a gene cluster align to multiple genomic locations
	if (exists $res{$gckey}) {
	    my $lastKey = $gckey;
	    while (exists $res{"$lastKey$suffix"}) { $lastKey .= $suffix; }
	    my $lastAligns = $res{$lastKey};

	    my $dist = $h->{target_start} - $lastAligns->[scalar(@$lastAligns)-1]->{target_end};
	    if ($dist <= $self->{location_merge}) {
		push @{ $res{$lastKey} }, $h;
	    } else {
		push @$sortedKeys, "$lastKey$suffix";
		$res{"$lastKey$suffix"} = [$h];
	    }
	} else {
            push @$sortedKeys, $gckey;
	    $res{$gckey} = [$h];
	}
    }

    print "#info ${SUB}::deleting status file $sf\n";
    unlink $sf or die "could not unlink status file $sf\n";

    $self->{geneClusters}->{clusters} = \%res;
    $self->{geneClusters}->{sortedKeys} = $sortedKeys;

    return $self->{geneClusters};
}

sub getFixedTarget {
    my $self = shift;
    return ($self->{taxon_id}, $self->{chrom_id}, $self->{start}, $self->{end});
}

sub getStatusFile {
    my $self = shift;
    return $self->{statusFile};
}

############ private subroutines ##############################################

# merge clusters with genomic alignment overlap
#
sub _doAlignmentMerge {
    my $self = shift;

    my $SUB = '_doAlignmentMerge';
    print "#debug ${SUB}::performing merge by genomic alignment overlap...\n"
	if $self->{debug}; 
    print "#debug ${SUB}::plus strand gene counts before merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(1, $MERGE_BY_ALIGNMENT, 0);
    print "#debug ${SUB}::plus strand gene counts after merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};

    print "#debug ${SUB}::minus strand gene counts before merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(0, $MERGE_BY_ALIGNMENT, 0);
    print "#debug ${SUB}::minus strand gene counts after merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
}

# merge clusters within specified distance on genome
#
sub _doLocationMerge {
    my $self = shift;
    my $location_merge = $self->{location_merge};

    my $SUB = '_doLocationMerge';
    print "#debug ${SUB}::performing merge by location ($location_merge bp)...\n"
	if $self->{debug}; 
    print "#debug ${SUB}::plus strand gene counts before merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(1, $MERGE_BY_LOCATION, $location_merge);
    print "#debug ${SUB}::plus strand gene counts after merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};

    print "#debug ${SUB}::minus strand gene counts before merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(0, $MERGE_BY_LOCATION, $location_merge);
    print "#debug ${SUB}::minus strand gene counts after merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
}

# merge clusters within specified distance if supported by clone info
#
sub _doCloneMerge {
    my $self = shift;
    my $clone_merge = $self->{clone_merge};

    my $SUB = '_doCloneMerge';
    print "#debug ${SUB}::performing merge by clone info ($clone_merge bp)...\n"
	if $self->{debug}; 
    print "#debug ${SUB}::plus strand gene counts before merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(1, $MERGE_BY_CLONE, $clone_merge);
    print "#debug ${SUB}::plus gene counts after merge: ", scalar(@{ $self->{plusSeeds} }), "\n"
	if $self->{debug};

    print "#debug ${SUB}::minus strand gene counts before merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
    $self->_doGeneSeedMerge(0, $MERGE_BY_CLONE, $clone_merge);
    print "#debug ${SUB}::minus gene counts after merge: ", scalar(@{ $self->{minusSeeds} }), "\n"
	if $self->{debug};
}

# actual gene cluster merge algorithm, by clone info and/or genomic location
# NOTE: in place merge
#
sub _doGeneSeedMerge {
    my ($self, $strand, $merge_mode, $merge_param) = @_;

    my $old_seeds = $strand ? $self->{plusSeeds} : $self->{minusSeeds};
    my $old_tot =  scalar(@$old_seeds);

    my $SUB = '_doGeneSeedMerge';
    for (my $i=0; $i<$old_tot-1; $i++) {
	my $c_keys_1 = $old_seeds->[$i]->{cluster_keys};
	next unless $c_keys_1;

	#debug
	if ($self->{debug} >= 3) {
	    print "#debug \t\t${SUB}::i = $i, content of gene_seeds:\n";
	    print "#debug \t\t${SUB}::";
	    foreach (@$old_seeds) { 
		print "cluster_keys => ",
		defined($_->{cluster_keys}) ? join(', ', @{ $_->{cluster_keys} }) : 'undef', " :: ";
	    }
	    print "\n";
	}

	for (my $j=$i+1; $j<$old_tot; $j++) {
            $c_keys_1 = $old_seeds->[$i]->{cluster_keys};
	    next unless $c_keys_1;
            my $c_keys_2 = $old_seeds->[$j]->{cluster_keys};
	    # NOTE this can occur with MERGE_BY_ALIGNMENT in the situation of gene within gene!!!
	    next unless $c_keys_2;
	    print "#debug \t\t${SUB}::j = $j, cluster keys at $i: ", join(',', @$c_keys_1),
	          ", cluster keys at $j: ", join(',', @$c_keys_2), "\n" if $self->{debug} >= 3;

	    my ($merge, $stop) = $self->_detectGeneSeedMerge($c_keys_2, $c_keys_1, $merge_mode, $merge_param);
	    if($merge) {
		push @{ $old_seeds->[$i]->{cluster_keys} }, @{ $old_seeds->[$j]->{cluster_keys} };
		$old_seeds->[$j]->{cluster_keys} = undef;
		print "#debug \t\t${SUB}::appended gene seed at " . ($j) . " to $i\n" if $self->{debug} >= 3;
	    }
	    last if $stop;
	}
    }

    my $new_seeds = [];
    foreach (@$old_seeds) { push @$new_seeds, $_ if $_->{cluster_keys}; }
    if ($strand) { $self->{plusSeeds} = $new_seeds; } else { $self->{minusSeeds} = $new_seeds; }
}

# move cluster ids from one gene seed to another 
# from_key and to_keys are already one the same strand so no need to check that
#
sub _detectGeneSeedMerge {
    my ($self, $from_keys, $to_keys, $merge_mode, $merge_param) = @_;

    my $shouldMerge = 0;
    my ($c1s, $c1e) = $self->_getGeneSeedOuterBounds($to_keys);
    my ($c2s, $c2e) = $self->_getGeneSeedOuterBounds($from_keys);
    my $stopHere = ($c2s - $c1e > $merge_param);

    my $SUB = '_detectGeneSeedMerge';
    my ($alignmentOverlap, $sameClone);
    if ($merge_mode == $MERGE_BY_ALIGNMENT && !$stopHere) {
	my ($pct_from, $pct_to) = $self->_getGeneSeedAlignmentOverlap($from_keys, $to_keys);
	$alignmentOverlap = ($pct_from > 0 && $pct_to > 0);
    } elsif ($merge_mode == $MERGE_BY_CLONE) {
	# TRICKY: if two gene seeds has the same cluster aligning to two location, 
	# exclude that cluster in search for common clone info
	my ($cks1, $cks2, $exc) = &excludeCommonClusterKeys($from_keys, $to_keys);
	if ($self->{debug} >= 2) {
	    print "#debug \t${SUB}::gene seeds after excluding common clusters(", join(',', @$exc) . "): ";
	    print "from ", join(',', @$cks1), " to ", join(',', @$cks2), "\n";
	}

	my ($tot_seqs1, $tot_seqs2, $cloneLinks) = (-1, -1, []);
	if (scalar(@$cks1) >= 1 && scalar(@$cks2) >= 1) {
	    ($tot_seqs1, $tot_seqs2, $cloneLinks) = $self->_getGeneSeedCloneLinks($cks1, $cks2);
	}

#        if ($self->{debug} >= 2) {
        if (scalar(@$cloneLinks)>0) {
	    print "#debug \t${SUB}::gene seed " . join(',', @$cks1) . " has $tot_seqs1 seqs\n"
		. "#debug \t${SUB}::gene seed " . join(',', @$cks2) . " has $tot_seqs2 seqs\n"
	        . "#debug \t${SUB}::common clone info: " . (scalar(@$cloneLinks) > 0 ? "" : 'NONE') . "\n";
	    foreach (@$cloneLinks) { print "#debug \t$SUB::", $_. "\n"; }
	}

	if ($tot_seqs1 && $tot_seqs2) {
	    $sameClone = scalar(@$cloneLinks)*1.0/$tot_seqs1 > 0 && scalar(@$cloneLinks)*1.0/$tot_seqs2 > 0;
	} else { $sameClone = 0; }
    }

    if ($self->{'debug'} >= 2) {
	my $msg = 'genomic alignment overlap';
	$msg = 'clone' if $merge_mode == $MERGE_BY_CLONE;
	$msg = 'location' if $merge_mode == $MERGE_BY_LOCATION;
	print "#debug \t${SUB}::Just examined merge of gene seed ", join(',', @$from_keys), " at [$c2s\-$c2e] to ",
        join(',', @$to_keys), " at [$c1s\-$c1e] (", ($c2s-$c1e), " bp away) by $msg:";
    }

    my $alignDetected = ($merge_mode == $MERGE_BY_ALIGNMENT && $alignmentOverlap);
    my $cloneDetected = ($merge_mode==$MERGE_BY_CLONE && $c2e-$c1s<=$merge_param && $sameClone);
    my $locatDetected = ($self->{alignment_merge} ?
			 ($merge_mode==$MERGE_BY_LOCATION && $c2s-$c1e>=0 && $c2s-$c1e<=$merge_param) :
			 ($merge_mode==$MERGE_BY_LOCATION && $c2s-$c1e<=$merge_param));

    if ($alignDetected || $cloneDetected || $locatDetected) {
	print " merge detected!" if $self->{debug} >= 2;
        $shouldMerge = 1;
    } elsif ($self->{debug} >=2) {
	print " no merge because no alignment overlap" if $merge_mode==$MERGE_BY_ALIGNMENT && !$sameClone;
	print " no merge because no clone info support" if $merge_mode==$MERGE_BY_CLONE && !$sameClone;
	print " no merge because distance ", ($c2s-$c1e), " > $merge_param" 
	    if ($merge_mode==$MERGE_BY_CLONE || $merge_mode==$MERGE_BY_LOCATION) && $c2s-$c1e > $merge_param;
    }
    print "\n\n" if $self->{debug} >= 2;

    return ($shouldMerge, $stopHere);
}

# given a pair of gene seeds (lists of cluster ids), 
# find the overlap of their genomic alignments
#
sub _getGeneSeedAlignmentOverlap {
    my ($self, $cks1, $cks2) = @_;

    my $SUB = '_getGeneSeedAlignmentOverlap';
    print "#debug \t${SUB}::checking genomic alignment overlap between ", join(',', @$cks1),
          " and ", join(',', @$cks2) , "\n" if $self->{debug} >= 2;
    my $clusterAlignHash = $self->getGeneClusterGenomicAlignments->{clusters};
    my $gm1 = &makeOneGeneModel(1, -1, $cks1, $clusterAlignHash, $self->{debug});
    my $gm2 = &makeOneGeneModel(1, -1, $cks2, $clusterAlignHash, $self->{debug});

    my $geneSize1 = $gm1->getGeneSize;
    my $geneSize2 = $gm2->getGeneSize;
    my $exonCoords1 = $gm1->getExonCoords;
    my $exonCoords2 = $gm2->getExonCoords;

    my $overlap = 0;
    my $c1 = scalar(@$exonCoords1);
    my $c2 = scalar(@$exonCoords2);

    if ($self->{debug} >= 2 or $geneSize1 == 0 or $geneSize2 == 0) {
      if ($self->{debug} < 2) {
	print "#WARNING \t${SUB}::checking genomic alignment overlap between "
	  . join(',', @$cks1) . " and " . join(',', @$cks2) . "\n";
      }
      print "#debug \t\t${SUB}::exonCoords1(c1:$c1): ", join(',', map { $_->[0] . '-' . $_->[1] } @$exonCoords1), "\n";
      print "#debug \t\t${SUB}::exonCoords2(c2:$c2): ", join(',', map { $_->[0] . '-' . $_->[1] } @$exonCoords2), "\n";
    }

    my $lastJ = 0;
    for(my $i=0; $i<$c1; $i++) {
	my ($s1, $e1) = @{ $exonCoords1->[$i] };
	for(my $j = $lastJ; $j<$c2; $j++) {
	    my ($s2, $e2) = @{ $exonCoords2->[$j] };

            print "#debug \t\t${SUB}::i=$i, j=$j, lastJ=$lastJ, (s1,e1)=($s1,$e1), (s2,e2)=($s2,$e2)\n"
                if $self->{debug} >= 3;

            if ($e1 < $s2) { 
                 $lastJ = $j; 
                 # TRICKY: if the next inner loop is out of range, the next outer
                 # should begin with this inner loop index, not incremented one
                 $lastJ-- if $lastJ;
                 last;
            }

	    next if $s1 > $e2;
	    my $o = 0;
	    if ($s2 >= $s1 && $s2 <= $e1) {
		if ($e1 <= $e2) { $o = $e1 - $s2; } else { $o = $e2 - $s2; }
	    } elsif ($s1 >= $s2 && $s1 <= $e2) {
		if ($e1 <= $e2) { $o = $e1 - $s1; } else { $o = $e2 - $s1; }
	    }

	    print "#debug \t\t${SUB}::($s1, $e1) and ($s2, $e2) have overlap of $o bp\n" 
		if $self->{debug} >= 3;
	    $overlap += $o;
	}
    }

    print "#WARNING \t${SUB}::(geneSize1, geneSize2, overlap) = ($geneSize1, $geneSize2, $overlap)\n"
	if $self->{debug}>= 2 or $geneSize1 == 0 or $geneSize2 == 0;


    my ($pct_o1, $pct_o2);
    $pct_o1 = $overlap*1.0/$geneSize1 if $geneSize1;
    $pct_o2 = $overlap*1.0/$geneSize2 if $geneSize2;

    return ($pct_o1, $pct_o2);
}

# given a gene seed (a list of cluster ids), find the outer bounds
#
sub _getGeneSeedOuterBounds {
    my $self = shift;
    my $c_keys = shift;

    my $clusters = $self->getGeneClusterGenomicAlignments->{clusters};
    my ($gs, $ge, $is_rev) = (0, 0, undef);
    foreach my $c (@$c_keys) {
        my $alignments = $clusters->{$c};
        foreach my $a (@$alignments) {
            my $s = $a->{target_start};
            my $e = $a->{target_end};
            $gs = $s if $gs == 0 || $s < $gs;
            $ge = $e if $e > $ge;

            my $r = $a->{is_reversed};
            die "Error: gene seed contain both \+ and \- strand content!\n"
                if (defined($is_rev) && $is_rev != $r); 
            $is_rev = $r;
        }
    }
    ($gs, $ge, ($is_rev ? 0 : 1));
}

# given a gene seed (a list of cluster ids), find shared clone info.
#
sub _getGeneSeedCloneLinks {
    my ($self, $c_keys_1, $c_keys_2)  = @_;

    my $SUB = '_getGeneSeedCloneLinks';

    my @cids1 = &getUniqueIdsFromClusterKeys($c_keys_1);
    my @cids2 = &getUniqueIdsFromClusterKeys($c_keys_2);
    my $cls = [];
    my ($cis, $cws);
    my ($num_seq1, $num_seq2) = (0, 0);
    my $ccc_image_id = $self->_getCommonCloneCache('image_id');
    my $ccc_washu_name = $self->_getCommonCloneCache('washu_name');
 
    foreach my $cid1 (@cids1) {
	foreach my $cid2 (@cids2) {
	    print "#debug \t${SUB}::checking common clone info of cluster $cid1 in ",
                "gene seed 1 and $cid2 in gene seed 2" if $self->{debug} >= 2;

            my ($iR, $wR);
	    $cis = $ccc_image_id->{"$cid1$cid2"};
	    unless ($cis) {
                $cis = $ccc_image_id->{"$cid2$cid1"};
                $iR = 1;
            }
	    $cws = $ccc_washu_name->{"$cid1$cid2"};
            unless ($cws) {
	        $cws = $ccc_washu_name->{"$cid2$cid1"};
                $wR = 1;
            }

            my ($ns1, $ns2) = (0, 0); 
            if ($cis) {
                $ns1 = $ccc_image_id->{num_seqs_hash}->{$cid1};
                $ns2 = $ccc_image_id->{num_seqs_hash}->{$cid2};
                push @$cls, @$cis;
                print "#debug ${SUB}::gene clusters ",
                      ($iR ? "$cid2 and $cid1" : "$cid1 and $cid2"),
                      " share image_ids as in ", join(',', @$cis), "\n"
                      if $self->{debug} == 1;
            }
            if ($cws) {
                $ns1 = $ccc_washu_name->{num_seqs_hash}->{$cid1};
                $ns2 = $ccc_washu_name->{num_seqs_hash}->{$cid2};
                push @$cls, @$cws;
                print "#debug ${SUB}::gene clusters ", 
                      ($wR ? "$cid2 and $cid1" : "$cid1 and $cid2"),
                      " share washu_names as in ", join(',', @$cws), "\n"
                      if $self->{debug} == 1;
            }

            $num_seq1 += $ns1;
            $num_seq2 += $ns2;

	    if ($self->{debug} >= 2) {
		if ($cis && $cws) {
		    print ": found by image_ids as in ", join(',', @$cis),
                          " and washu_names as in ", join(',', @$cws), "\n";
		} elsif ($cis) {
		    print ": found by image_ids as in ", join(',', @$cis), "\n";
		} elsif ($cws) {
		    print ": found by washu_names as in ", join(',', @$cws), "\n";
		} else {
		    print ": not found\n";
		}
	    }
	}
    }

    ($num_seq1, $num_seq2, $cls);
}

# take each gene cluster as a gene seed, one list for each strand
#
sub _initializeGeneSeeds {
    my $self = shift;
    my $sortedKeys = $self->getGeneClusterGenomicAlignments->{sortedKeys};
    my $plusSeeds = [];
    my $minusSeeds = [];
    foreach my $ck (@$sortedKeys) {
	if ($ck =~ /\-$/) {
	    push @$minusSeeds, { cluster_keys => [$ck], };
	} else {
	    push @$plusSeeds, { cluster_keys => [$ck], };
	}
    }
    $self->{plusSeeds} = $plusSeeds;
    $self->{minusSeeds} = $minusSeeds;
}

# flatten genomic alignments for each gene to create gene models
# Note: in place model population
#
sub _createGeneModels {
    my $self = shift;

    my $SUB = '_createGeneModels';
    print "#debug ${SUB}::create gene models on plus strand\n" if $self->{debug};
    $self->_createOneStrandGeneModels(1);

    print "#debug ${SUB}::create gene models on minus strand\n" if $self->{debug};
    $self->_createOneStrandGeneModels(0);

    print "#debug ${SUB}::combine ordered genes from both strand\n" if $self->{debug};
    $self->_combineOrderedGeneLists;
}

# flatten genomic alignments for each gene to create gene models
# Note: in place model population
#
sub _createOneStrandGeneModels {
    my ($self, $strand) = @_;

    my $gene_seeds = ($strand ? $self->{plusSeeds} : $self->{minusSeeds});
    my $gid = ($strand ? $self->{gene_id_base_plus} : $self->{gene_id_base_minus});
    my $gcga = $self->getGeneClusterGenomicAlignments;

    foreach my $gene (@$gene_seeds) {
        my $cluster_keys = $gene->{cluster_keys};
	my $clusterAlignHash = $gcga->{clusters};
	$gene->{model} = &makeOneGeneModel($gid++, $strand, $cluster_keys, $clusterAlignHash, $self->{debug});
    }
}

# combine ordered genes from both strand
sub _combineOrderedGeneLists {
    my $self = shift;

    my $geneList1 = $self->{plusSeeds};
    my $geneList2 = $self->{minusSeeds};

    my $c = scalar(@$geneList2);
    my $j = 0;
    foreach my $g1 (@$geneList1) {
	my $g1s = $g1->{model}->getGenomicStart;
	my $g1e = $g1->{model}->getGenomicEnd;
	for (; $j < $c; $j++) {
	    my $g2 = $geneList2->[$j];
	    my $g2s = $g2->{model}->getGenomicStart;
	    my $g2e = $g2->{model}->getGenomicEnd;
	    if ($g2s >= $g1s) {
		push @{ $self->{genes} }, $g1;
		last;
	    }
	    push @{ $self->{genes} }, $g2;
	}
	push @{ $self->{genes} }, $g1 if $j == $c;
    }
}

# get clone information common to pairs of gene clusters from files,
# and cache sequence counts for each gene cluster involved
#
sub _getCommonCloneCache {
    my ($self, $key) = @_;

    my $SUB = "_getCommonCloneCache";

    my $ccc = $self->{commonCloneCache}->{$key};
    if ($ccc) {
	print "#debug\t${SUB}::found common clone info in cache\n" if $self->{debug} >= 2;
	return $ccc;
    }

    my $f = $self->{commonCloneCacheFile}->{$key};
    die "File $f not found!" unless -f $f;

    $ccc = { num_seqs_hash => {} };
    print "#debug ${SUB}::reading from file $f\n" if $self->{debug};
    open CCCF, $f;
    my $total = 0;
    while (my $line = <CCCF>) {
	if ($line =~ /^(\d+)\((\d+)\):(\d+):(\d?) (\S+) (\d?):(\d+):\((\d+)\)(\d+)$/) {
	    my ($c1, $n1, $s1, $p1, $c, $p2, $s2, $n2, $c2) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
            if ($c1 >= 0 && $c2 >= 0 and $c1 != $c2) {
	        $total++;
	        if (exists $ccc->{"$c1$c2"}) {
                    push @{ $ccc->{"$c1$c2"} }, "$s1:$p1 $c $s2:$p2";
                } elsif (exists $ccc->{"$c2$c1"}) {
                    push @{ $ccc->{"$c2$c1"} }, "$s2:$p2 $c $s1:$p1";
                } else {
	            $ccc->{"$c1$c2"} = ["$s1:$p1 $c $s2:$p2"];
                }
                $ccc->{num_seqs_hash}->{$c1} = $n1;
                $ccc->{num_seqs_hash}->{$c2} = $n2;
            }
	}
    }
    close CCCF;
    print "#debug ${SUB}::found $total distinct gene cluster pairs with common $key \n"
        if $self->{debug};
    $self->{commonCloneCache}->{$key} = $ccc;

    return $self->{commonCloneCache}->{$key};
}

################ file scoped subroutines ################################

sub reversedTranscript {
  my ($dbh, $bid, $aggressive) = @_;
  my $sql =<<SQL
SELECT count(*) FROM ygan.BlatAlignmentSignals
WHERE blat_alignment_id = $bid
AND ((splice_signal_score_opposite > splice_signal_score and
      splice_signal_score_opposite > 1) or
     (splice_signal_score_opposite > splice_signal_score and
      polya_signal_score_opposite > polya_signal_score))
SQL
;
  if ($aggressive) {
    $sql =<<AGSQL
SELECT count(*) FROM ygan.BlatAlignmentSignals
WHERE blat_alignment_id = $bid
AND (splice_signal_score_opposite > splice_signal_score or
     (splice_signal_score_opposite <= splice_signal_score and
      polya_signal_score_opposite > polya_signal_score));
AGSQL
;
  }
  my $sth = $dbh->prepare($sql) or die "bad sql $sql: $!\n";
  $sth->execute or die "could not execute $sql: $!\n";
  my ($c) = $sth->fetchrow_array;
  $c;
}

sub makeOneGeneModel {
    my ($gid, $strand, $cluster_keys, $clusterAlignHash, $debug) = @_;

    my $SUB = 'makeOneGeneModel';
    print "#debug \t${SUB}::creating gene mode for clusters: ",
          join(',', @$cluster_keys), "\n" if $debug >= 2;
    my $coords = [];
    foreach my $cid (@$cluster_keys) {
	my $alignments = $clusterAlignHash->{$cid};
	# process each genomic alignments
	# TODO maybe refactor to an object
	foreach my $a (@$alignments) {
	    my $num_spans = $a->{number_of_spans};
	    my $tstarts = $a->{tstarts};
	    my $blocksizes = $a->{blocksizes};
	    my @ts = split(/,/, $tstarts);
	    my @bs = split(/,/, $blocksizes);
	    die "Error splitting tstarts $tstarts or blocksizes $blocksizes "
                . "(expecting $num_spans elements)\n" 
		    unless scalar(@ts) == $num_spans and scalar(@bs) == $num_spans;
	    for(my $i=0; $i<$num_spans; $i++) { 
		my @p = ($ts[$i], $ts[$i] + $bs[$i]);
		push @$coords, \@p;
	    }
	}
    }
    return new GeneModel($gid++, $strand, $coords, $debug);
}

sub excludeCommonClusterKeys {
    my ($ori_cks1, $ori_cks2) = @_;
    my (%ck1, %ck2, %c1, %c2, %exc);

    foreach my $k (@$ori_cks1) { 
	$ck1{$k} = "";
	my $c = $k; $c =~ s/\-+|\++//g; $c1{$c} = "";
    }

    foreach my $k (@$ori_cks2) { 
	$ck2{$k} = "";
	my $c = $k; $c =~ s/\-+|\++//g; $c2{$c} = "";
    }

    foreach my $k (keys %ck1) {
	my $c = $k; $c =~ s/\-+|\++//g;
	if (exists $c2{$c}) { delete $ck1{$k}; $exc{$k} = ""; } 
    }

    foreach my $k (keys %ck2) {
	my $c = $k; $c =~ s/\-+|\++//g;
	if (exists $c1{$c}) { delete $ck2{$k}; $exc{$k} = ""; } 
    }

    my @cks1 = keys(%ck1);
    my @cks2 = keys(%ck2);
    my @excs = keys(%exc);
    return (\@cks1, \@cks2, \@excs);
}

sub getUniqueIdsFromClusterKeys {
    my ($c_keys) = @_;

    my %unique_cids;
    foreach my $ck (@$c_keys) {
        my $c = $ck;
        $c =~ s /[\-|\+]//g;
        $unique_cids{$c} = "";
    }
    return keys %unique_cids;
}

1;
