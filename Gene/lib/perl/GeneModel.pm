#!/usr/bin/perl

# -------------------------------------------------------
# Gene Model from a set of genomic alignment coordinates
#
# Yongchang Gan July 28, 2002
# -------------------------------------------------------

package GeneModel;

use strict;

# constructor
#
sub new {
    my($class, $id, $strand, $coords, $debug) = @_;

    my $self = { id => $id, strand => $strand, };
    my $gene = &createGeneModel($coords, $debug);
    $self->{model} = $gene->{model};
    $self->{model_string} = $gene->{model_string};
    $self->{gene_size} = $gene->{gene_size};
    $self->{genomic_start} = $gene->{genomic_start};
    $self->{genomic_end} = $gene->{genomic_end};
    $self->{genomic_size} = $gene->{genomic_size};
    $self->{min_exon_size} = $gene->{min_exon_size};
    $self->{max_exon_size} = $gene->{max_exon_size};
    $self->{min_intron_size} = $gene->{min_intron_size};
    $self->{max_intron_size} = $gene->{max_intron_size};
    die "# ERROR: GeneModel::new: no gene size or genomic_end for initial coords: ",
      join(', ', map { $_->[0] . '-' . $_->[1] } @$coords), "\n" unless $self->{gene_size} && $self->{genomic_end};

    bless $self, $class;
    return $self;
}

# subroutine
#
sub getId {
    my $self = shift;
    return $self->{id};
}

sub getStrand {
    my $self = shift;
    return $self->{strand};
}

sub getGeneModel {
    my $self = shift;
    return $self->{model};
}

sub getExonCoords {
    my $self = shift;
    return $self->getGeneModel;
}

sub getNumberOfExons {
    my $self = shift;
    return scalar(@{ $self->getExonCoords });
}

sub getMinExonSize {
    my $self = shift;
    return $self->{min_exon_size};
}

sub getMaxExonSize {
    my $self = shift;
    return $self->{max_exon_size};
}

sub getMinIntronSize {
    my $self = shift;
    return $self->{min_intron_size};
}

sub getMaxIntronSize {
    my $self = shift;
    return $self->{max_intron_size};
}

sub getModelString {
    my $self = shift;
    return $self->{model_string};
}

sub getGeneSize {
    my $self = shift;
    return $self->{gene_size};
}

sub getGenomicStart {
    my $self = shift;
    return $self->{genomic_start};
}

sub getGenomicEnd {
    my $self = shift;
    return $self->{genomic_end};
}

sub getGenomicSize {
    my $self = shift;
    return $self->{genomic_end} - $self->{genomic_start};
}
   
######## file scoped subroutines ###################
sub createGeneModel {
    my ($coords, $debug) = @_;
    my $sortedCoords = &sortCoords($coords, $debug);
    my $flatSorted = &flatSortedCoords($sortedCoords, $debug);
    my ($gene_size, $gStart, $gEnd, $gSize);

    my $s0 = $flatSorted->[0]->[0];
    my $e0 = $flatSorted->[0]->[1];
    die "# ERROR: GeneModel::createGeneModel: bad model: " 
      . join(', ', map { $_->[0] . '-' . $_->[1] } @$flatSorted) . "(original coords="
      . join(', ', map { $_->[0] . '-' . $_->[1] } @$coords) . "; sorted coords="
      . join(', ', map { $_->[0] . '-' . $_->[1] } @$sortedCoords) . ")\n" unless $e0;
    my $es = $e0 - $s0;
    $gene_size += $es;

    my $model_string = $s0 . '-'. $e0 . ' .. ';
    my $minE = $es;
    my $maxE = $es;
    my ($minI, $maxI) = (0,0);
    my $c = scalar(@$flatSorted);
    for (my $i=1; $i<$c; $i++) {
	my ($s, $e) = @{ $flatSorted->[$i] };
	$es = $e - $s;
        my $is = $s - $flatSorted->[$i-1]->[1];

	$model_string .= $s . '-' . $e . ' .. ';
	$gene_size += $es;
	$minE = $es if $es < $minE;
	$maxE = $es if $es > $maxE;
	$minI = $is if $is < $minI or $minI == 0;
	$maxI = $is if $is > $maxI;
    }
    $model_string =~ s/ .. $//;
    $gStart = $flatSorted->[0]->[0];
    $gEnd = $flatSorted->[scalar(@$flatSorted)-1]->[1];
    $gSize = $gEnd - $gStart;

    return { model => $flatSorted, model_string => $model_string,
             gene_size => $gene_size, genomic_start => $gStart,
             genomic_end => $gEnd, genomic_size => $gSize, 
	     min_exon_size => $minE, max_exon_size => $maxE,
	     min_intron_size => $minI, max_intron_size => $maxI,
	 };
}

sub sortCoords {
    my ($coords, $debug) = @_;

    my $SUB = 'sortCoords';
    my %coord_hash;
    foreach my $c (@$coords) { 
        my ($s, $e) = @$c;
        if (exists $coord_hash{$s}) {
            push @{ $coord_hash{$s} }, $e;
        } else {
            $coord_hash{$s} = [$e];
        }
    }
    my @sorted_starts = sort { $a <=> $b } keys %coord_hash;

    my @sorted_coords;
    foreach my $ss (@sorted_starts) {
        my @es = @{ $coord_hash{$ss} };
        my @ses = sort { $a <=> $b } @es;
        foreach my $se (@ses) { push @sorted_coords, [$ss, $se]; }
    }

    if ($debug == 3) {
        print "#debug \t\t${SUB}::x coords before sort: ";
        foreach (@$coords) { print $_->[0], ", "; }
        print "\n";
        print "#debug \t\t${SUB}::x coords after sort: ", join(', ', @sorted_starts), "\n"; 
    }

    return \@sorted_coords;
}

# TODO implement efficient flattening algorithm (Union-Find or Bit-op?).
sub flatSortedCoords {
    my ($old_coords, $debug) = @_;

    my $SUB = 'flatSortedCoords';
    my $old_tot = scalar(@$old_coords);

    if ($debug >= 3) {
	print "#debug \t\t${SUB}::content of old_coords: ";
	foreach my $oc (@$old_coords) { print "(", join(', ', @$oc), "), "; }
	print "\n";
    }

    for (my $i=0; $i<$old_tot-1; $i++) {
	my ($s1,$e1) = @{ $old_coords->[$i] };
	next unless defined($s1) and defined($e1);
	unless (($s1 || $s1 == 0) && $e1) {
	    print "#WARNING \t\t${SUB}: ($s1, $e1) at $i!\n" if $debug >= 3;
	    next;
	}

	for (my $j=$i+1; $j<$old_tot; $j++) {
	    ($s1, $e1) = @{ $old_coords->[$i] };
	    my ($s2,$e2) = @{ $old_coords->[$j] };

	    print "#debug \t\t\t${SUB}::examine join between ($s1,$e1) at $i with ($s2,$e2) at $j\n" 
		if $debug >= 4;

	    die "coordinate must be none descending by x" if $s1 > $s2;

	    my $tiny_gap = 4;
	    if ($s2 <= $e1 + $tiny_gap) {
		my ($s, $e) = ($s1, $e1);
		$s = $s2 if $s2 < $s1;
		$e = $e2 if $e2 > $e1;
		$old_coords->[$i] = [$s, $e];
		$old_coords->[$j] = [];
		print "#debug \t\t${SUB}::merging coord at $j to $i\n" if $debug >= 3;
	    } else { last; }
	}
    }

    my $new_coords = [];
    foreach my $coord (@$old_coords) {
	my ($s, $e) = @$coord;
	if (defined($s) && $e) {
	    if ($s >= 0 && $e > 0 && $e > $s) { 
		push @$new_coords, $coord;
	    } else {
		print "#WARNING\t${SUB}::ignore unexpected exon coordinate ($s, $e)\n";
	    }
	}
    }

    if ($debug >= 3) {
        print "#debug \t\t${SUB}::content of new_coords: ";
        foreach my $nc (@$new_coords) { print "(", join(', ', @$nc), "), "; }
        print "\n";
    }

    return $new_coords;
}

1;
