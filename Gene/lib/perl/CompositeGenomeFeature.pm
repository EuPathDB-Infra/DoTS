package DoTS::Gene::CompositeGenomeFeature;
@ISA = qw( DoTS::Gene::GenomeFeature );

use strict;
use DoTS::Gene::GenomeFeature;

sub new {
    my($class, $gf, $constituents) = @_;

    unless ($gf->{coords}) {
	my @coord_sets; foreach (@$constituents) { push @coord_sets, $_->{coords}; }
	$gf->{coords} = &mergeCoordinateSets(@coord_sets);
    }

    my $self = DoTS::Gene::GenomeFeature->new($gf);
    $self->{cs} = $constituents;

    bless $self, $class;
    return $self;
}

sub mergeCoordinateSets {
    my @coord_sets = @_;
    my @coords;
    foreach (@coord_sets) { push @coords, @$_; }

    my $sortedCoords = &_sortCoords(\@coords);
    &_flatSortedCoords($sortedCoords);    
}

sub getConstituents {
    my $self = shift;

    my $genome_id = $self->getGenomeId;
    my $chr = $self->getChrom;
    my $strand = $self->getStrand;

    my @res;
    foreach my $c (@{ $self->{cs} }) {
	my $id = $c->{id};
	my $coords = $c->{coords};
	my $gf_args = { id=>$id, genome_id=>$genome_id, chr=>$chr, coords=>$coords };
	my $gf = DoTS::Gene::GenomeFeature->new($gf_args);
	foreach my $nam (keys %$c) {
	    if ($nam ne 'id' && $nam ne 'coords') {
		$gf->setAnnotationProperty($nam, $c->{$nam});
	    }
	}
	push @res, $gf;
    }
    @res;
}

sub _sortCoords {
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
        print STDERR "#debug \t\t${SUB}::x coords before sort: ";
        foreach (@$coords) { print STDERR $_->[0], ", "; }
        print STDERR "\n";
        print STDERR "#debug \t\t${SUB}::x coords after sort: ", join(', ', @sorted_starts), "\n"; 
    }

    return \@sorted_coords;
}

# TODO implement efficient flattening algorithm (Union-Find or Bit-op?).
sub _flatSortedCoords {
    my ($old_coords, $debug) = @_;

    my $SUB = 'flatSortedCoords';
    my $old_tot = scalar(@$old_coords);

    if ($debug >= 3) {
	print STDERR "#debug \t\t${SUB}::content of old_coords: ";
	foreach my $oc (@$old_coords) { print STDERR "(", join(', ', @$oc), "), "; }
	print STDERR "\n";
    }

    for (my $i=0; $i<$old_tot-1; $i++) {
	my ($s1,$e1) = @{ $old_coords->[$i] };
	unless (($s1 || $s1 == 0) && $e1) {
	    print STDERR "#WARNING \t\t${SUB}: ($s1, $e1) at $i!\n" if $debug >= 3;
	    next;
	}

	for (my $j=$i+1; $j<$old_tot; $j++) {
	    ($s1, $e1) = @{ $old_coords->[$i] };
	    my ($s2,$e2) = @{ $old_coords->[$j] };

	    print STDERR "#debug \t\t\t${SUB}::examine join between ($s1,$e1) at $i with ($s2,$e2) at $j\n" 
		if $debug >= 4;

	    die "coordinate must be none descending by x" if $s1 > $s2;

	    my $tiny_gap = 5;
	    if ($s2 <= $e1 + $tiny_gap) {
		my ($s, $e) = ($s1, $e1);
		$s = $s2 if $s2 < $s1;
		$e = $e2 if $e2 > $e1;
		$old_coords->[$i] = [$s, $e];
		$old_coords->[$j] = [];
		print STDERR "#debug \t\t${SUB}::merging coord at $j to $i\n" if $debug >= 3;
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
		print STDERR "#WARNING\t${SUB}::ignore unexpected exon coordinate ($s, $e)\n";
	    }
	}
    }

    if ($debug >= 3) {
        print STDERR "#debug \t\t${SUB}::content of new_coords: ";
        foreach my $nc (@$new_coords) { print STDERR "(", join(', ', @$nc), "), "; }
        print STDERR "\n";
    }

    return $new_coords;
}

1;
