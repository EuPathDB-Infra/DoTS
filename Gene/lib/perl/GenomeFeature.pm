#!/usr/bin/perl

package DoTS::Gene::GenomeFeature;

use strict;
use DoTS::Gene::Util;
use Carp;

# constructor
#
sub new {
    my($class, $gf) = @_;

    my $coords = $gf->{coords};
    confess "required coordinate set is not found" unless defined $coords;
    my $self = { id => $gf->{id}, genome_id => $gf->{genome_id},
		 chr=> $gf->{chr}, strand => $gf->{strand} };
    my $feat = &createFeature($coords);
    $self->{model} = $feat->{model};
    $self->{model_string} = $feat->{model_string};
    $self->{total_span_size} = $feat->{total_span_size};
    $self->{chrom_start} = $feat->{chrom_start};
    $self->{chrom_end} = $feat->{chrom_end};
    $self->{genomic_size} = $feat->{genomic_size};
    $self->{min_span_size} = $feat->{min_span_size};
    $self->{max_span_size} = $feat->{max_span_size};
    $self->{min_interspan_size} = $feat->{min_interspan_size};
    $self->{max_interspan_size} = $feat->{max_interspan_size};
    $self->{number_of_spans} = $feat->{number_of_spans};
    die "# GeneFeature::new: bad coords: " . join(',', map { $_->[0] . '-' . $_->[1] } @$coords)
        unless $self->{total_span_size} && $self->{chrom_end};
    $self->{ap} = {};

    bless $self, $class;
    return $self;
}

# subroutine
#
sub getId { $_[0]->{id}; }

sub getGenomeId { $_[0]->{genome_id}; }

sub getChrom { $_[0]->{chr}; }

sub getStrand { $_[0]->{strand}; }

sub getModel { $_[0]->{model}; }

sub getSpanCoordinates { $_[0]->getModel; }

sub getSpanStarts {
    my $coords = $_[0]->getSpanCoordinates();
    return map { $_->[0] } @$coords;
}

sub getSpanEnds {
    my $coords = $_[0]->getSpanCoordinates();
    return map { $_->[1] } @$coords;
}

sub getInterspanCoordinates {
  my $self = shift;
  my $spans = $self->getSpanCoordinates;
  my $sc = $self->getNumberOfSpans;
  my @interspans;
  for (my $i=0; $i<$sc-1; $i++) {
    my $is = $spans->[$i]->[1];
    my $ie = $spans->[$i+1]->[0];
    push @interspans, [$is, $ie];
  }
  \@interspans;
}

sub getInterspanStarts {
    my $coords = $_[0]->getInterspanCoordinates();
    return map { $_->[0] } @$coords;
}

sub getInterspanEnds {
    my $coords = $_[0]->getInterspanCoordinates();
    return map { $_->[1] } @$coords;
}

sub getNumberOfSpans { $_[0]->{number_of_spans}; }

sub getMinSpanSize { $_[0]->{min_span_size}; }

sub getMaxSpanSize { $_[0]->{max_span_size}; }

sub getMinInterspanSize { $_[0]->{min_interspan_size}; }

sub getMaxInterspanSize { $_[0]->{max_interspan_size}; }

sub getModelString { $_[0]->{model_string}; }

sub getTotalSpanSize { $_[0]->{total_span_size}; }

sub getChromStart { $_[0]->{chrom_start}; }

sub getChromEnd { $_[0]->{chrom_end}; }

sub getGenomicBoundaries { ($_[0]->getChromStart, $_[0]->getChromEnd); }

sub getGenomicSize {
    my $self = shift;
    return $self->getChromEnd - $self->getChromStart;
}

sub getSpanSequences {
  my ($self, $dbh) = @_;
  my $coords = $self->getSpanCoordinates;
  $self->getSeqs($dbh, $coords);
}

sub getIntronSequences {
  my ($self, $dbh) = @_;
  my $coords = $self->getInterspanCoordinates;
  $self->getSeqs($dbh, $coords);
}

sub getFlankSequence5 {
  my ($self, $dbh, $fi, $fo) = @_;
  $self->getFlankSeq($dbh, $fi, $fo, 5);
}

sub getFlankSequence3 {
  my ($self, $dbh, $fi, $fo) = @_;
  $self->getFlankSeq($dbh, $fi, $fo, 3);
}

sub getSeqs {
  my ($self, $dbh, $coords) = @_;

  my $genome_id = $self->getGenomeId;
  my $chr = $self->getChrom;
  my $str = $self->getStrand;
  my @seqs;
  foreach (@$coords) {
    my ($cs, $ce) = @$_;
    my ($s, $len) = ($cs+1, $ce - $cs);
    die "bad coordinate [$cs, $ce]\n" if $len < 0;
    if ($len == 0) { push @seqs, '';  next; }
    my $sql = "select substr(sequence, $s, $len) "
            . "from DoTS.VirtualSequence "
	    . "where external_database_release_id = $genome_id "
	    . "and chromosome = '$chr'";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
    $sth->execute or die "could not run $sql:!\n";
    my $seq;
    if (($seq) = $sth->fetchrow_array) { ; }
    die "could not get $len bp from chr$chr:$s\n" unless $seq;
    if ($str =~ /\-/) { $seq = &DoTS::Gene::Util::reverseComplement($seq); }
    push @seqs, $seq;
    $sth->finish;
  }
  @seqs = reverse @seqs if $str =~ /\-/;
  \@seqs;
}

sub getFlankSeq {
  my ($self, $dbh, $fi, $fo, $p) = @_;

  die "flank offsets must be >=0" if $fi < 0 || $fo < 0;

  my $genome_id = $self->getGenomeId;
  my $chr = $self->getChrom;
  my $str = $self->getStrand;
  my $gs = $self->getChromStart;
  my $ge = $self->getChromEnd;
  my $isPlus = ($str =~ /\+/);
  my $isP5 = ($p == 5);
  my $s;
  $s = $gs - $fo if ($isPlus && $isP5) or (!$isPlus && !$isP5);
  $s = $ge - $fi if (!$isPlus && $isP5) or ($isPlus && !$isP5);
  $s = 0 if $s < 0;
  my $len = $fo+$fi;
  return undef if $len < 1;

  my $sql = "select substr(sequence, $s, $len) "
          . "from DoTS.VirtualSequence "
          . "where external_database_release_id = $genome_id "
          . "and chromosome = '$chr'";
  my $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
  $sth->execute or die "could not run $sql:!\n";
  my $seq;
  if (($seq) = $sth->fetchrow_array) { ; }
  die "could not get $len bp from chr$chr:$s\n" unless $seq;
  unless ($isPlus) { $seq = &DoTS::Gene::Util::reverseComplement($seq); }
  $sth->finish;
  $seq;
}

sub getAnnotationProperties { $_[0]->{ap}; }

sub setAnnotationProperty {
    my ($self, $nam, $val) = @_;
    $self->{ap}->{$nam} = $val;
}

sub getSpanOverlapWith {
    my ($self, $coords) = @_;
    &getSpanOverlap($self->getSpanCoordinates(), $coords);
}

######## file scoped subroutines ###################

sub createFeature {
    my ($coords) = @_;

    # assume the coords contain sorted non-overlap spans!

    my ($total_span_size, $gStart, $gEnd, $gSize);

    my $s0 = $coords->[0]->[0];
    my $e0 = $coords->[0]->[1];
    my $ss = $e0 - $s0;
    $total_span_size += $ss;

    my $model_string = $s0 . '-'. $e0 . ' .. ';
    my $minS = $ss;
    my $maxS = $ss;
    my ($minI, $maxI) = (0,0);
    my $c = scalar(@$coords);
    for (my $i=1; $i<$c; $i++) {
	my ($s, $e) = @{ $coords->[$i] };
	$ss = $e - $s;
        my $is = $s - $coords->[$i-1]->[1];

	$model_string .= $s . '-' . $e . ' .. ';
	$total_span_size += $ss;
	$minS = $ss if $ss < $minS;
	$maxS = $ss if $ss > $maxS;
	$minI = $is if $is < $minI or $minI == 0;
	$maxI = $is if $is > $maxI;
    }
    $model_string =~ s/ .. $//;
    $gStart = $coords->[0]->[0];
    $gEnd = $coords->[scalar(@$coords)-1]->[1];
    $gSize = $gEnd - $gStart;

    return { model => $coords, model_string => $model_string, number_of_spans => $c,
             total_span_size => $total_span_size, chrom_start => $gStart,
             chrom_end => $gEnd, genomic_size => $gSize, 
	     min_span_size => $minS, max_span_size => $maxS,
	     min_interspan_size => $minI, max_interspan_size => $maxI,
	 };
}

sub getSpanOverlap {
    my ($coords1, $coords2) = @_;

    my $overlap = 0;
    my $c1 = scalar(@$coords1);
    my $c2 = scalar(@$coords2);

    my $lastJ = 0;
    for(my $i=0; $i<$c1; $i++) {
	my ($s1, $e1) = @{ $coords1->[$i] };
	for(my $j = $lastJ; $j<$c2; $j++) {
	    my ($s2, $e2) = @{ $coords2->[$j] };

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

	    $overlap += $o;
	}
    }

    return $overlap;
}

1;
