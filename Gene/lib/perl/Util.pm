package DoTS::Gene::Util;

use strict;

sub reverseComplement {
    my ($seq) = @_;
    $seq =~ tr/ATCGatcg/TAGCtagc/;
    return reverse $seq;
}

sub arraySum {
  my @es = @_;
  my $res;
  foreach (@es) { $res += $_; }
  return $res;
}

sub sharedHashKeys {
    my ($h1, $h2) = @_;
    my %ks;
    foreach (keys %$h1) {
	$ks{$_} = 1 if exists $h2->{$_};
    }
    return %ks;
}

sub getAvgStddev {
  my @es = @_;
  my $tot = scalar(@es);

  (undef, undef) if $tot < 1;
  return ($es[0], 0) if $tot < 2;

  my ($sum, $sum_sq) = (0, 0);
  foreach my $e (@es) { $sum += $e; $sum_sq += $e**2; }
  my $avg = $sum / $tot;
  my $stddev = sqrt(($tot * $sum_sq - $sum**2)/($tot*($tot-1)));

  ($avg, $stddev);
}

1;

