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

sub getChromIdAndLen {
  my ($dbh, $ext_db_rel_id, $chr) = @_;

  my $sql = "select na_sequence_id, length(sequence) from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id and chromosome = '$chr'";
  my $sth = $dbh->prepareAndExecute($sql);
  my ($chr_id, $len) = $sth->fetchrow_array;

  die "no chr_id found for $chr and external database release $ext_db_rel_id\n" unless $chr_id;
  ($chr_id,$len);
}

sub getChromsInfo {
  my ($dbh, $ext_db_rel_id, $skipChrs) = @_;

  my $sql = "select na_sequence_id, chromosome, length(sequence) from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id order by chromosome_order_num";
  if (scalar(@$skipChrs) > 0) {   
    $sql .= " and chromosome not in (" . join(', ', map { "'" . $_ . "'" } @$skipChrs). ")";
  }
  my $sth = $dbh->prepareAndExecute($sql);

  my @chroms;
  while (my ($chr_id, $chr, $len) = $sth->fetchrow_array) {
      push @chroms, { chr_id => $chr_id, chr => $chr, len => $len, start => 0, end => $len };
  }

  return @chroms;
}

sub getCoordSelectAndSkip {
    my ($dbh, $genomeId, $opt) = @_;

    if ($opt->{test}) {
	my $chr = '1';
	my ($chr_id) = &getChromIdAndLen($dbh, $genomeId, $chr);
	return [{chr=>$chr, chr_id=>$chr_id, start=>5e6, end=>10e6}];
    }

    my $chr = $opt->{chr}; $chr =~ s/^chr//i;
    my $start = $opt->{start};
    my $end = $opt->{end};

    die "start or end set while chr is not defined" if (defined($start) || $end) && !$chr;
    if ($chr) {
	my ($chr_id, $chr_len) = &getChromIdAndLen($dbh, $genomeId, $chr);
	my $s = 0 unless $start; 
	my $e = $chr_len unless $end < $chr_len;
	return ([{ chr_id=> $chr_id, chr => $chr, start => $s, end => $e }], []);
    } else {
	my $skipChrs = $opt->{skip_chrs}; $skipChrs =~ s/chr//gi;
	my @skipChrs = split(/,/, $skipChrs);
	my @chroms = &getChromsInfo($dbh, $genomeId, \@$skipChrs);
	return (\@chroms, \@skipChrs);
    }
}

1;

