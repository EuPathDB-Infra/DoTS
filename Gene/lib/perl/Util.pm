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

sub getAnalysisId {
    my ($dbh, $taxonId, $genomeId) = @_;

    my $sql = "select aligned_gene_analysis_id from Allgenes.AlignedGeneAnalysis "
	. "where parameters like '%--t $taxonId --dbr %' "
	. "order by aligned_gene_analysis_id desc";
    my $sth = $dbh->prepareAndExecute($sql);
    my @aids = ();
    if (my $aid = $sth->fetchrow_array) { push @aids, $aid; }
    my $c = scalar(@aids);
    die("expecting one analysis id but found $c") unless $c == 1;

    return $aids[0];
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

sub getChromToId {
  my @chroms = &getChromsInfo(@_);

  my $res = {};
  foreach (@chroms) {
      my ($c, $cid) = ($_->{chr}, $_->{chr_id});
      $res->{$c} = $cid;
  }
  return $res;
}

sub getIdToChrom {
  my @chroms = &getChromsInfo(@_);

  my $res = {};
  foreach (@chroms) {
      my ($c, $cid) = ($_->{chr}, $_->{chr_id});
      $res->{$cid} = $c;
  }
  return $res;
}


sub getChromsInfo {
  my ($dbh, $ext_db_rel_id, $skipChrs) = @_;

  my $sql = "select na_sequence_id, chromosome, length(sequence) from DoTS.VirtualSequence "
    . "where external_database_release_id = $ext_db_rel_id";
  if ($skipChrs && scalar(@$skipChrs) > 0) {   
    $sql .= " and chromosome not in (" . join(', ', map { "'" . $_ . "'" } @$skipChrs). ")";
  }

  $sql .= " order by chromosome_order_num";

  my $sth = $dbh->prepareAndExecute($sql);

  my @chroms;
  while (my ($chr_id, $chr, $len) = $sth->fetchrow_array) {
      push @chroms, { chr_id => $chr_id, chr => $chr, len => $len, start => 0, end => $len };
  }

  return @chroms;
}

sub getCoordSelectAndSkip {
    my ($dbh, $genomeId, $opt) = @_;

    my $chr = $opt->{chr}; $chr =~ s/^chr//i;
    my $start = $opt->{start};
    my $end = $opt->{end};

    if ($opt->{test}) {	
	$chr = '1' unless $chr;
	$start = 5e6 unless $start;
	$end = 10e6 unless $end;
	my ($chr_id) = &getChromIdAndLen($dbh, $genomeId, $chr);
	return [{chr=>$chr, chr_id=>$chr_id, start=>$start, end=>$end}];
    }

    die "start or end set while chr is not defined" if (defined($start) || $end) && !$chr;
    if ($chr) {
	my ($chr_id, $chr_len) = &getChromIdAndLen($dbh, $genomeId, $chr);
	$start = 0 unless $start; 
	$end = $chr_len unless $end && $end < $chr_len;
	return ([{ chr_id=> $chr_id, chr => $chr, start => $start, end => $end }], []);
    } else {
	my $skipChrs = $opt->{skip_chrs};
	if(ref($skipChrs) eq 'ARRAY') { $skipChrs = join(',', @$skipChrs); }
        $skipChrs =~ s/chr//gi;
	my @skipChrs = split(/,/, $skipChrs);
	my @chroms = &getChromsInfo($dbh, $genomeId, \@skipChrs);
	return (\@chroms, \@skipChrs);
    }
}

sub delete {
    my ($dbh, $tab, $idCol, $ids, $batchSize) = @_;

    $batchSize = 1000 unless $batchSize;
    my @ids = @{ $ids };

    my $tot = scalar(@ids);
    my $ct = 0;
    for(my $i=0;$i<$tot;$i += $batchSize){
	my $j = ($i + $batchSize - 1 < $tot ? $i + $batchSize - 1 : $tot - 1);
	$ct += $dbh->do("delete from $tab where $idCol in (".join(', ',@ids[$i..$j]).")");
	$dbh->commit();
    }
    if ($ct != $tot) { 
	print STDERR "WARNING: only deleted $ct entries from $tab for $tot $idCol\n";
    }
    return $ct;
}

1;

