package DoTS::Gene::ConfidenceScore::Composition;

use strict;

1;

=item setComposition:

set input EST/mRAN composition info for given gDG

=cut

sub setComposition {
    my ($dbh, $gdg_tab, $gdt_tab, $gid, $cs_info, $debug) = @_;

    # get number of EST clones and number of 5'-3' pairs
    my ($libs, $clones, $pairs) = &_getEstLibCloneAndPairCounts($dbh, $gdg_tab, $gdt_tab, $gid);
    $cs_info->{number_of_est_libraries} = $libs;
    $cs_info->{number_of_est_clones} = $clones;
    $cs_info->{number_of_est_p53pairs} = $pairs;

    # get number of RNAs and whether contains mRNA
    my ($num_rnas, $contains_mrna) = &_getRNAInfo($dbh, $gdt_tab, $gid);
    $cs_info->{number_of_rnas} = $num_rnas;
    $cs_info->{contains_mrna} = $contains_mrna;
}

sub _getRNAInfo {
  my ($dbh, $gdt_tab, $gid) = @_;
  my $sql = "select count(*), max(a.contains_mrna) "
	 . "from $gdt_tab gdt, DoTS.Assembly a "
	 . "where gdt.genome_dots_gene_id = $gid "
         . "and gdt.na_sequence_id = a.na_sequence_id ";
  my $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
  $sth->execute or die "could not run $sql:$!\n";
  my ($rnas, $has_mrna) = (0, 0);
  unless (($rnas, $has_mrna) = $sth->fetchrow_array) {
    die "no information found for gene $gid!\n";
  }
  $has_mrna = 0 unless $has_mrna;
  ($rnas, $has_mrna);
}

sub _getEstLibCloneAndPairCounts {
    my ($dbh, $gdg_tab, $gdt_tab, $gid) = @_;

    my $sql = "select c.library_id, c.clone_id, "
         . "c.washu_name, c.image_id, est.na_sequence_id, est.p_end "
         . "from $gdg_tab gdg, $gdt_tab gdt, \n"
         . "DoTS.AssemblySequence aseq, DoTS.EST est, DoTS.Clone c \n"
         . "where gdg.genome_dots_gene_id = $gid \n"
         . "and gdg.genome_dots_gene_id = gdt.genome_dots_gene_id \n"
         . "and gdt.na_sequence_id = aseq.assembly_na_sequence_id \n"
         . "and aseq.na_sequence_id = est.na_sequence_id \n"
         . "and est.clone_id = c.clone_id";

    my $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
    $sth->execute or die "could not run $sql:$!\n";

    my (@rows, %libs, %clones);
    while (my @row = $sth->fetchrow_array) {
        push @rows, \@row;
        my ($lid, $cid) = @row;
        $libs{$lid} = '';
        $clones{"$lid:$cid"} = '';
    }
    $sth->finish;

    my $libs = scalar(keys %libs);
    my $clones = scalar(keys %clones);

    my $pairs = 0;
    my $num_rows = scalar(@rows);
    for (my $i=0; $i<$num_rows-1; $i++) {
        my ($l1, $c1, $w1, $i1, $e1, $p1) = @{ $rows[$i] };
        for (my $j=$i+1; $j<$num_rows; $j++) {
            my ($l2, $c2, $w2, $i2, $e2, $p2) = @{ $rows[$j] };
            my $sameLib = ($l1 == $l2);
            my $sameClone = ($w1 eq $w2 && $w1) || ($i1 eq $i2 && $i1);
            my $pairEnds = ($p1 eq '5' && $p2 eq '3') || ($p1 eq '3' && $p2 eq '5');
            $pairs++ if $sameLib && $sameClone && $pairEnds;
        }
    }

    ($libs, $clones, $pairs);
}
