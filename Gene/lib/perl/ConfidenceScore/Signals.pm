package DoTS::Gene::ConfidenceScore::Signals;

use strict;

1;

=item setSignals:

set splice signal, polya signal, polya track info for given gDG

=cut

sub setSignals {
    my ($dbh, $gdt_tab, $bs_tab, $gid, $cs_info, $debug) = @_;

    # Simplistic approach for polya signal/polya track, but should be ok for splice signal
    #
    my $sql = "select max(bs.splice_signal_score), max(bs.splice_signal_score_opposite), "
            . " max(bs.polya_signal_score), max(bs.polya_signal_score_opposite), "
            . " max(bs.has_polya_track), max(bs.has_polya_track_opposite) "
            . "from $bs_tab bs, $gdt_tab gdt "
            . "where bs.blat_alignment_id = gdt.blat_alignment_id "
            . "and gdt.genome_dots_gene_id = $gid";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
    $sth->execute or die "could not execute $sql: $!\n";
    my ($spp, $spm, $pap, $pam, $hpp, $hpm) = $sth->fetchrow_array;

    my ($splsig, $pas_type, $has_pat) = (0, 0, 0);
    if ($spp > $spm) { $splsig = $spp; }
    if ($spm > $spp) { $splsig = 0-$spm; }
    if ($pap > $pam) { $pas_type = $pap; }
    if ($pam > $pap) { $pas_type = 0-$pam; }
    if ($hpp) { $has_pat = 1; } elsif ($hpm) { $has_pat = -1; } else { $has_pat = 0; }

    $cs_info->{number_of_splice_signals} = $splsig;
    $cs_info->{polya_signal_type} = $pas_type;
    $cs_info->{has_polya_track} = $has_pat;
}
