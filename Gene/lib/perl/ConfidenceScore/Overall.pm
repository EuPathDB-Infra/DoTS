package DoTS::Gene::ConfidenceScore::Overall;

use strict;
1;

sub setScore {
    my ($dbh, $gdg_tab, $gid, $cs_info, $stats, $debug) = @_;

    my ($sql, $sth);

    # info already in the gene table
    $sql = "select max_intron, number_of_exons, deprecated "
	 . "from $gdg_tab where genome_dots_gene_id = $gid";
    $sth = $dbh->prepare($sql) or die "bad sql $sql:$!\n";
    $sth->execute or die "could not run $sql:$!\n";

    my ($max_intron, $num_exons, $dep);
    unless (($max_intron, $num_exons) = $sth->fetchrow_array) {
	die "no information found for gene $gid!\n";
    }
    $sth->finish;

    my $score = 0;
    # see top of script for calculation of overall score
    if ($max_intron >= 47) {
	$score += 1; $stats->{spliced}++; }
    if ($num_exons >= 3) {
	$score += 1; $stats->{exon3}++; }
    if (abs($cs_info->{number_of_splice_signals}) >= 1) {
	$score += 1; $stats->{splice_signals1}++; }
    if (abs($cs_info->{number_of_splice_signals}) >= 2) {
	$score += 1; $stats->{splice_signals2}++; }
    if ($cs_info->{contains_mrna}) {
	$score += 3; $stats->{contains_mrna}++; }
    if ($cs_info->{number_of_rnas} >= 2) {
	$score += 1; $stats->{rnas}++; }
    if ($cs_info->{number_of_est_libraries} >= 2) {
	$score += 1; $stats->{est_libs}++; }
    if ($cs_info->{number_of_est_clones} >= 2) {
	$score += 1; $stats->{est_clones}++; }
    if ($cs_info->{number_of_est_p53pairs} >= 1) {
	$score += 1; $stats->{est_p53pairs}++; }
    if (abs($cs_info->{polya_signal_type}) == 2) {
	$score += 2; $stats->{polya2}++; }
    if (abs($cs_info->{polya_signal_type}) == 1) {
	$score += 1; $stats->{polya1}++; }
    if (abs($cs_info->{est_plot_score}) > 0) {
	$score += 1; $stats->{plot_score}++ }
    if ($cs_info->{est_plot_score} == 100) {
	$score += 1; $stats->{plot_score100}++; }
    if ($cs_info->{est_plot_score} == -100) {
	$score -= 1; $stats->{plot_score_neg100}++; }

    if (abs($cs_info->{has_polya_track}) == 1) {
	$score = -$score; $stats->{polya_track}++;
    } elsif ($dep) {
        $score = -$score; $stats->{anti_sense}++;
    }

    $cs_info->{confidence_score} = $score;
}
