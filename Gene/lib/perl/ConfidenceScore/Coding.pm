package DoTS::Gene::ConfidenceScore::Coding;

use strict;

1;

=item setCoding:

set coding potential attributes for given gDG

=cut

sub setCoding {
  my ($dbh, $gdg_tab, $gdt_tab, $gid, $cs_info, $debug) = @_;

  my $sql = "select max(tas.length), max(taf.translation_score), min(taf.p_value) "
    . "from $gdg_tab gdg, $gdt_tab gdt, dots.NAFeature naf, "
    . "dots.TranslatedAAFeature taf, dots.TranslatedAASequence tas "
    . "where gdg.genome_dots_gene_id = $gid "
    . "and gdg.genome_dots_gene_id = gdt.genome_dots_gene_id "
    . "and gdt.na_sequence_id = naf.na_sequence_id "
    . "and naf.na_feature_id = taf.na_feature_id "
    . "and taf.aa_sequence_id = tas.aa_sequence_id";

  my $sth = $dbh->prepare($sql) or die "bad sql $sql: $!";
  $sth->execute or die "could not run $sql: $!";

  my ($max_orf_len, $max_orf_score, $min_orf_pval);
  if (($max_orf_len, $max_orf_score, $min_orf_pval) = $sth->fetchrow_array) {}

  $cs_info->{max_orf_length} = ($max_orf_len ? $max_orf_len : 0);
  $cs_info->{max_orf_score} = ($max_orf_score ? $max_orf_score : 0);
  $cs_info->{min_orf_pval} = ($min_orf_pval ? $min_orf_pval : 9999);
}
