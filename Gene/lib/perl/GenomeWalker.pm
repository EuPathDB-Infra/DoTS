package DoTS::Gene::GenomeWalker;

use strict;
use DoTS::Gene::Util;
use DoTS::Gene::GenomeAlignmentSelector;
use DoTS::Gene::GenomeAlignmentMerger;
use DoTS::Gene::GenomeFeature;
use DoTS::Gene::CompositeGenomeFeature;

# constructor
#
sub new {
    my($class, $db, $select_criteria, $merge_criteria) = @_;

    my $self = { db => $db, sc => $select_criteria, mc => $merge_criteria };

    bless $self, $class;
    return $self;
}

# get merged alignments in this genomic region,
# return array ref of DoTS::Gene::CompositeGenomeFeature
#
sub getCompositeGenomeFeatures {
    my $self = shift;
    return $self->{cgf} if defined $self->{cgf};

    $self->{cgf} = $self->walkGenome();
    return $self->{cgf};
}

#####################################################################################

sub walkGenome {
    my $self = shift;
    my $db = $self->{db}; my $dbh = $db->getQueryHandle();

    my $sc = $self->{sc};
    my ($bq, $bt, $oq, $ot) = ($sc->{baseQ}, $sc->{baseT}, $sc->{optQ}, $sc->{optT});

    # TODO: break into bins softly (with no composite genome feature crossing boundary)
    #       so that each bin can be processed independently
    my @bins = &_breakIntoBins($ot->{start}, $ot->{end}, 10e6);

    my @bin_res;
    foreach (@bins) {
	my ($s, $e) = @$_; my $otb = { start=>$s, end=>$e };
	my $aln_sel = DoTS::Gene::GenomeAlignmentSelector->new($db,$bq,$bt,$oq,$otb);
	my $sas = $aln_sel->getSortedAlignments();
	push @bin_res, $sas;
    }

    my $srt_alns = [];
    my $qRevs = &_getQueryReversed($dbh, $bt->{target_external_db_release_id}, $oq->{blat_signals_cache});
    foreach (@bin_res) {
	foreach my $aln (@$_) {
	    my $id = $aln->getBlatAlignmentId();
	    $aln->retrieveFromDB();
	    my $q_seq_id = $aln->getQueryNaSequenceId();
	    if (!$q_seq_id) { warn("query seq id not found for blat alignment $id, skip"); next; }

	    my $isRev = $aln->getIsReversed();
	    $aln->setIsReversed(1-$isRev) if $qRevs->{$id};
	    push @$srt_alns, $aln;
	    $db->undefPointerCache();
	}
    }

    my $aln_mgr = DoTS::Gene::GenomeAlignmentMerger->new($db, $srt_alns, $self->{mc});
    $aln_mgr->getCompositeGenomeFeatures();
}

sub _breakIntoBins {
    my ($start, $end, $bin_size) = @_;

    my @res;
    for (my $s=$start; $s<$end; $s += $bin_size) {
	my $e = $s+$bin_size; $e = $end if $e > $end;
	push @res, [$s, $e];
    }
    return @res;
}

sub _getQueryReversed {
  my ($dbh, $genome_id, $blat_signals_cache, $aggressive) = @_;
  my $sql =<<SQL
SELECT blat_alignment_id FROM $blat_signals_cache
WHERE target_external_db_release_id = $genome_id
AND ((splice_signal_score_opposite > splice_signal_score and
      splice_signal_score_opposite > 1) or
     (splice_signal_score_opposite > splice_signal_score and
      polya_signal_score_opposite > polya_signal_score))
SQL
;
  if ($aggressive) {
    $sql =<<AGSQL
SELECT blat_alignment_id FROM $blat_signals_cache
WHERE target_external_db_release_id = $genome_id
AND (splice_signal_score_opposite > splice_signal_score or
     (splice_signal_score_opposite <= splice_signal_score and
      polya_signal_score_opposite > polya_signal_score));
AGSQL
;
  }
  my $sth = $dbh->prepareAndExecute($sql);

  my %res;
  while (my ($bid) = $sth->fetchrow_array) { $res{$bid} = 1; }

  \%res;
}

1;
