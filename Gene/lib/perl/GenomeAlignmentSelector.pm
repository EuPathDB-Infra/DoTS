#!/usr/bin/perl

# ------------------------------------------------------------------
# select a set of genome alignments satisfying given set of criteria
#
# created July 26, 2002 by ygan@pcbi.upenn.edu
#
# $Version$ $Date$ $Author$
# ------------------------------------------------------------------

package DoTS::Gene::GenomeAlignmentSelector;

use strict;
use GUS::Model::DoTS::BLATAlignment;
require Exporter;

# constructor
#
sub new {
    my($class, $dbh, $qSel, $tSel, $optQSel) = @_;

    croak(__PACKAGE__ . "::new: query handle required") unless $dbh;

    my $self = { dbh => $dbh, qsel => $qSel, tsel => $tSel, oqsel => $optQSel };

    bless $self, $class;
    return $self;
}

########## subroutine ###################################################

# return selected genomic alignment in this genomic region,
# sorted by start and end locations on the genome
# return type: array ref of GUS::Model::DoTS::BLATAlignment
#
sub getSortedAlignments { $_[0]->getAlignments(1); }

sub getAlignments {
    my ($self, $sort) = @_;

    my $res = $self->{alignments};
    return $res if defined $res;

    my $sql = $self->getSelectSql($sort);
    my $dbh = $self->{dbh};
    my $sth = $dbh->prepareAndExecute($sql);

    $res = [];
    while(my $id = $sth->fetchrow_array) {
	my $ba = GUS::Model::DoTS::BLATAlignment->new($id);
	my $basc = $self->{oqsel}->{blat_signals_cache};
	if ($basc) {
	    my $isRev = $ba->getIsReversed();
	    if (&isQueryReversed($dbh, $id, $basc)) {
		$ba->setIsReversed(1-$isRev);
	    }
	}
	push @$res, $ba;
    }

    $self->{alignments} = $res;
    return $self->{alignments};
}

sub getSelectSql{
    my ($self, $sort) = @_;

    my $qSel = $self->{qsel};
    my $tSel = $self->{tsel};
    my $oqSel = $self->{oqsel};

    my $sql = "select blat_alignment_id from DoTS.BlatAlignment";
    my @sqlWs;
    foreach (keys %$qSel) { push @sqlWs, "$_ = $qSel->{$_}"; }
    foreach (keys %$tSel) { push @sqlWs, "$_ = $tSel->{$_}"; }
    my $sp = $oqSel->{splice_minimum};
    my $mrna = $oqSel->{contains_mrna};
    my $xs = $oqSel->{exclude_singleton};
    if ($sp) { push @sqlWs, "max_target_gap >= $sp";  }
    if ($mrna || $xs) {
	my $sub = "select na_sequence_id from DoTS.Assembly where";
	if ($mrna && $xs) {
	    $sub .= " contains_mrna = 1 and number_of_contained_sequences > 1";
	} elsif ($xs) {
	    $sub .= "number_of_contained_sequences > 1";
	} else { 
	    $sub .= " contains_mrna = 1";
	}
	push @sqlWs, "query_na_sequence_id in ($sub)";
    }

    if (scalar(@sqlWs) > 0) {
	$sql .= " where " . join(' and ',  @sqlWs);
    }

    $sql .= " order by target_start asc, target_end asc" if $sort;

    return $sql;
}

sub isQueryReversed {
  my ($dbh, $bid, $blat_signals_cache, $aggressive) = @_;
  my $sql =<<SQL
SELECT count(*) FROM $blat_signals_cache
WHERE blat_alignment_id = $bid
AND ((splice_signal_score_opposite > splice_signal_score and
      splice_signal_score_opposite > 1) or
     (splice_signal_score_opposite > splice_signal_score and
      polya_signal_score_opposite > polya_signal_score))
SQL
;
  if ($aggressive) {
    $sql =<<AGSQL
SELECT count(*) FROM $blat_signals_cache
WHERE blat_alignment_id = $bid
AND (splice_signal_score_opposite > splice_signal_score or
     (splice_signal_score_opposite <= splice_signal_score and
      polya_signal_score_opposite > polya_signal_score));
AGSQL
;
  }
  my $sth = $dbh->prepareAndExecute($sql);
  my ($c) = $sth->fetchrow_array;
  $c;
}

1;
