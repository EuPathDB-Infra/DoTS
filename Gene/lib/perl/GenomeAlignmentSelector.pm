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
    my($class, $db, $qSel, $tSel, $optQSel, $optTSel) = @_;

    croak(__PACKAGE__ . "::new: query handle required") unless $db;

    my $dbh = $db->getQueryHandle();
    my $self = { db=>$db, dbh=>$dbh, qsel=>$qSel, tsel=>$tSel, oqsel=>$optQSel, otsel=>$optTSel };

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

    my $dbh = $self->{dbh};
    my $cntSql = $self->getSelectSql($sort, 1);
    print STDERR "running cntSql: $cntSql ...\n";
    my $cntSth = $dbh->prepareAndExecute($cntSql);
    my $cnt = $cntSth->fetchrow_array();
    if ($cnt > 4000) {
	print STDERR "increasing object cache size to 2 * $cnt + 2000 (>10000)\n";
	$self->{db}->setMaximumNumberOfObjects(2 * $cnt + 2000);
    }

    my $sql = $self->getSelectSql($sort);
    print STDERR "running sql: $sql ...\n";
    my $sth = $dbh->prepareAndExecute($sql);
    $res = [];
    while(my $id = $sth->fetchrow_array) {
	push @$res, GUS::Model::DoTS::BLATAlignment->new({blat_alignment_id=>$id});
    }
    print STDERR "tot alignments selected: " . scalar(@$res) . "\n";

    $self->{alignments} = $res;
    return $self->{alignments};
}

sub getSelectSql{
    my ($self, $sort, $countOnly) = @_;

    my $qSel = $self->{qsel};
    my $tSel = $self->{tsel};
    my $oqSel = $self->{oqsel};
    my $otSel = $self->{otsel};

    my $sql = "select " . ($countOnly ? 'count(*)' : 'blat_alignment_id') . " from DoTS.BlatAlignment";
    my @sqlWs;
    foreach my $k (keys %$qSel) { push @sqlWs, "$k = $qSel->{$k}"; }
    foreach my $k (keys %$tSel) { push @sqlWs, "$k = $tSel->{$k}"; }

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

    push @sqlWs, "target_start >= $otSel->{start}" if $otSel->{start};
    push @sqlWs, "target_start < $otSel->{end}" if $otSel->{end};

    if (scalar(@sqlWs) > 0) {
	$sql .= " where " . join(' and ',  @sqlWs);
    }

    $sql .= " order by target_start asc, target_end asc" if $sort;

    return $sql;
}

1;
