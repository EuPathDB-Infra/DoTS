package DoTS::Gene::ESTOrientationRecord;

# -----------------------------------------------------------
# data container for ESTOrientationPlot.pm
#
# Yongchang Gan, Jan 24, 03
# ---------------------------------------------------------

use strict;
use DoTS::Gene::ESTOrientationEntry;

my $CN = 'ESTOrientationRecord';

# constructor
#
sub new {
    my($class, $args) = @_;

    my $TAG = $CN . "::new:";

    my $type = lc $args->{type};
    die "$TAG Type must of DT or DG" unless $type eq 'dt' or $type eq 'dg';

    my $self = { type => $type,
		 id => $args->{id},
		 debug => $args->{debug}
	     };

    # other data members
    $self->{dbh} = $args->{dbh};
    $self->{gdg_tab} = $args->{gdg_tab};
    $self->{gdt_tab} = $args->{gdt_tab};
    $self->{strand} = undef;
    $self->{chr_s} = undef;
    $self->{chr_e} = undef;
    $self->{entries} = undef;

    bless $self, $class;
    return $self;
}

sub getType {
    my $self = shift;
    return $self->{type};
}

sub getId {
    my $self = shift;
    return $self->{id};
}

sub getStrand {
    my $self = shift;
    return $self->{strand};
}

sub getChromosomeStart {
    my $self = shift;
    return $self->{chr_s};
}

sub getChromosomeEnd {
    my $self = shift;
    return $self->{chr_e};
}

# this is intended for aggregated processing
# to avoid querying db separately for every record
sub setRecord {
    my ($self, $info) = @_;

    $self->{type} = $info->{type} unless $self->{type};
    $self->{id} = $info->{id} unless $self->{id};

    if ($self->type eq 'dg') {
	$self->{strand} = $info->{strand};
	$self->{chr_s} = $info->{chr_s};
	$self->{chr_e} = $info->{chr_e};
	$self->{bs} = $info->{bs};
	$self->{qs} = $info->{qs};
	$self->{ts} = $info->{ts};
    }

    my @entries;
    my $ies = $info->{entries};
    foreach (@$ies) {
	push @entries, DoTS::Gene::ESTOrientationEntry->new($_);
    }
    $self->{entries} = \@entries;
}

sub getEntries {
    my $self = shift;

    my $TAG = $CN . "::getEntries:";
    my $dbg = $self->{debug};

    my $entries = $self->{entries};
    return $entries if $entries;

    # get entries from db
    my $type = $self->{type};
    my $id = $self->{id};
    my $dbh = $self->{dbh};
    my $sql = &_getESTOrientationQuery($type, $id, $self->{gdg_tab}, $self->{gdt_tab});
    my $sth = $dbh->prepare($sql) or die "bad sql $sql: $!\n";
    print "$TAG sql=\n$sql\n" if $dbg;
    $sth->execute;

    my @entries;
    while(my $h = $sth->fetchrow_hashref('NAME_lc')) {
	if ($dbg > 1) {
	    print $TAG . " ";
	    my @hks = sort keys %$h;
	    foreach (@hks) {
		print $_ . '=>' . $h->{$_} . '; ' unless $dbg < 3 && $_ eq 'gap_con';
	    }
	    print "\n";
	}
	if ($type eq 'dg') {
	    $self->{strand} = $h->{strand} unless $self->{strand};
	    $self->{chr_s} = $h->{chromosome_start} unless $self->{chr_s};
	    $self->{chr_e} = $h->{chromosome_end} unless $self->{chr_e};
	    $self->{bs} = $h->{blocksizes} unless $self->{bs};
	    $self->{qs} = $h->{qstarts} unless $self->{qs};
	    $self->{ts} = $h->{tstarts} unless $self->{ts};
	}
	push @entries, DoTS::Gene::ESTOrientationEntry->new($h);
    }

    $self->{entries} = \@entries;

    return $self->{entries};
}


#### file scoped subroutines

sub _getESTOrientationQuery {
    my ($type, $id, $gdg_tab, $gdt_tab) = @_;

    my $sql;
    if (lc $type eq 'dt') {
	$sql = <<DT_QUERY_END
SELECT est.p_end, length(s.sequence_gaps) as gap_seq_len, asm.length as asm_len,
       asm.gapped_consensus as gap_con, s.assembly_offset as asm_offset,
       s.na_sequence_id as seq_id
FROM DoTS.Assembly asm, DoTS.AssemblySequence s, DoTS.EST est
WHERE asm.na_sequence_id = $id
and asm.na_sequence_id = s.assembly_na_sequence_id
and s.na_sequence_id = est.na_sequence_id
DT_QUERY_END
    } else {
	$sql = <<DG_QUERY_END
SELECT g.strand, g.chromosome_start, g.chromosome_end,
       asm.na_sequence_id as asm_id, asm.length as asm_len, asm.gapped_consensus as gap_con,
       b.blat_alignment_id, b.blocksizes, b.qstarts, b.tstarts,
       s.assembly_offset as asm_offset,
       s.na_sequence_id as seq_id, est.p_end, length(s.gapped_sequence) as gap_seq_len
FROM $gdg_tab g, $gdt_tab a, DoTS.blatalignment b, DoTS.assembly asm,
     DoTS.assemblysequence s, DoTS.EST est, DoTS.virtualsequence v
WHERE g.genome_dots_gene_id = $id
and g.genome_dots_gene_id = a.genome_dots_gene_id
and a.na_sequence_id = b.query_na_sequence_id
and b.query_table_id = 56
and b.target_na_sequence_id = v.na_sequence_id
and v.chromosome = g.chromosome
and b.is_best_alignment = 1
and b.target_start >= g.chromosome_start and b.target_end <= g.chromosome_end
and asm.na_sequence_id = s.assembly_na_sequence_id
and a.na_sequence_id = s.assembly_na_sequence_id
and s.na_sequence_id = est.na_sequence_id
order by b.target_start, b.target_end
DG_QUERY_END
    }
}

1;
