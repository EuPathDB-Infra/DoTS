package DoTS::Gene::ESTOrientEntry;

# -----------------------------------------------------------
# data container for ESTOrientGroup.pm
#
# Yongchang Gan, Jan 24, 03
# ---------------------------------------------------------

use strict;

# constructor
#
sub new {
    my($class, $args) = @_;

    my $self = {};

    $self->{p_end} = $args->{p_end};
    $self->{asm_offset} = $args->{asm_offset};
    $self->{gap_seq_len} = $args->{gap_seq_len};
    $self->{asm_id} = $args->{asm_id};
    $self->{asm_len} = $args->{asm_len};
    $self->{blat_id} = $args->{blat_alignment_id};
    $self->{bs} = $args->{blocksizes};
    $self->{qs} = $args->{qstarts};
    $self->{ts} = $args->{tstarts};
    $self->{gap_con} = $args->{gap_con};
    $self->{seq_id} = $args->{seq_id};

    bless $self, $class;
    return $self;
}

sub getPEnd {
    my $self = shift;
    return $self->{p_end};
}

sub getAssemblyOffset {
    my $self = shift;
    return $self->{asm_offset};
}

sub getGappedSeqLength {
    my $self = shift;
    return $self->{gap_seq_len};
}

sub getAssemblyId {
    my $self = shift;
    return $self->{asm_id};
}    

sub getAssemblyLength {
    my $self = shift;
    return $self->{asm_len};
}    

sub getBLATAlignmentId {
    my $self = shift;
    return $self->{blat_id};
}    

sub getBlocksizes {
    my $self = shift;
    my @bss = split(/,/, $self->{bs});
    return \@bss;
}

sub getQstarts {
    my $self = shift;
    my @qss = split(/,/, $self->{qs});
    return \@qss;
}

sub getTstarts {
    my $self = shift;
    my @tss = split(/,/, $self->{ts});
    return \@tss;
}

sub getGappedConsensus {
    my $self = shift;
    return $self->{gap_con};
}

sub getSequenceId {
    my $self = shift;
    return $self->{seq_id};
}

1;
