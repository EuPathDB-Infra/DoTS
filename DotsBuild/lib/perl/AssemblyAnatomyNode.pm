package DoTS::DotsBuild::AssemblyAnatomyNode;

sub new {
  my ($class, $anatomyId, $parentNode) = @_;

  my $self = {};

  bless($self, $class);

  $self->{anatomyId} = $anatomyId;
  $self->{parentNode} = $parentNode;
  $self->{kids} = [];
  $self->{estCount} = 0;
}

sub addKid {
  my ($self, $kid) = @_;

  push(@{$self->{kids}}, $kid);
}

# zero out counts pertinent to a dt
sub clearDTValues {
  my ($self) = @_;

  return if $self->{dtRawPercolated} == 0;

  $self->{dtRaw} = 0;
  $self->{dtRawPercolated} = 0;
  $self->{dtEffectivePercolated} = 0;

  foreach my $kid (@{$self->{kids}}) {
    $kid->clearDTValues();
  }
}

# bottom up through the tree, percolate raw and effective. write each node
sub percolateAndWrite {
  my ($self, $dtId, $sumEffective, $sumRaw, $taxonId) = @_;

  my $raw_percolated = $self->{dtRaw};
  my $effective_percolated = $self->getDTEffective();

  foreach my $kid (@{$self->{kids}}) {
    ($kid_raw_percolated, $kid_effective_percolated)
      = $kid->percolateAndWrite();

    $raw_percolated += $kid_raw_percolated;
    $effective_percolated += $kid_effective_percolated;
  }

  # write
  my $percent = $effective_percolated / $sumEffective;
  my $anatomyESTs = $raw_percolated;

  my $args = {percent => $effective_percolated / $sumEffective * 100,
	      anatomy_ests => $raw_percolated,
	      anatomy_id => $self->{anatomyId},
	      EST_count => $sumRaw,
	      taxon_id => $taxonId,
	      na_sequence_id => $dtId
	     }
  my $assemblyAnatPercent = new GUS::Model::DoTS::AssemblyAnatomyPercent($args);
  $assemblyAnatPercent->submit();
}


sub setESTCount {
  my ($self, $count) = @_;

  $self->{ESTCount} = $count;
}


sub setDTRaw {
  my ($self, $count) = @_;

  $self->{dtRaw} = $count;
}

sub getDTEffective {
  return $self->{dtRaw}/$self->{ESTCount};

}
