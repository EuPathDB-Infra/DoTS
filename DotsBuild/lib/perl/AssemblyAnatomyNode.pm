package DoTS::DotsBuild::AssemblyAnatomyNode;


sub new {
  my ($class, $anatomyId, $parentNode) = @_;

  my $self = {};

  bless($self, $class);

  $self->{anatomyId} = $anatomyId;
  $self->{parentNode} = $parentNode;
  $self->{kids} = [];
  $self->{estCount} = 0;
  return $self;
}

sub addKid {
  my ($self, $kid) = @_;

  push(@{$self->{kids}}, $kid);
}

# zero out counts pertinent to a dt
sub clearDTValues {
  my ($self) = @_;

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
  $self->{dtRawPercolated} = $self->{dtRaw};
  $self->{dtEffectivePercolated} = $self->getDTEffective();

  foreach my $kid (@{$self->{kids}}) {
    ($kid_raw_percolated, $kid_effective_percolated)
      = $kid->percolateAndWrite($dtId, $sumEffective, $sumRaw, $taxonId);

    $self->{dtRawPercolated} += $kid_raw_percolated;
    $self->{dtEffectivePercolated} += $kid_effective_percolated;
  }
  return ($self->{dtRawPercolated},$self->{dtEffectivePercolated}) if $self->{dtRawPercolated} == 0;
  return if $self->{anatomyId} == 0;
  # write
  my $percent = $self->{dtEffectivePercolated} / $sumEffective * 100;
  my $anatomyESTs = $self->{dtRawPercolated};
  
  my $args = {percent => $percent,
	      anatomy_ests => $anatomyESTs,
	      anatomy_id => $self->{anatomyId},
	      EST_count => $sumRaw,
	      taxon_id => $taxonId,
	      na_sequence_id => $dtId
	     };
  my $assemblyAnatPercent = new GUS::Model::DoTS::AssemblyAnatomyPercent($args);
  $assemblyAnatPercent->submit();
  return ($self->{dtRawPercolated},$self->{dtEffectivePercolated});
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
  my ($self) = @_;
  my $numEST;
  if ($self->{ESTCount} == 0) {
    $numEST = 0;
  }
  else {
    $numEST = $self->{dtRaw}/$self->{ESTCount};
  }
  return $numEST;
}

sub printNode {
  my ($self,$indentation) = @_;
  print  "$indentation $self->{anatomyId}\t$self->{ESTCount}\t$self->{dtRaw}\t$self->{dtRawPercolated}\t$self->{dtEffectivePercolated}\n" if $self->{dtRawPercolated} > 0;
  $indentation = "$indentation  " ;
  foreach my $kid (@{$self->{kids}}) {
    $kid->printNode($indentation);
  }
}


1;
