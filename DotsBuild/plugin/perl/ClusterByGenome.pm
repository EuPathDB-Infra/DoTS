package DoTS::DotsBuild::Plugin::ClusterByGenome;



@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::ObjRelP::DbiDatabase;

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);
  
  my $usage = 'DoTS Clustering using genome alignments';
  
  my $easycsp =
      [{o => 'stage',
	t => 'string',
	h => 'stage of clustering in DoTS build pipeline'
	},
       {o => 'taxon_id',
	t => 'int',
	h => 'taxon id'
	},
       {o => 'query_db_rel_id',
	t => 'int',
	h => 'query database release id'
	},
       {o => 'target_db_rel_id',
	t => 'int',
	h => 'target external database release id'
	},
       {o => 'out',
	t => 'string',
	h => 'output file for clustering result'
	},
       {o => 'test_chr',
	t => 'string',
	h => 'chromosome for test'
	},
       {o => 'sort',
	t => 'boolean',
	h => 'whether to sort the output by cluster size (ascending)'
	}
       ];
  
  $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$ $',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     revisionNotes => ' ',
		     easyCspOptions => $easycsp,
		     usage => $usage
		    });
  
  return $self;
}


$| = 1;

sub run {
  my $self   = shift;
  
  $self->log ($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");

  my $dbh = $self->getQueryHandle();
  my $cla = $self->getCla();

  my $stage = $cla->{'stage'};
  my $taxon_id = $cla->{'taxon_id'};
  my $genome_id = $cla->{'target_db_rel_id'};
  my $query_dbid = $cla->{'query_db_rel_id'};
  my $test_chr = $cla->{'test_chr'};
  my $out_file = $cla->{'out'};
  my $sort = $cla->{'sort'};

  my $aid = &getAlignedGeneAnalysisId($dbh,$taxon_id, $genome_id);

  print "seeding clusters ...\n";
  my $clusters = &getClusterSeeds($dbh, $aid, $test_chr);
  print "number of clusters seeded: " . scalar(keys %$clusters) . ".\n";

  print "adding new members to clusters ...\n";
  my ($new_seqs, $changed_dgs, $no_dg_seqs) =
      &addNewMembers($dbh, $aid, $taxon_id, $query_dbid, $genome_id, $stage, $clusters);
  print "$new_seqs new seqs, enriched $changed_dgs existing clusters, $no_dg_seqs do not overlap existing clusters\n";

  print "get links between clusters ...\n";
  my $clnks = &getClusterLinks($clusters);
  print "number of links found: " . scalar(keys %$clnks) . ".\n";

  print "get linked cluster groups ...\n";
  my $cgrps = &getLinkedClusterGroups($clnks);
  print "number of linked cluster groups: " . scalar(keys %$cgrps) . ".\n";

  print "merging linked clusters ...\n";
  &mergeLinkedClusters($clusters, $cgrps);
  print "number of final clusters: " . scalar(keys %$clusters) . ".\n";

  print "writing into $out_file ...\n";
  open O, ">$out_file" or die "could not write $out_file";
  my $c = 0;
  my $biggest = 0;
  foreach (keys %$clusters) {
      my @seqs = keys %{ $clusters->{$_} };
      my $csize = scalar(@seqs);
      print O ">Cluster_" . (++$c) . " ($csize sequences): (" . join(', ', @seqs) . ")\n";
      $biggest = $csize if $csize > $biggest;
  }
  close O;
  print "biggest cluster has $biggest sequences\n";
}

####################
sub mergeLinkedClusters {
    my ($clusters, $cgrps) = @_;

    foreach (keys %$cgrps) {
	my @sids = keys %{ $cgrps->{$_} };
	my $fst = $sids[0];
	my $count = scalar(@sids);
	for (my $i=1; $i<$count; $i++) {
	    my $nxt = $sids[$i];
	    my $nxt_c = $clusters->{$nxt};
	    foreach (keys %$nxt_c) {
		$clusters->{$fst}->{$_} = '';
	    }
	    delete $clusters->{$nxt};
	}
    }
}

sub getLinkedClusterGroups {
    my ($clnks) = @_;

    return $clnks unless $clnks;
    return $clnks unless scalar(keys %$clnks) > 1;

    my %cid_lid_map;
    foreach my $lid (keys %$clnks) {
	my $c1= $clnks->{$lid};
	my $g;
	foreach my $cid (keys %$c1) {
	    $g = $cid_lid_map{$cid} unless $g;
	}
	foreach my $cid (keys %$c1) {
	    $cid_lid_map{$cid} = ($g ? $g : $lid);
	}
    }

    $clnks = {};
    foreach my $cid (keys %cid_lid_map) {
	my $lid = $cid_lid_map{$cid};
	$clnks->{$lid} = {} unless $clnks->{$lid};
	$clnks->{$lid}->{$cid} = '';
    }

    $clnks;
}

sub isLinked {
    my ($h1, $h2) = @_;

    foreach (keys %$h1) {
	return 1 if exists $h2->{$_};
    }
    return 0;
}

sub getClusterLinks {
    my ($clusters) = @_;

   # clusters linked by a member seq
    my $clnks = {};
    foreach my $cid (keys %$clusters) {
	my $c = $clusters->{$cid};
	foreach my $mid (keys %$c) {
	    $clnks->{$mid} = {} unless $clnks->{$mid};
	    $clnks->{$mid}->{$cid} = '';
	}
    }

    # delete member seqs that only "link" one cluster
    # when >1 member seqs link the same set of clusters, keep any link
    my %unique_clnks;
    foreach my $lid (keys %$clnks) {
	my $l = $clnks->{$lid};
	delete $clnks->{$lid} unless scalar(keys %$l) > 1; 
	my $uk = join(',', sort keys %$l);
	$unique_clnks{$uk} = $lid;
    }

    $clnks = {};
    foreach my $k (keys %unique_clnks) {
	my @cids = split(/,/, $k);
	my $mid = $unique_clnks{$k};
	foreach my $cid (@cids) {
	    $clnks->{$mid} = {} unless $clnks->{$mid};
	    $clnks->{$mid}->{$cid} = '';
	}
    }

    $clnks;
}

sub addNewMembers {
    my ($dbh, $aid, $taxon_id, $query_dbid, $genome_id, $stage, $clusters) = @_;

    my $new_seqs = 0;
    my %changed_dgs;
    my $no_dg_seqs = 0;
    my $sql = "select s.assembly_sequence_id, b.target_start, b.target_end, b.target_na_sequence_id "
	. "from DoTS.BlatAlignment b, DoTS.AssemblySequence s"
	. "where b.query_na_sequence_id = s.na_sequence_id "
	. "b.query_table_id = 57 and b.query_taxon_id = $taxon_id "
	. "and b.query_external_db_release_id = $query_dbid "
	. "and b.target_table_id = 245 and b.target_taxon_id = $taxon_id "
	. "and b.target_external_db_release_id = $genome_id and b.is_best_alignment = 1";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    while (my ($id, $ts, $te, $chr_id) = $sth->fetchrow_array) {
	my $sql1 = "select ag.aligned_gene_id from Allgenes.AlignedGene ag, Dots.VirtualSequence vs "
	    . "where ag.chr = vs.chromosome and vs.na_sequence_id = $chr_id "
	    . "and not (ag.chromosome_start > $te or ag.chromosome_end < $ts) "
	    . "and ag.aligned_gene_analysis_id = $aid";
	my $sth1 = $dbh->prepare($sql1) or die "bad sql $sql1";
	$sth1->execute or die "could not run $sql1";
	$new_seqs++;
	my $has_dg_overlap = 0;
	while (my ($agid) = $sth->fetchrow_array) {
	    # TODO: maybe evaluation exon overlaps by comparing block coordinates
	    $clusters->{$agid}->{$id} = '';
	    $changed_dgs{$agid} = '';
	    $has_dg_overlap++;
	}
	$no_dg_seqs++ unless $has_dg_overlap;
	$sth1->finish;
    }
    $sth->finish;

    return ($new_seqs, scalar(keys %changed_dgs), $no_dg_seqs);
}

sub getAlignedGeneAnalysisId {
    my ($dbh,$taxon_id, $genome_id) = @_;
    my $sql = "select aligned_gene_analysis_id from Allgenes.AlignedGeneAnalysis "
	. "where target_external_db_id = $genome_id and parameters like '%--t $taxon_id %' "
	. "order by aligned_gene_analysis_id desc";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    my $aid;
    if (($aid) = $sth->fetchrow_array) { ; }
    $sth->finish;
    die "could not get latest aligned_gene_analysis_id for taxon $taxon_id" unless defined($aid);
    $aid;
}

sub getClusterSeeds {
    my ($dbh, $aid, $testChr) = @_;

    my $sql = "select ag.aligned_gene_id, 'DT.' || aga.na_sequence_id "
	. "from Allgenes.AlignedGene ag, Allgenes.AlignedGeneAssembly aga "
	. "where ag.aligned_gene_analysis_id = $aid and ag.aligned_gene_id = aga.aligned_gene_id";
    $sql .= " and ag.chromosome = '$testChr'" if $testChr;
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    my $seed_clusters = {};
    while (my ($agid, $dt) = $sth->fetchrow_array) {
	$seed_clusters->{$agid} = {} unless exists $seed_clusters->{$agid};
	$seed_clusters->{$agid}->{$dt} = '';
    }
    $sth->finish;
    $seed_clusters;
}


1;

