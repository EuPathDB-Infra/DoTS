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

  $self->logCommit;
  $self->logArgs;

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

  $self->logData("seeding clusters ...");
  my $clusters = &getClusterSeeds($dbh, $aid, $genome_id, $taxon_id, $test_chr, $query_dbid);
  $self->logData("number of clusters seeded: " . scalar(keys %$clusters) . '.');

  $self->logData("adding new members to clusters ...");
  my ($new_seqs, $changed_dgs, $no_dg_seqs) =
      &addNewMembers($dbh, $aid, $taxon_id, $query_dbid, $genome_id, $stage, $clusters);
  $self->logData("$new_seqs new seqs, enriched $changed_dgs existing clusters, $no_dg_seqs do not overlap existing clusters");

  $self->logData("get links between clusters ...");
  my $clnks = &getClusterLinks($clusters);
  $self->logData("number of links found: " . scalar(keys %$clnks) . ".");

  $self->logData("get linked cluster groups ...");
  my $cgrps = &getLinkedClusterGroups($clnks);
  $self->logData("number of linked cluster groups: " . scalar(keys %$cgrps) . ".");

  $self->logData("merging linked clusters ...");
  &mergeLinkedClusters($clusters, $cgrps);
  $self->logData("number of final clusters: " . scalar(keys %$clusters) . ".");

  $self->logData("writing into $out_file ...");
  open O, ">$out_file" or die "could not write $out_file";

  my $cKeys = &getClusterKeyList($clusters, $sort);

  my $c = 0;
  my $biggest = 0;
  foreach (@$cKeys) {
      my @seqs = keys %{ $clusters->{$_} };
      my $csize = scalar(@seqs);
      print O "Cluster_" . (++$c) . " ($csize sequences): (" . join(', ', @seqs) . ")\n";
      $biggest = $csize if $csize > $biggest;
  }
  close O;
  $self->logData("biggest cluster has $biggest sequences.");
}

####################

sub getClusterKeyList {
    my ($clusters, $sort) = @_;

    if (!$sort) { my @res = keys %$clusters; return \@res; }

    my %countHash;
    foreach my $ck (keys %$clusters) {
	my $c = $clusters->{$ck};
	my $sz = scalar(keys %$c);
	$countHash{$sz} = [] unless $countHash{$sz};
	push @{ $countHash{$sz} }, $ck;
    }

    my $res = [];
    foreach (sort {$a <=> $b} keys %countHash) {
	push @$res, @{ $countHash{$_} };
    }

    $res;
}

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
    # print "****  TotLnks: " . scalar(keys %$clnks) . "\n";

    # delete member seqs that only "link" one cluster
    # when >1 member seqs link the same set of clusters, keep any link
    my %unique_clnks;
    foreach my $lid (keys %$clnks) {
	my $l = $clnks->{$lid};
	if (scalar(keys %$l) == 1) {
	    delete $clnks->{$lid};
	} else {
	    my $uk = join(',', sort keys %$l);
	    $unique_clnks{$uk} = $lid;
	}
    }
    # print "****  Lnks w >1 DGs: " . scalar(keys %$clnks) . "\n";
    # print "****  UniqLnks: " . scalar(keys %unique_clnks) . "\n";

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

    return if ! defined $aid;

    my $new_seqs = 0;
    my %changed_dgs;
    my %assSeqs;
    my $no_dg_seqs = 0;

    my $sql = "select s.assembly_sequence_id, b.target_start, b.target_end, b.target_na_sequence_id "
	. "from DoTS.BlatAlignment b, DoTS.AssemblySequence s "
	. "where b.query_na_sequence_id = s.na_sequence_id "
	. "and b.query_table_id = 57 and b.query_taxon_id = $taxon_id "
	. "and b.query_external_db_release_id = $query_dbid "
	. "and b.target_table_id = 245 and b.target_taxon_id = $taxon_id "
	. "and b.target_external_db_release_id = $genome_id and b.is_best_alignment = 1";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";

    while (my ($id, $ts, $te, $chr_id) = $sth->fetchrow_array) {
      $assSeqs{$id} = [$ts,$te,$chr_id];
    }
    $sth->finish();

    my $sql1 = "select ag.aligned_gene_id from Allgenes.AlignedGene ag, Dots.VirtualSequence vs "
      . "where ag.chromosome = vs.chromosome and vs.na_sequence_id = ? "
	. "and not (ag.chromosome_start > ? or ag.chromosome_end < ?) "
	  . "and ag.aligned_gene_analysis_id = $aid";
    my $sth1 = $dbh->prepare($sql1) or die "bad sql $sql1";
    foreach my $assSeq (keys  %assSeqs) {
      $sth1->execute($assSeqs{$assSeq}[2],$assSeqs{$assSeq}[1],$assSeqs{$assSeq}[0]) or die "could not run $sql1";
      $new_seqs++;
      my $has_dg_overlap = 0;
      while (my ($agid) = $sth1->fetchrow_array) {
	# TODO: maybe evaluation exon overlaps by comparing block coordinates
	$clusters->{$agid}->{$assSeq} = '';
	$changed_dgs{$agid} = '';
	$has_dg_overlap++;
      }
      $sth1->finish;

      # HACK: put seqs on each chr in a group if they do not overlap any dgs
      unless ($has_dg_overlap) {
	$no_dg_seqs++;
	$clusters->{"c$assSeqs{$assSeq}[2]"} = {} unless ($clusters->{"c$assSeqs{$assSeq}[2]"});
	$clusters->{"c$assSeqs{$assSeq}[2]"}->{$assSeq} = '';
	# print "*** adding seq $assSeqs on chr $chr_id to a group\n";
      }
    }

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
    # die "could not get latest aligned_gene_analysis_id for taxon $taxon_id" unless defined($aid);
    $aid;
}

sub getClusterSeeds {
    my ($dbh, $aid, $genome_id, $taxon_id, $testChr, $query_dbid) = @_;

    return &makeClusters($dbh, $genome_id, $taxon_id, $testChr, $query_dbid) if !defined $aid;

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

sub makeClusters {
    my ($dbh, $genome_id, $taxon_id, $testChr, $query_dbid) = @_;

    my @seqs = &getSequencePieces($dbh, $genome_id, $testChr);

    my $clusters = {};
    my $cinfo = {id=>0, start=>0, end=>0};
    foreach my $seqId (@seqs) {
        my $sql = "select '' || s.assembly_sequence_id as sid, b.target_start, b.target_end "
        . "from DoTS.BlatAlignment b, DoTS.AssemblySequence s "
        . "where b.query_na_sequence_id = s.na_sequence_id "
        . "and b.query_table_id = 57 and b.query_taxon_id = $taxon_id "
        . "and b.query_external_db_release_id = $query_dbid "
        . "and b.target_table_id = 245 and b.target_taxon_id = $taxon_id "
        . "and b.target_na_sequence_id = $seqId and b.is_best_alignment = 1 "
        . "union "
        . "select 'DT.' || a.na_sequence_id as sid, b.target_start, b.target_end "
        . "from Dots.BlatAlignment b, Dots.Assembly a "
        . "where b.query_na_sequence_id = a.na_sequence_id "
        . "and b.query_table_id = 56 and b.query_taxon_id = $taxon_id "
        . "and b.target_table_id = 245 and b.target_taxon_id = $taxon_id "
        . "and b.target_na_sequence_id = $seqId and b.is_best_alignment = 1 ";
        $sql = "select * from ($sql) order by target_start asc, target_end asc";

        my $sth = $dbh->prepare($sql) or die "bad sql $sql";
        $sth->execute or die "could not run $sql";

        while (my ($sid, $s, $e) = $sth->fetchrow_array) {
            # make this a new command line arg
            my $CDIST = 100000;
            if ($cinfo->{end} == 0 || $s - $cinfo->{start} > $CDIST) {
                my $newCid = $cinfo->{id} + 1;
                $cinfo = {id=>$newCid, start=>$s, end=>$e};
                $clusters->{$newCid}->{$sid} = '';
            } else {
                my $cid = $cinfo->{id};
                $clusters->{$cid}->{$sid} = '';
                $cinfo->{end} = $e if $e > $cinfo->{end};
            }
        }
        $sth->finish;
    }

    return $clusters;
}

sub getSequencePieces {
    my ($dbh, $genome_id, $testChr) = @_;

    my $sql = "select na_sequence_id from dots.VirtualSequence where external_database_release_id = $genome_id";
    $sql .= " where chromosome = '$testChr'" if $testChr;
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    my @seqs = ();
    while (my ($id) = $sth->fetchrow_array) {
        push @seqs, $id;
    }
    @seqs;
}

1;

