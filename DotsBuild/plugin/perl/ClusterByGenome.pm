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
  my $out_file = $cla->{'out'};
  my $sort = $cla->{'sort'};

  my $aid = &getAlignedGeneAnalysisId($dbh,$taxon_id, $genome_id);
  my $clusters = &getClusterSeeds($dbh, $aid);

  &addNewMembers($dbh, $aid, $taxon_id, $query_dbid, $genome_id, $stage, $clusters);

  &mergeLinkedClusters($clusters);
}

####################
sub mergeLinkedClusters {
    my ($clusters) = @_;

    my $prev_size;
    do {
	my @cids = keys %$clusters;
	$prev_size = scalar(@cids);
	for (my $i=0; $i<$prev_size-1; $i++) {
	    my $cid1 = $cids[$i];
	    my $c1 = $clusters->{$cid1};
	    for (my $j=$i+1; $j<$prev_size; $j++) {
		my $cid2 = $cids[$j];
		my $c2 = $clusters->{$cid2});
		if (&isLinked($c1, $c2) {
		    foreach (keys %$c2) {
			$c1->{$_} = '';
		    }
		    delete $clusters->{$c2};
		}
	    }
	}
    } while (sclalar(keys %$clusters) < $prev_size);
}

sub isLinked {
    my ($c1, $c2) = @_;

    foreach (keys %$c1) {
	return 1 if exists $c2->{$_};
    }
    return 0;
}

sub addNewMembers {
    my ($dbh, $aid, $taxon_id, $query_dbid $genome_id, $stage, $clusters) = @_;

    my $sql = "select query_na_sequence_id, target_start, target_end, target_na_sequence_id "
	. "from DoTS.BlatAlignment where query_table_id = 89 and query_taxon_id = $taxon_id "
	. "and query_external_db_release_id = $query_dbid "
	. "and target_table_id = 245 and target_taxon_id = $taxon_id "
	. "and target_external_db_release_id = $genome_id and is_best_alignment = 1";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    while (my ($id, $ts, $te, $chr_id) = $sth->fetchrow_array) {
	my $sql1 = "select ag.aligned_gene_id from Allgenes.AlignedGene ag, Dots.VirtualSequence vs "
	    . "where ag.chr = vs.chromosome and vs.na_sequence_id = $chr_id "
	    . "and not (ag.chromosome_start > $te or ag.chromosome_end < $ts) "
	    . "and ag.aligned_gene_analysis_id = $aid";
	my $sth1 = $dbh->prepare($sql1) or die "bad sql $sql1";
	$sth1->execute or die "could not run $sql1";
	while (my ($agid) = $sth->fetchrow_array) {
	    # TODO: maybe evaluation exon overlaps by comparing block coordinates
	    $clusters->{$agid}->{$id} = '';
	}
	$sth1->finish;
    }
    $sth->finish;
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
    my ($dbh, $aid) = @_;

    my $sql = "select ag.aligned_gene_id, 'DT.' || aga.na_sequence_id "
	. "from Allgenes.AlignedGene ag, Allgenes.AlignedGeneAssembly aga "
	. "where ag.aligned_gene_analysis_id = $aid and ag.aligned_gene_id = aga.aligned_gene_id";
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    my $seed_clusters = {};
    while (my ($agid, $dt) = $sth->fetchrow_array) {
	$clusters->{$agid} = {} unless exists $clusters->{$agid};
	$clusters->{$agid}->{$dt} = '';
    }
    $sth->finish;
    $seed_clusters;
}


1;

