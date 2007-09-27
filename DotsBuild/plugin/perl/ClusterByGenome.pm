package DoTS::DotsBuild::Plugin::ClusterByGenome;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::PluginMgr::Plugin;
use GUS::ObjRelP::DbiDatabase;

my $purposeBrief = <<PURPOSEBRIEF;
DoTS Clustering using genome alignments
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
DoTS Clustering using genome alignments
PLUGIN_PURPOSE

#check the documentation for this
my $tablesAffected = [];

my $tablesDependedOn = [];

my $howToRestart = <<PLUGIN_RESTART;
PLUGIN_RESTART

my $failureCases = <<PLUGIN_FAILURE_CASES;
PLUGIN_FAILURE_CASES

my $notes = <<PLUGIN_NOTES;
PLUGIN_NOTES

my $documentation = {
             purposeBrief => $purposeBrief,
		     purpose => $purpose,
		     tablesAffected => $tablesAffected,
		     tablesDependedOn => $tablesDependedOn,
		     howToRestart => $howToRestart,
		     failureCases => $failureCases,
		     notes => $notes
		    };


my $argsDeclaration =
      [
      stringArg({
          name => 'stage',
          descr => 'stage of clustering in DoTS build pipeline',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      integerArg({
          name => 'taxon_id',
          descr => 'taxon id',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
      integerArg({
          name => 'query_db_rel_id',
          descr => 'query database release id',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      integerArg({
          name => 'target_db_rel_id',
          descr => 'target external database release id',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
      stringArg({
          name => 'target_table_name',
          descr => 'name of table containing target sequences',
          default => 'VirtualSequence',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
          name => 'out',
          descr => 'output file for clustering result',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
          name => 'test_chr',
          descr => 'chromosome for test',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      booleanArg({
          name => 'mixedESTs',
          descr => 'query seqs are ESTs from multiple sources, use alternative sql to make clusters',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }), 
      booleanArg({
          name => 'sort',
          descr => 'whether to sort the output by cluster size (ascending)',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      integerArg({
          name => 'distanceBetweenStarts',
          descr => 'alignments are included in the same clusters if their starts are less than this number, default is 100,000',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      })
];

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);

  $self->initialize({requiredDbVersion => 3.5,
		     cvsRevision => '$Revision$', # cvs fills this in!
		     name => ref($self),
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation
		     });

  return $self;
}

$| = 1;

sub run {
  my $self   = shift;

  $self->logCommit;
  $self->logArgs;

  my $dbh = $self->getQueryHandle();

  my $stage = $self->getArg('stage') if $self->getArg('stage');
  my $taxon_id = $self->getArg('taxon_id');
  my $genome_id = $self->getArg('target_db_rel_id');
  my $query_dbid = $self->getArg('query_db_rel_id');
  my $test_chr = $self->getArg('test_chr');
  my $out_file = $self->getArg('out');
  my $sort = $self->getArg('sort');
  my $target_table_name = $self->getArg('target_table_name');
  my $target_table_id = $self->getTableId($dbh, $target_table_name);
  
  my $aid = $self->getAlignedGeneAnalysisId($dbh,$taxon_id, $genome_id);

  $self->logData("seeding clusters ...");
  my $clusters = $self->getClusterSeeds($dbh, $aid, $genome_id, $target_table_name, $taxon_id, $test_chr, $query_dbid);
  $self->logData("number of clusters seeded: " . scalar(keys %$clusters) . '.');

  $self->logData("adding new members to clusters ...");
  my ($new_seqs, $changed_dgs, $no_dg_seqs) =
      $self->addNewMembers($dbh, $aid, $taxon_id, $target_table_id, $query_dbid, $genome_id, $stage, $clusters);
  $self->logData("$new_seqs new seqs, enriched $changed_dgs existing clusters, $no_dg_seqs do not overlap existing clusters");

  $self->logData("get links between clusters ...");
  my $clnks = $self->getClusterLinks($clusters);
  $self->logData("number of links found: " . scalar(keys %$clnks) . ".");

  $self->logData("get linked cluster groups ...");
  my $cgrps = $self->getLinkedClusterGroups($clnks);
  $self->logData("number of linked cluster groups: " . scalar(keys %$cgrps) . ".");

  $self->logData("merging linked clusters ...");
  $self->mergeLinkedClusters($clusters, $cgrps);
  $self->logData("number of final clusters: " . scalar(keys %$clusters) . ".");

  $self->logData("writing into $out_file ...");
  open O, ">$out_file" or die "could not write $out_file";

  my $cKeys = $self->getClusterKeyList($clusters, $sort);

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
    my ($self,$clusters, $sort) = @_;

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
    my ($self,$clusters, $cgrps) = @_;

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
    my ($self,$clnks) = @_;

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
    my ($self,$clusters) = @_;

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
    my ($self,$dbh, $aid, $taxon_id, $target_table_id, $query_dbid, $genome_id, $stage, $clusters) = @_;

    return if ! defined $aid;

    my $new_seqs = 0;
    my %changed_dgs;
    my %assSeqs;
    my $no_dg_seqs = 0;

    my $sql = <<"EOSQL";
    SELECT s.assembly_sequence_id, 
           b.target_start,
           b.target_end, 
           b.target_na_sequence_id 
    FROM   DoTS.BlatAlignment b, 
           DoTS.AssemblySequence s,
           core.tableinfo ti
    WHERE  b.query_na_sequence_id = s.na_sequence_id 
      AND  b.query_table_id = ti.table_id  
      AND  lower(ti.name) = 'assemblysequence' 
      AND  b.query_taxon_id = $taxon_id 
      AND  b.query_external_db_release_id = $query_dbid 
      AND  b.target_table_id = $target_table_id 
      AND  b.target_taxon_id = $taxon_id 
      AND  b.target_external_db_release_id = $genome_id 
      AND  b.is_best_alignment = 1
EOSQL

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
    my ($self,$dbh,$taxon_id, $genome_id) = @_;

    GUS::PluginMgr::Plugin->logData("no AllGenes schema, skipping aligned_gene_analysis_id lookup.") and return unless $self->haveAllgenesSchema($dbh);

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
    my ($self, $dbh, $aid, $genome_id, $target_table_name, $taxon_id, $testChr, $query_dbid) = @_;

    return $self->makeClusters($dbh, $genome_id, $target_table_name, $taxon_id, $testChr, $query_dbid) if !defined $aid;

    my $sql = "select ag.aligned_gene_id, 'DT.' || aga.na_sequence_id "
	. "from Allgenes.AlignedGene ag, Allgenes.AlignedGeneAssembly aga "
	. "where ag.aligned_gene_analysis_id = $aid and ag.aligned_gene_id = aga.aligned_gene_id";
    $sql .= " and ag.chromosome = '$testChr'" if $testChr;
    
    print STDERR "getClsterSeeds sql: $sql\n\n"; # DEBUG
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
    my ($self, $dbh, $genome_id, $target_table_name, $taxon_id, $testChr, $query_dbid) = @_;

    my @seqs = $self->getSequencePieces($dbh, $target_table_name, $genome_id, $testChr); 

    my $targetNumber = @seqs;

    print STDERR "Total number of target pieces = $targetNumber\n";

    my $target_table_id = $self->getTableId($dbh, $target_table_name);

    my $clusters = {};
    my $cinfo = {id=>0, start=>0, end=>0};
    foreach my $seqId (@seqs) {
        my $sql = <<"EOSQL";
        SELECT '' || s.assembly_sequence_id AS sid, 
               b.target_start, b.target_end 
        FROM   DoTS.BlatAlignment b, 
               DoTS.AssemblySequence s, 
               core.tableinfo ti 
        WHERE  b.query_na_sequence_id = s.na_sequence_id 
          AND  b.query_table_id = ti.table_id  
          AND  lower(ti.name) = 'assemblysequence' 
          AND  b.query_taxon_id = $taxon_id 
          AND  b.query_external_db_release_id = $query_dbid 
          AND  b.target_table_id = $target_table_id  
          AND b.target_taxon_id = $taxon_id 
          AND  b.target_na_sequence_id = $seqId  
          AND b.is_best_alignment = 1 
        UNION
        SELECT 'DT.' || a.na_sequence_id AS sid, 
               b.target_start, b.target_end 
        FROM   Dots.BlatAlignment b, 
               Dots.Assembly a,
               core.tableinfo ti 
        WHERE  b.query_na_sequence_id = a.na_sequence_id 
          AND  b.query_table_id = ti.table_id  
          AND  lower(ti.name) = 'assembly'  
          AND  b.query_taxon_id = $taxon_id 
          AND  b.target_table_id = $target_table_id  
          AND  b.target_taxon_id = $taxon_id 
          AND  b.target_na_sequence_id = $seqId
          AND  b.is_best_alignment = 1
EOSQL

	my $altSql = <<"EOSQL";
	SELECT s.assembly_sequence_id AS sid, b.target_start, b.target_end 
        FROM DoTS.BlatAlignment b,
             DoTS.AssemblySequence s,
             core.tableinfo ti,
             dots.externalnasequence x,
             sres.sequenceontology so  
        WHERE  b.query_na_sequence_id = s.na_sequence_id 
          AND  b.query_table_id = ti.table_id
          AND  lower(ti.name) = 'assemblysequence'
          AND  b.query_taxon_id = $taxon_id 
          AND  b.target_table_id = $target_table_id
          AND  b.target_na_sequence_id = $seqId
          AND  b.is_best_alignment = 1
	  AND  s.na_sequence_id = x.na_sequence_id
          AND  x.sequence_ontology_id = so.sequence_ontology_id
          AND  so.term_name = 'EST'
EOSQL


	$sql = $altSql if $self->getArg('mixedESTs');

        $sql = "select * from ($sql) order by target_start asc, target_end asc";

	print STDERR "SQL to make clusters: $sql\n";

        my $sth = $dbh->prepare($sql) or die "bad sql $sql";
        $sth->execute or die "could not run $sql";

        while (my ($sid, $s, $e) = $sth->fetchrow_array) {
            # make this a new command line arg
            my $CDIST = $self->getArg('distanceBetweenStarts') ? $self->getArg('distanceBetweenStarts') : 100000;
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


    my $clusterNum = scalar (keys %{$clusters});

    print STDERR "Number of clusters from make clusters: $clusterNum\n";
    return $clusters;
}

sub getSequencePieces {
    my ($self,$dbh, $target_table_name, $genome_id, $testChr) = @_;

    my $sql = "select na_sequence_id from dots.${target_table_name} where external_database_release_id = $genome_id";
    $sql .= " where chromosome = '$testChr'" if $testChr;
   
    my $sth = $dbh->prepare($sql) or die "bad sql $sql";
    $sth->execute or die "could not run $sql";
    my @seqs = ();
    while (my ($id) = $sth->fetchrow_array) {
        push @seqs, $id;
    }
    @seqs;
}

sub getTableId {
    my ($self,$dbh, $target_table_name) = @_;
    my $sth = $dbh->prepare("select table_id from core.tableinfo where lower(name) = lower('$target_table_name')");
    $sth->execute(); 
    my ($id) = $sth->fetchrow();
    $sth->finish();     
    return $id;
}

sub haveAllgenesSchema {
    my ($self,$dbh) = @_;
    my $sth = $dbh->prepare("select count(*) from all_tables where owner = 'ALLGENES'");
    $sth->execute();
    my ($ct) = $sth->fetchrow();
    $sth->finish();
    return $ct;
}

1;

