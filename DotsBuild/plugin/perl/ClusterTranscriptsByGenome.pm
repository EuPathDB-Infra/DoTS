package DoTS::DotsBuild::Plugin::ClusterTranscriptsByGenome;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::PluginMgr::Plugin;
use GUS::ObjRelP::DbiDatabase;

my $purposeBrief = <<PURPOSEBRIEF;
Clustering ESTs and mRNAs using genome alignments
PURPOSEBRIEF

my $purpose = <<PLUGIN_PURPOSE;
Transcript clustering using genome alignments
PLUGIN_PURPOSE

my $tablesAffected = [];

my $tablesDependedOn = [['DoTS::BlatAlignment', 'data for clustering derived from the blatalignment table'],['SRes::ExternalDatabaseRelease', 'must contain version for both query and target sequences'],['SRes::ExternalDatabase', 'must contain name for both query and target sequences'],['SRes::Taxon', 'must contain taxon information for ']];

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
      integerArg({
          name => 'ncbiTaxId',
          descr => 'tax_id from NCBI',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
      stringArg({
          name => 'queryDbName',
          descr => 'sres.externaldatabase.name for the query sequences in the blatalignment table',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
       stringArg({
          name => 'queryDbRlsVer',
          descr => 'sres.externaldatabaserelease.version for the query sequences in the blatalignment table',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
       stringArg({
          name => 'targetDbName',
          descr => 'sres.externaldatabase.name for the target sequences in the blatalignment table',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
       stringArg({
          name => 'targetDbRlsVer',
          descr => 'sres.externaldatabaserelease.version for the target sequences in the blatalignment table',
          constraintFunc => undef,
          reqd => 1,
          isList => 0
      }),
      stringArg({
          name => 'targetTableName',
          descr => 'name of table containing target sequences',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
          name => 'outFile',
          descr => 'output file for clustering result',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
	  name => 'gffDir',
          descr => 'output file for clustering result',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
          name => 'testChr',
          descr => 'chromosome for test',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      booleanArg({
          name => 'sortClusterSize',
          descr => 'whether to sort the output by cluster size (ascending)',
          constraintFunc => undef,
          reqd => 0,
          isList => 0
      }),
      stringArg({
	  name => 'ucscScoreTableSpace',
          descr => 'tablespace where ucscscore table is located',
          constraintFunc => undef,
          reqd => 1,
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

  $self->logData("making clusters");

  my $clusters = $self->makeClusters($dbh);

  my $outFile = $self->getArg('outFile');

  my $sort = $self->getArg('sortClusterSize');

  $self->logData("writing into $outFile ...");
  open (OUT, ">$outFile") || die "could not write $outFile\n";

  my $cKeys = &getClusterKeyList($clusters, $sort);

  my $c = 0;
  my $biggest = 0;
  foreach (@$cKeys) {
      my @seqs = keys %{ $clusters->{$_} };
      my $csize = scalar(@seqs);
      print OUT "Cluster_" . (++$c) . " ($csize sequences): (" . join(', ', @seqs) . ")\n";
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


sub makeClusters {
    my ($self, $dbh) = @_;

    my $taxonId = $self->getTaxonId($dbh);

    my $targetTableId = $self->getTableId($dbh);

    my $qDbRlsId = $self->getExtDbRlsId($self->getArg('queryDbName'),$self->getArg('queryDbRlsVer'));

    my $tDbRlsId = $self->getExtDbRlsId($self->getArg('targetDbName'),$self->getArg('targetDbRlsVer'));

    my $seqs = $self->getSequencePieces($dbh, $tDbRlsId);

    my $space = $self->getArg('ucscScoreTableSpace');

    my $clusters = {};

    my $cinfo = {id=>0, start=>0, end=>0};

    foreach my $seqId (keys %{$seqs}) {
       my $sql = <<"EOSQL";
         SELECT s.assembly_sequence_id, x.source_id
                b.target_start, b.target_end,
                b.tstarts, b.blocksizes
         FROM   DoTS.BlatAlignment b,
                DoTS.AssemblySequence s,
                core.tableinfo ti,
		$space.ucscscore u,
                Dots.ExternalNASequence x
         WHERE  b.query_na_sequence_id = s.na_sequence_id
           AND  x.na_sequence_id = s.na_sequence_id
           AND  b.query_table_id = ti.table_id
           AND  lower(ti.name) = 'assemblysequence'
	   AND  b.blat_alignment_id = u.blat_align_id
           AND  b.query_taxon_id = ?
           AND  b.query_external_db_release_id = ?
           AND  b.target_table_id = ?
           AND  b.target_taxon_id = ?
           AND  b.target_na_sequence_id = ?
           AND  b.target_external_db_release_id = ?
           AND  u.is_best_score = 1
EOSQL


       $sql = "select * from ($sql) order by target_start asc, target_end asc";

       my $sth = $dbh->prepare($sql) or die "bad sql $sql";
       $sth->execute($taxonId,$qDbRlsId,$targetTableId,$taxonId,$seqId,$tDbRlsId) or die "could not run $sql";

       my %gffFile;

       while (my ($sid, $qSourceId, $s, $e,$tstarts,$blocksizes) = $sth->fetchrow_array) {
	 if ($cinfo->{'end'} == 0 || $s > $cinfo->{'end'}) { #make new cluster, just starting or no overlap
	   $self->makeGffFile(\%gffFile,$seqs->{$seqId},$cinfo->{start},$cinfo->{end}) if $cinfo->{end} > 0;
	   undef %gffFile;
	   my $newCid = $cinfo->{'id'} + 1;
	   $cinfo = {'id'=>$newCid, 'start'=>$s, 'end'=>$e};
	   $clusters->{$newCid}->{$sid} = '';   #clusterHash{cluster number}->{assembly_sequence_id}=''
           $gffFile{$cinfo-{'id'}}->{$qSourceId}->{'tstarts'} = $tstarts;
	   $gffFile{$cinfo-{'id'}}->{$qSourceId}->{'blocksizes'} = $blocksizes;
	 } 
	 else {
	   my $cid = $cinfo->{'id'};   #otherwise, put the sequence in the same cluster and set the end to the larger value
	   $gffFile{$cinfo-{'id'}}->{$qSourceId}->{'tstarts'} = $tstarts;
	   $gffFile{$cinfo-{'id'}}->{$qSourceId}->{'blockSizes'} = $blocksizes;
	   $clusters->{$cid}->{$sid} = '';
	   $cinfo->{'end'} = $e if $e > $cinfo->{'end'};
	 }
       }
       $sth->finish;
     }

    return $clusters;
}

sub getSequencePieces {
   my ($self, $dbh, $tDbRlsId) = @_;

   my $targetTableName = $self->getArg('targetTableName');

   my $testChr = $self->getArg('testChr');

   $testChr = ",$testChr";

   my $sql = "select na_sequence_id, source_id from $targetTableName where external_database_release_id = ?";
   $sql .= " where chromosome = ?" if $testChr;

   my $sth = $dbh->prepare($sql) or die "bad sql $sql";

   $sth->execute($tDbRlsId,$testChr) || die "could not run $sql";

   my %seqs;

   while (my ($naSeqId,$sourceId) = $sth->fetchrow_array) {
     $seqs{$naSeqId}=$sourceId;
   }

   return \%seqs;
}

sub getTableId {
    my ($self, $dbh) = @_;

    my $targetTableName = $self->getArg('targetTableName');

    $targetTableName =~ s/\S+\.//;

    my $sth = $dbh->prepare("select table_id from core.tableinfo where lower(name) = lower(?)");

    $sth->execute($targetTableName);

    my ($tableId) = $sth->fetchrow();

    $sth->finish();

    return $tableId;
}

sub getTaxonId {
  my ($self, $dbh) = @_;

  my $ncbiTaxId = $self->getArg('ncbiTaxId');

  my $sth = $dbh->prepare("select taxon_id from sres.taxon where ncbi_tax_id = ?");

  $sth->execute($ncbiTaxId);

  my ($taxonId) = $sth->fetchrow();

  $sth->finish();

  return $taxonId;
}

sub makeGffFile {
  my ($self, $gffFile, $tSourceId, $clusterStart, $clusterEnd) = @_;

  my $clusterId = keys %{$gffFile};

  my $dir = $self->getArg('gffDir');

  open (GFF,"> $dir/${clusterId}.gff") || die "Can't open $dir/${clusterId}.gff for writing\n";

  print GFF "browser position ${tSourceId}:${clusterStart}-$clusterEnd
             track name=Cluster$clusterId description='transcript alignments for Cluster$clusterId' visibility=2";

  foreach my $qSourceId (keys %{$gffFile->{$clusterId}}) {
    my $score = $gffFile->{$clusterId}->{$qSourceId}->{'score'};
    my @tStarts = split (/,/, $gffFile->{$clusterId}->{$qSourceId}->{'tstarts'});
    my @blockSizes = split(/,/, $gffFile->{$clusterId}->{$qSourceId}->{'blocksizes'});
    for (my $i=0;$i<@tStarts;$i) {
      print GFF "$tSourceId\tBLAT\texon\t$tStarts[$i]\t$tStarts[$i]+$blockSizes[$i]\t$score\t.\t.\t$qSourceId\n";
    }
  }

    close GFF;
}

1;

