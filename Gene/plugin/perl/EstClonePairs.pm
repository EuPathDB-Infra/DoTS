package DoTS::Gene::Plugin::EstClonePairs;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::Disp;
use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'look for pairs of ESTs linked image_id or washu_name';
  
    my $purpose = $purposeBrief;

    my $tablesDependedOn = [['DoTS::Est','est'],['DoTS::Clone', 'est clone'], ['DoTS::Library', 'est library'], ['DoTS.AssemblySequence', 'DoTS assembly seqs'], ['DoTS.NAFeature', 'na feature'], ['DoTS.RNAInstance', 'RNA instance'], ['DoTS.RNA', 'RNA']];

    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
Cannot be restarted, yet. 
RESTART

    my $failureCases;

    my $notes = <<NOTES;
=pod

=head2 F<General Description>

Plugin find EST clones linked by common image_id or washu_name, and cach info such as what sDGs the ests belong to

=head1 AUTHORS

Y Thomas Gan

=head1 COPYRIGHT

Copyright, Trustees of University of Pennsylvania 2004. 

=cut

NOTES
    
    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
    my $argsDeclaration  =
    [
     stringArg({name => 'temp_login',
		descr => 'login for temp table usage',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0, 
	    }),

     stringArg({name => 'join_column',
		descr => 'column of DoTS.Clone table on which to join for Est pairs',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0, 
	    }),

    integerArg({name => 'taxon_id',
		descr => 'taxon id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		})
    ];

    $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$',
		     cvsTag => '$Name$',
		     name => ref($self),
		     revisionNotes => '',
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation
		    });
    return $self;
}

sub run {
  my $self = shift;
    
  $self->logAlgInvocationId();
  $self->logCommit();
  $self->logArgs();

  my $dbh = $self->getQueryHandle();
  my $tempLogin = $self->getArgs->{'temp_login'};
  my $taxonId = $self->getArgs->{'taxon_id'};
  my $joinCol = $self->getArgs->{'join_column'};
  my $sql = &getSql($taxonId, $joinCol);

  print STDERR "# running sql:\n# $sql\n";
  my $sth = $dbh->prepareAndExecute($sql);

  my (%pairs, %lib_clone_counts, %cluster_pairs);
  my ($total, $sub_total, $linked_diff_clusters, $no_cluster_1, $no_cluster_2);
  my @clone_pairs;
  while (my @res = $sth->fetchrow_array) {
      my ($seq1, $p1, $lib_clone, $p2, $seq2) = @res;
      $lib_clone_counts{$lib_clone}++;
      $total++;

      my ($k1, $k2) = ("$seq1$p1", "$seq2$p2");
      unless (exists $pairs{"$k1$k2"} or exists $pairs{"$k2$k1"}) {
	  $sub_total++;
	  my ($cid1, $num_seqs1) = &getGeneClusterInfo($dbh, $seq1);
	  my ($cid2, $num_seqs2) = &getGeneClusterInfo($dbh, $seq2);
	  $pairs{"$k1$k2"} = "";
	  push @clone_pairs, [$cid1, $num_seqs1, $seq1, $p1, $lib_clone, $p2, $seq2, $num_seqs2, $cid2];

	  if ($cid1 == -1 and $cid2 == -1) {
	      $no_cluster_2++;
	  } elsif ($cid1 == -1 or $cid2 == -1) {
	      $no_cluster_1++;
	  } else {
	      $cluster_pairs{"$cid1$cid2"} = ""; 
	      $linked_diff_clusters++ if $cid1 != $cid2;
	  }
      }
  }

  # create the result table if it is not there yet, clean old result if there
  my @cols = &cleanOrCreateTempTable($dbh, $tempLogin, $taxonId, $joinCol);
  print STDERR "saving result into db...\n";
  &saveClonePairs($dbh, $tempLogin, $taxonId, lc($joinCol), \@clone_pairs, \@cols);

  my @counts = sort { $a <=> $b } values %lib_clone_counts;
  my $sz = scalar(@counts); 
  my $min = $counts[0];
  my $med = $counts[$sz/2];
  my $max = $counts[$sz-1];
  my $avg = $total/$sz;

  # print statistics
  select STDERR;
  print "# $total taxon $taxonId EST pairs have common $joinCol\n";
  print "# $sub_total pairs printed out after removing symmetric duplicates\n";
  print "# distribution of number of ESTs with common $joinCol:\n";
  print "# min = $min, med = $med, avg = $avg, max = $max\n"; 
  print "\n";
  print "# ", scalar(keys %cluster_pairs), " taxon $taxonId gene cluster pairs have common $joinCol\n";
  print "# $linked_diff_clusters taxon $taxonId pairs invoving different gene clusters have common $joinCol\n";
  print "# $no_cluster_1 taxon $taxonId EST pairs with common $joinCol do not have gene cluster for one EST\n";
  print "# $no_cluster_2 taxon $taxonId EST pairs with common $joinCol do not have gene cluster for either EST\n";
  select STDOUT;
}

##########################

sub cleanOrCreateTempTable {
  my ($dbh, $tmpLogin, $taxonId, $joinCol, $cols) = @_;

  my $tmp = uc($tmpLogin);
  my $lcJoinCol = lc($joinCol);

  my @cols = ('taxon_id',
	      'join_column',
	      'gene_id_1',
              'num_seqs_1',
              'seq_id_1',
              'pend_1',
              'lib_clone',
              'pend_2',
              'seq_id_2',
              'num_seqs_2',
              'gene_id_2');

  my $sql = "select count(*) from all_tables "
          . "where table_name = 'ESTCLONEPAIR' and owner = '${tmp}'";
  my $sth = $dbh->prepare($sql);
  $sth->execute;
  my ($c) = $sth->fetchrow_array;
  if ($c) {
      print STDERR "# cleaning up ${tmp}.EstClonePair table...\n";
      $dbh->sqlexec("delete ${tmp}.EstClonePair where taxon_id = $taxonId and join_column = '$lcJoinCol'");
  } else {
    $sql = <<BAS_SQL;
CREATE TABLE ${tmp}.EstClonePair(
    $cols[0] NUMBER(12),
    $cols[1] VARCHAR(32),
    $cols[2] NUMBER(12),
    $cols[3] NUMBER(12),
    $cols[4] NUMBER(12),
    $cols[5] VARCHAR(2),
    $cols[6] VARCHAR(64),
    $cols[7] VARCHAR(2),
    $cols[8] NUMBER(12),
    $cols[9] NUMBER(12),
    $cols[10] NUMBER(12)
    /* FOREIGN KEY ($cols[2]) REFERENCES DoTS.BlatAlignment */
)
BAS_SQL

    print STDERR "# creating ${tmp}.EstClonePair table...\n";
    $dbh->sqlexec($sql);
    $dbh->sqlexec("GRANT SELECT ON ${tmp}.EstClonePair to PUBLIC");
    $dbh->sqlexec("create index ESTCLONEPAIR_IND01 on ${tmp}.EstClonePair($cols[4])");
    $dbh->sqlexec("create index ESTCLONEPAIR_IND02 on ${tmp}.EstClonePair($cols[8])");
  }

  @cols;
}

sub saveClonePairs {
    my ($dbh, $tmp, $taxonId, $lcJoinCol, $clone_pairs, $cols) = @_;

    foreach my $cp (@$clone_pairs) {
	my $sql = "insert into ${tmp}.EstClonePair (" . join(',', @$cols) . ")"
	    . " values($taxonId, '$lcJoinCol', $cp->[0], $cp->[1], $cp->[2],"
	    . "'$cp->[3]', '$cp->[4]', '$cp->[5]', $cp->[6], $cp->[7], $cp->[8])";
	$dbh->sqlexec($sql);
    }
}

sub getSql {
    my ($tid, $join_col) = @_;

    my $sql = "select s1.na_sequence_id as seq1, s1.p_end as p1, "
        . "c1.library_id || \':\' || c1.$join_col as lib_clone, "
        . "s2.p_end as p2, s2.na_sequence_id as seq2 "
        . "from DoTS.Library l1, DoTS.Clone c1, DoTS.EST s1, "
        . "     DoTS.Library l2, DoTS.Clone c2, DoTS.EST s2 "
        . "where l1.taxon_id = $tid and l2.taxon_id = $tid "
        . "and l1.library_id = c1.library_id and c1.clone_id = s1.clone_id "
        . "and l2.library_id = c2.library_id and c2.clone_id = s2.clone_id "
        . "and c1.$join_col = c2.$join_col and c1.library_id = c2.library_id "
        . "and not s1.na_sequence_id = s2.na_sequence_id ";

    return $sql;
}

sub getGeneClusterInfo {
    my ($dbh, $seqId) = @_;

    my ($cid, $num_seqs) = (-1, -1);

    my ($sql, $sth);

    # get the gene cluster id for this sequence
    $sql = "select r.gene_id "
         . "from DoTS.AssemblySequence aseq, DoTS.NAFeature naf, "
         . "  DoTS.RNAInstance rs, DoTS.RNA r "
         . "where aseq.na_sequence_id = $seqId " 
         . "and naf.na_sequence_id = aseq.assembly_na_sequence_id "
         . "and rs.na_feature_id = naf.na_feature_id "
         . "and r.rna_id = rs.rna_id ";
    $sth = $dbh->prepare($sql);

    $sth->execute();
    if (my ($c) = $sth->fetchrow_array()) { $cid = $c; } 
  
    # get the number of sequences in this gene cluster, if there is such a cluster
    if ($cid && $cid >= 0) {
	$sql = "select count(aseq.na_sequence_id) "
	    . "from DoTS.AssemblySequence aseq, DoTS.NAFeature naf, DoTS.RNAInstance rs, DoTS.RNA r "
	    . "where naf.na_sequence_id = aseq.assembly_na_sequence_id "
	    . "and rs.na_feature_id = naf.na_feature_id "
            . "and r.rna_id = rs.rna_id "
            . "and r.gene_id = $cid";
	$sth = $dbh->prepare($sql) or die "bad sql (seqId=$seqId) $sql:$!\n";
	$sth->execute() or die "could not run sql $sql:$!\n";
	if (($num_seqs) = $sth->fetchrow_array()) { }
    }

    return ($cid, $num_seqs);
}

1;
