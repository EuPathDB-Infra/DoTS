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

    $self->initialize({requiredDbVersion => 3.5,
		     cvsRevision => '$Revision$',
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

  # create the result table if it is not there yet, clean old result if there
  my $tab = &cleanOrCreateTempTable($dbh, $tempLogin, $taxonId, $joinCol);

  my $sql = &getInsertSql($taxonId, $tab, $joinCol);

  print STDERR "# running sql:\n# $sql\n";
  $dbh->sqlexec($sql);

  print "done\n";
}

##########################

sub cleanOrCreateTempTable {
  my ($dbh, $tmpLogin, $taxonId, $joinCol, $cols) = @_;

  my $tmp = uc($tmpLogin);
  my $lcJoinCol = lc($joinCol);

  my @cols = ('taxon_id',
	      'join_column',
	      'dt_id_1',
              'ass_seq_id_1',
              'pend_1',
              'lib_clone',
              'pend_2',
              'ass_seq_id_2',
              'dt_id_2');

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
    $cols[4] VARCHAR(2),
    $cols[5] VARCHAR(64),
    $cols[6] VARCHAR(2),
    $cols[7] NUMBER(12),
    $cols[8] NUMBER(12)
    /* FOREIGN KEY ($cols[2]) REFERENCES DoTS.BlatAlignment */
)
BAS_SQL

    print STDERR "# creating ${tmp}.EstClonePair table...\n";
    $dbh->sqlexec($sql);
    $dbh->sqlexec("GRANT SELECT ON ${tmp}.EstClonePair to PUBLIC");
    $dbh->sqlexec("create index ESTCLONEPAIR_IND04 on ${tmp}.EstClonePair($cols[2],$cols[8],$cols[5])");
  }

  "${tmp}.EstClonePair";
}

sub getInsertSql {
    my ($tid, $tab, $join_col) = @_;

    my $sql = "select $tid as taxon_id, '$join_col' as join_column, "
	. "aseq1.assembly_na_sequence_id as dt_id_1, "
	. "s1.na_sequence_id as ass_seq1, s1.p_end as p1, "
        . "c1.library_id || \':\' || c1.$join_col as lib_clone, "
        . "s2.p_end as p2, s2.na_sequence_id as ass_seq2, "
	. "aseq2.assembly_na_sequence_id as dt_id_2 "
        . "from DoTS.Library l1, DoTS.Clone c1, DoTS.EST s1, DoTS.AssemblySequence aseq1, "
        . "     DoTS.Library l2, DoTS.Clone c2, DoTS.EST s2, DoTS.AssemblySequence aseq2 "
        . "where l1.taxon_id = $tid and l2.taxon_id = $tid "
        . "and l1.library_id = c1.library_id and c1.clone_id = s1.clone_id "
        . "and l2.library_id = c2.library_id and c2.clone_id = s2.clone_id "
        . "and c1.$join_col = c2.$join_col and c1.library_id = c2.library_id "
	. "and s1.na_sequence_id = aseq1.na_sequence_id "
	. "and s2.na_sequence_id = aseq2.na_sequence_id "
        . "and not s1.na_sequence_id = s2.na_sequence_id ";

    return "insert into $tab ($sql)";
}

1;
