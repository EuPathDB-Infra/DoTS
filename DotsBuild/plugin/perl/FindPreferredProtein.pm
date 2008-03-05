package DoTS::DotsBuild::Plugin::FindPreferredProtein;

@ISA = qw(GUS::PluginMgr::Plugin);
use GUS::PluginMgr::Plugin;
use strict;



$| = 1;

sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self,$class);
    
    my $usage = 'Plug_in to determine the preferred assembly translation, GenBank assigned protein sequence for a contained mRNA preferred to a FrameFinder translation.';

 my  $argsDeclaration =
    [ integerArg({name => 'refseq_db_rel_id',
		  descr => 'external_database_release_id for RefSeq in dots.ExternalNASequence',
		  constraintFunc => undef,
		  reqd => 1,
		  isList => 0
		 }),
      integerArg({name => 'taxon_id',
		  descr => 'taxon_id',
		  constraintFunc => undef,
		  reqd => 1,
		  isList => 0
		 }),
      integerArg({name => 'swissprot_db_rel_id',
		 descr => 'external_database_release_id for SwissProt in dots.NRDBEntry',
		 constraintFunc => undef,
		 reqd => 1,
		 isList => 0
		}),
    ];
 my $purposeBrief = <<PURPOSEBRIEF;
Determine the preferred assembly translation
PURPOSEBRIEF

  my $purpose = <<PLUGIN_PURPOSE;
Determine the preferred assembly translation, GenBank assigned protein sequence for a contained mRNA preferred to a FrameFinder translation.
PLUGIN_PURPOSE

  #check the documentation for this
  my $tablesAffected = [['DoTS::ProteinInstance', '']
		       ];

  my $tablesDependedOn = [
			  ['DoTS::ProteinInstance', '']
			 ];

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

  $self->initialize({requiredDbVersion => 3.5,
		     cvsRevision => '$Revision$',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation,
		    });
  return $self;
}


sub run {
    my $self   = shift;

    my $dbh = $self->getQueryHandle();

    $self->resetIsReferenceToFalse($dbh);

    $self->setIsReferenceToTrue($dbh); #set proteininstance.is_reference = 1 for assemblies without mRNA

    my $assemblies = $self->getAssemblies($dbh); #get assemblies with mRNAs

    $self->updateProteinInstance($dbh,$assemblies);#set proteininstance.is_reference = 1 for assemblies with mRNA:RefSeq (longest protein)>swissprot(longest protein)>the longest protein>FF translation if mRNA doesn't have a translation in the db

    my $returnStmt = 'set proteininstance.is_reference = 1 for assemblies with mRNA:RefSeq>swissprot>the longest protein>FF translation';

    return $returnStmt;
}

sub resetIsReferenceToFalse {
  my ($self,$dbh) = @_;
  my $taxon = $self->getArg('taxon_id');

  my $updateSQL = "update dots.proteininstance set is_reference = 0 where protein_instance_id in (select /*+ RULE */ pi.protein_instance_id from dots.proteininstance pi, dots.protein p, dots.rnainstance ri,dots.rnafeature rf, dots.assembly a where a.taxon_id = $taxon and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = ri.na_feature_id and ri.rna_id = p.rna_id and p.protein_id = pi.protein_id)";

  $self->log ("Updating is_reference to 0 for assemblies with taxon_id = $taxon\n");

  my $num = $dbh->do($updateSQL) || die "Insert failed.\nSQL: $updateSQL\n";

  $self->log ("Committing $num inserts\n") if $self->getArg('commit');

  $dbh->commit() if $self->getArg('commit');
}


sub setIsReferenceToTrue {
  my ($self,$dbh) = @_;
  my $taxon = $self->getArg('taxon_id');

  my $updateSQL = "update dots.proteininstance set is_reference = 1 where protein_instance_id in (select /*+ RULE */ p.protein_instance_id from dots.proteininstance p, dots.translatedaafeature tf, dots.rnafeature rf, dots.assembly a where a.taxon_id = $taxon and a.contains_mrna = 0 and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_feature_id = p.aa_feature_id)";

  $self->log ("Updating is_reference for assemblies with no mRNA\n");

  my $num = $dbh->do($updateSQL) || die "Insert failed.\nSQL: $updateSQL\n";

  $self->log ("Committing $num inserts\n") if $self->getArg('commit');

  $dbh->commit() if $self->getArg('commit');

}

sub getAssemblies {
  my ($self,$dbh) = @_;

  my %assemblies;

  my $taxon = $self->getArg('taxon_id');

  my $query = "select na_sequence_id from dots.assembly where taxon_id = $taxon and contains_mrna = 1";

  my $stmt = $dbh->prepareAndExecute($query);

  while(my ($assembly) = $stmt->fetchrow_array()){
    $assemblies{$assembly} = 1;
  }
  my $num = scalar (keys %assemblies);
  $self->log ("$num DTs contain mRNA and will have preferred proteins assigned\n");
  return \%assemblies;
}

sub updateProteinInstance {

  my ($self,$dbh,$assemblies) = @_;
  my $refseq_db_rel = $self->getArg('refseq_db_rel_id');
  my $swissprot_db_rel = $self->getArg('swissprot_db_rel_id');
  my $num;
  my $refseqquery = "select p.protein_instance_id from dots.assemblysequence ass, dots.rnafeature rf, dots.nrdbentry n, dots.translatedaafeature tf, dots.externalaasequence xa, dots.proteinInstance p  where ass.assembly_na_sequence_id = ? and ass.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_sequence_id = xa.aa_sequence_id and xa.aa_sequence_id = n.aa_sequence_id and n.external_database_release_id = $refseq_db_rel and tf.aa_feature_id = p.aa_feature_id order by xa.length desc";

  my $stmt =  $dbh->prepare($refseqquery);
  foreach my $id (keys %{$assemblies}) {
    $stmt->execute($id);
    my ($proteininstance) = $stmt->fetchrow_array();
    $stmt->finish();
    if ($proteininstance) {
      $dbh->do("update dots.proteininstance set is_reference = 1 where protein_instance_id = $proteininstance")       || die "Can't update dots.proteininstance row : protein_instance_id = $proteininstance\n";
      $num++;
      delete($assemblies->{$id});
    }
  }

  $self->log("Committing updates for $num RefSeq containing assemblies\n") if $self->getArg('commit');
  $dbh->commit() if $self->getArg('commit');
  $num = 0;
  my $swissprotquery = "select p.protein_instance_id from dots.assemblysequence ass, dots.rnafeature rf, dots.nrdbentry n, dots.translatedaafeature tf, dots.externalaasequence xa, dots.proteinInstance p  where ass.assembly_na_sequence_id = ? and ass.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_sequence_id = xa.aa_sequence_id and xa.aa_sequence_id = n.aa_sequence_id and n.external_database_release_id = $swissprot_db_rel and tf.aa_feature_id = p.aa_feature_id order by xa.length desc";

  my $stmt =  $dbh->prepare($swissprotquery);

  foreach my $id (keys %{$assemblies}) {
    $stmt->execute($id);
    my ($proteininstance) = $stmt->fetchrow_array();
    $stmt->finish();
     if ($proteininstance) {
       $dbh->do("update dots.proteininstance set is_reference = 1 where protein_instance_id = $proteininstance")       || die "Can't update dots.proteininstance row : protein_instance_id = $proteininstance\n";
       $num++;
       delete($assemblies->{$id});
     }
  }

  $self->log("Committing updates for $num assemblies having mRNA with translations in SwissProt\n") if $self->getArg('commit');
  $dbh->commit() if $self->getArg('commit');
  $num = 0;
  my $lengthquery = "select p.protein_instance_id from dots.assemblysequence ass, dots.rnafeature rf, dots.translatedaafeature tf, dots.externalaasequence xa, dots.proteinInstance p  where ass.assembly_na_sequence_id = ? and ass.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_sequence_id = xa.aa_sequence_id and tf.aa_feature_id = p.aa_feature_id order by xa.length desc";

  my $stmt =  $dbh->prepare($lengthquery);

  foreach my $id (keys %{$assemblies}) {
    $stmt->execute($id);
    my ($proteininstance) = $stmt->fetchrow_array();
    $stmt->finish();
    if ($proteininstance) {
      $dbh->do("update dots.proteininstance set is_reference = 1 where protein_instance_id = $proteininstance")       || die "Can't update dots.proteininstance row : protein_instance_id = $proteininstance\n";
      $num++;
      delete($assemblies->{$id});
    }
  }

  $self->log("Committing updates for $num assemblies having mRNA with the longest translations\n") if $self->getArg('commit');
  $dbh->commit() if $self->getArg('commit');
  $num = 0;
  my $query = "select p.protein_instance_id from dots.assembly a, dots.rnafeature rf, dots.translatedaafeature tf, dots.proteinInstance p  where a.na_sequence_id = ? and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_feature_id = p.aa_feature_id";

  my $stmt =  $dbh->prepare($query);

  foreach my $id (keys %{$assemblies}) {
    $stmt->execute($id);
    my ($proteininstance) = $stmt->fetchrow_array();
    $stmt->finish();
    if ($proteininstance) {
      $dbh->do("update dots.proteininstance set is_reference = 1 where protein_instance_id = $proteininstance")       || die "Can't update dots.proteininstance row : protein_instance_id = $proteininstance\n";
      $num++;
      delete($assemblies->{$id});
    }
  }

  $self->log("Committing updates for $num assemblies having mRNA without translations\n") if $self->getArg('commit');
  $dbh->commit() if $self->getArg('commit');

  $num = scalar (keys %{$assemblies});

  if ($num > 0) {
    $self->log ("The following DTs containing mRNA have not been assigned a preferred protein:\n");
    foreach my $id (keys %{$assemblies}) {
      $self->log ("    $id\n");
    }
  }
}

1;

