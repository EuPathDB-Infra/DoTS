package DoTS::DotsBuild::Plugin::FindPreferredProtein;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;



$| = 1;

sub new {
    my ($class) = @_;
    
    my $self = {};
    bless($self,$class);
    
    my $usage = 'Plug_in to determine the preferred assembly translation, GenBank assigned protein sequence for a contained mRNA preferred to a FrameFinder translation.';
    
    my $easycsp =
	[{o => 'taxon_id',
	  t => 'int',
	  h => 'the taxon_id of the assemblies to be used',
         },
	 {o => 'refseq_db_rel_id',
	  t => 'string',
	  h => 'external_database_release_id for RefSeq in dots.ExternalNASequence'
         },
	 {o => 'swissprot_db_rel_id',
	  t => 'string',
	  h => 'external_database_release_id for SwissProt in dots.NRDBEntry'
         }];
    
    $self->initialize({requiredDbVersion => {},
		       cvsRevision => '$Revision$',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		       name => ref($self),
		       revisionNotes => 'make consistent with GUS 3.0',
		       easyCspOptions => $easycsp,
		       usage => $usage
		       });
    
    return $self;
}


sub run {
    my $self   = shift;

    $self->log ($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");

    unless ($self->getArgs()->{'taxon_id'}) {
	die "you must provide a taxon_id\n";
    }
    
    my $dbh = $self->getQueryHandle();

    unless ($self->getArgs()->{'swissprot_db_rel_id'}) { 
      die "you must provide an ext_db_rel_id for SwissProt in dots.nrdbentry";
    }
    unless ($self->getArgs()->{'refseq_db_rel_id'}) { 
      die "you must provide an ext_db_rel_id for RefSeq mRNA";
    } 

    $self->setIs_Reference($dbh); #set proteininstance.is_reference = 1 for assemblies without mRNA

    my $assemblies = $self->getAssemblies($dbh); #get assemblies with mRNAs

    $self->updateProteinInstance($dbh,$assemblies);#set proteininstance.is_reference = 1 for assemblies with mRNA:RefSeq (longest protein)>swissprot(longest protein)>the longest protein>FF translation if mRNA doesn't have a translation in the db
} 
    

sub setIs_Reference {
  my ($self,$dbh) = @_;
  my $taxon = $self->getArgs()->{'taxon_id'};

  my $updateSQL = "update dots.proteininstance set is_reference = 1 where protein_instance_id in (select /*+ RULE */ p.protein_instance_id from dots.proteininstance p, dots.translatedaafeature tf, dots.rnafeature rf, dots.assembly a where a.taxon_id = $taxon and a.contains_mrna = 0 and a.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_feature_id = p.aa_feature_id)";

  $self->log ("Updating is_reference for assemblies with no mRNA\n");

  $dbh->do($updateSQL) || die "Insert failed.\nSQL: $updateSQL\n";

  $self->log ("Committing insert\n") if $self->getArgs()->{'commit'} ;

  $dbh->commit() if $self->getArgs()->{'commit'};

}

sub getAssemblies {
  my ($self,$dbh) = @_;

  my %assemblies;

  my $taxon = $self->getArgs()->{'taxon_id'};

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
  my $refseq_db_rel = $self->getArgs()->{'refseq_db_rel_id'};
  my $swissprot_db_rel = $self->getArgs()->{'swissprot_db_rel_id'}; 
  my $num;
  my $refseqquery = "select p.protein_instance_id from dots.assemblysequence ass, dots.rnafeature rf, dots.externalnasequence x, dots.translatedaafeature tf, dots.externalaasequence xa, dots.proteinInstance p  where ass.assembly_na_sequence_id = ? and ass.na_sequence_id = x.na_sequence_id and x.external_database_release_id = $refseq_db_rel and ass.na_sequence_id = rf.na_sequence_id and rf.na_feature_id = tf.na_feature_id and tf.aa_sequence_id = xa.aa_sequence_id and tf.aa_feature_id = p.aa_feature_id order by xa.length desc";

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

  $self->log("Committing updates for $num RefSeq containing assemblies\n") if $self->getArgs()->{'commit'} ;
  $dbh->commit() if $self->getArgs()->{'commit'} ;
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

  $self->log("Committing updates for $num assemblies having mRNA with translations in SwissProt\n") if $self->getArgs()->{'commit'} ;
  $dbh->commit() if $self->getArgs()->{'commit'} ;
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

  $self->log("Committing updates for $num assemblies having mRNA with the longest translations\n") if $self->getArgs()->{'commit'} ;
  $dbh->commit() if $self->getArgs()->{'commit'} ;
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

  $self->log("Committing updates for $num assemblies having mRNA without translations\n") if $self->getArgs()->{'commit'} ;
  $dbh->commit() if $self->getArgs()->{'commit'} ;

  $num = scalar (keys %{$assemblies});

  if ($num > 0) {
    $self->log ("The following DTs containing mRNA have not been assigned a preferred protein:\n");
    foreach my $id (keys %{$assemblies}) {
      $self->log ("    $id\n");
    }
  }
}

1;

