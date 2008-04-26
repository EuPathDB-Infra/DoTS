package DoTS::DotsBuild::Plugin::AssignSequenceDescription;

## Updates the RNA.name and RNA.description fields..

## makes evidence links to contained mRNA or Similarity...
## need to deal with deleting existing evidence....

##for now since updating all rows, just delete with sql

##  delete Evidence where target_table_id = 192

## Brian Brunk 7/21/00

@ISA = qw(GUS::PluginMgr::Plugin);
use GUS::PluginMgr::Plugin;
use strict;
use GUS::Model::DoTS::RNA;
use GUS::Model::DoTS::Similarity;
use GUS::Model::DoTS::ExternalNASequence;

$| = 1;

sub new {
  my ($class) = @_;

  my $self = {};

  bless($self,$class);;

  my $usage = 'Assigns description to sequence based on similarity to nrdb neighbors';

  my $argsDeclaration =
    [integerArg({name => 'testnumber',
		 descr => 'number of iterations for testing',
		 constraintFunc => undef,
		 reqd => 0,
		 isList => 0
		}),
     stringArg({name => 'restart',
		descr => 'sql string that returns primary_key list from --table to ignore',
		constraintFunc => undef,
		reqd => 0,
		isList => 0
	       }),
     booleanArg({name => 'deleteEvidence',
		 descr => 'delete evidence foreach RNA if update (more efficient to version/delete with sql if doing all entries)',
		 constraintFunc => undef,
		 reqd => 0
		}),
     booleanArg({name => 'doNotVersion',
		 descr => 'Sets globalDoNotVersion to 1',
		 constraintFunc => undef,
		 reqd => 0
		}),
     booleanArg({name => 'addEvidence',
		 descr => 'evidence to support description',
		 constraintFunc => undef,
		 reqd => 0
		}),
     booleanArg({name => 'use_mrna',
		 descr => 'use contained_mrna to assign description if Assembly',
		 constraintFunc => undef,
		 reqd => 0
		}),
     stringArg({name => 'idSQL',
		descr => 'SQL statement:  should return table,query_table primary_key list from --table',
		constraintFunc => undef,
		reqd => 1,
		isList => 0
	       }),
     stringArg({name => 'ignoreAlgInv',
		descr => 'includes " and s.row_alg_invocation_id not in (ignoreAlgInv)" in similarity query',
		constraintFunc => undef,
		reqd => 0,
		isList => 0
	       }),
     stringArg({name => 'sim_mod_date',
		descr => 'includes " and s.modification_date >= \'sim_mod_date\'" in similarity query',
		constraintFunc => undef,
		reqd => 0,
		isList => 0
	       }),
     stringArg({name => 'query_table',
		descr => 'query table which contains similarities (full className)',
		constraintFunc => undef,
		reqd => 1,
		isList => 0
	       }),
     stringArg({name => 'table',
		descr => 'table to  update description (full className)',
		constraintFunc => undef,
		reqd => 1,
		isList => 0
	       }),
     integerArg({name => 'nrdb_ext_db_rls_id',
		 descr => 'id of external db release',
		 constraintFunc => undef,
		 reqd => 1,
		 isList => 0
		}),
     stringArg({name => 'attribute',
		descr => 'name of attribute to udpate',
		constraintFunc => undef,
		reqd => 0,
		isList => 0
	       }),
     stringArg({name => 'dots_mgi_file',
		descr => 'name of a file to output iformation for Carol Bult...hack',
		constraintFunc => undef,
		reqd => 1,
		isList => 0
	       }),
     booleanArg({name => 'update_rna_descriptions',
		 descr => 'Update the RNA.description field',
		 constraintFunc => undef,
		 reqd => 0
		 }),
     integerArg({name => 'taxon_id',
		 descr => 'taxon_id',
		 constraintFunc => undef,
		 reqd => 1,
		 isList => 0
		}),

     booleanArg({name => 'copy_manual_descriptions',
		 descr => 'copy the  manually reviewed descriptions back to Assembly',
		 constraintFunc => undef,
		 reqd => 0
		})];
my $purposeBrief = <<PURPOSEBRIEF;
Assigns description to sequence based on similarity to nrdb neighbors
PURPOSEBRIEF

  my $purpose = <<PLUGIN_PURPOSE;
Assigns description to sequence based on similarity to nrdb neighbors
PLUGIN_PURPOSE

  #check the documentation for this'Assigns description to sequence based on similarity to nrdb neighbors'
  my $tablesAffected = [];

  my $tablesDependedOn = [];

  my $howToRestart = <<PLUGIN_RESTART;
Use string argumnet restart
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



my $simStmt;

sub run {
  my ($self)   = @_;

  my $query_table = $self->getArg('query_table');

  print "Testing on " . $self->getArg('testnumber') if $self->getArg('testnumber');

  my $dbh = $self->getQueryHandle();

  $self->setGlobalNoVersion(1) if $self->getArg('doNotVersion');

  my $stmt;

  my %ignore;
  if ($self->getArg('restart')) {
    $stmt = $dbh->prepare($self->getArg('restart'));
    $stmt->execute();
    while (my($id) = $stmt->fetchrow_array()) {
      $ignore{$id} = 1;
    }
    print "Restarting:  ignoring ".scalar(keys%ignore). " rnas modified since ". $self->getArg('restart') . "\n";
  }

  if($self->getArg('dots_mgi_file')){
    open(MGI,">>" . $self->getArg('dots_mgi_file'));
  }

  $stmt = $dbh->prepare($self->getArg('idSQL'));
  $stmt->execute();
  my @todo;
  my $ctout = 0;
  while (my($id,$qid) = $stmt->fetchrow_array()) {
    next if exists $ignore{$id};
    $ctout++;
    if ($self->getArg('testnumber') && $ctout > $self->getArg('testnumber')) {
      $stmt->finish();
      last;
    }	
    print STDERR "getting ids: $ctout\n" if $ctout %10000 == 0;
    push(@todo,[$id,$qid]); 
  }
  my $totalToDo = scalar(@todo);
  print "Processing $totalToDo " . $self->getArg('table') . "objects\n";
  undef %ignore;                ##free memory.

  my $query_table_id = $self->getTableIdFromTableName($query_table);
  my $table_pk = $self->getTablePKFromTableId($self->getTableIdFromTableName($self->getArg('table')));
  eval("require " . $self->getArg('table'));

  my $protQuery = "select s.similarity_id,ea.aa_sequence_id,edr.external_database_release_id,ea.source_id,ea.name,ea.description,s.number_identical,s.total_match_length,s.pvalue_mant,s.pvalue_exp,ea.length,s.min_subject_start,s.max_subject_end,s.number_of_matches
      from DoTS.Similarity s, DoTS.ExternalAASequence ea,
           SRes.ExternalDatabaseRelease edr
      where query_table_id = $query_table_id 
      and ea.external_database_release_id = edr.external_database_release_id
      and query_id = ?
      and s.subject_table_id = 228";

  $protQuery .= " and s.row_alg_invocation_id not in " . $self->getArg('ignoreAlgInv') if $self->getArg('ignoreAlgInv');
  $protQuery .= " and s.modification_date >= " . $self->getArg('sim_mod_date') if $self->getArg('sim_mod_date');

  $protQuery .= " and s.subject_id = ea.aa_sequence_id
      and ea.external_database_release_id = ". $self->getArg('nrdb_ext_db_rls_id') .
      " order by s.pvalue_exp,s.pvalue_mant,s.number_identical/s.total_match_length desc";

  print STDERR "debug: protQuery: $protQuery\n" if $self->getArg('verbose');

  my $protStmt = $dbh->prepare($protQuery);

  ##need a query to get Similarities if protein as there could be duplicate rows...
  my $simQuery = "select * from DoTS.Similarity where
    query_table_id = $query_table_id 
    and query_id = ?
    and subject_table_id = 228
    and subject_id = ?
    order by similarity_id desc";  ##want the newest first...

      print STDERR "debug: simQuery: $simQuery\n" if $self->getArg('verbose');

  $simStmt = $dbh->prepare($simQuery);

  my $mRNAQuery = "select edr.external_database_release_id,e.source_id,e.name,e.description 
  from Dots.AssemblySequence aseq, Dots.ExternalNASequence e,
       SRes.ExternalDatabaseRelease edr
  where aseq.assembly_na_sequence_id = ? 
  and e.external_database_release_id = edr.external_database_release_id
  and e.na_sequence_id =  aseq.na_sequence_id 
  and e.sequence_type_id in (2,7) 
  order by e.length desc";

  print STDERR "debug: mRNAQuery: $mRNAQuery\n" if $self->getArg('verbose');
  my $mRNAStmt = $dbh->prepare($mRNAQuery);

  ##main loop
  my $count = 0;
  foreach my $r (@todo) {
    my($id,$qid) = @{$r};
    ##following must be in loop to allow garbage collection...
    $self->undefPointerCache();
    $count++;

    print STDERR "debug: Processing $id: $qid\n" if $self->getArg('verbose');
		
    print "Processing $id: $count finished, " . ($totalToDo - $count) . " remaining\n" if $count % 10 == 0;

    ##first get the object...
    my $className = "GUS::Model::" . $self->getArg('table');
    my $obj = $className->new({ $table_pk => $id });
    if (!$obj->retrieveFromDB()) {
      print "ERROR: unable to retrieve $id from database\n";
      next;
    }

    ##need to get any current evidence and mark deleted...
    if ($self->getArg('deleteEvidence')) {
      $obj->retrieveAllEvidenceFromDB();
      foreach my $e ($obj->getAllEvidence()) {
        next unless $e->getAttributeName() eq $self->getArg('attribute');
        print STDERR "debug: marking evidence ".$e->getId()." deleted\n" if $self->getArg('verbose');
        $e->markDeleted();
      }
    }

    my $haveName = 0;

    ##if is an assembly then want to set name from contained mRNA if possible.
    if($self->getArg('use_mrna') && $self->getArg('table') eq "Assembly" && $obj->getContainsMrna()){
      print STDERR "Contains mRNA...getting name of longest one..\n" if $self->getArg('verbose');
      $mRNAStmt->execute($id);
      while (my($ext_db_rel_id,$accession,$gusName,$gusDesc) = $mRNAStmt->fetchrow_array()) {
        print STDERR "debug: mrnaQuery:($ext_db_rel_id,$accession,$gusName,$gusDesc)\n" if $self->getArg('verbose');
        $gusDesc =~ s/^\s*(\S.*\S)\s*$/$1/;
        if($gusDesc){
          $self->updateAssemblyFromRNA($obj,$ext_db_rel_id,$accession,$gusName,$gusDesc);
          $mRNAStmt->finish();
          $haveName = 1;
          print MGI "DT.$id: $ext_db_rel_id|$accession - contained mRNA\n" if $self->getArg('dots_mgi_file');
          last;
        }
      }
      
    }
    if ($haveName) { next; }

    $protStmt->execute($qid);
    print STDERR "debug: Retrieving protein Similarities\n" if $self->getArg('verbose');
    while (my($sim_id,$aaSeqid,$ext_db_rel_id,$accession,$pname,$pdescription,$ident,$ml,$mant,$exp,$aa_length,$sub_start,$sub_end,$num_matches) = $protStmt->fetchrow_array()) {

      print STDERR "debug:   protSim($aaSeqid,$ext_db_rel_id,$accession,$pname,$pdescription,$ident,$ml)\n" if $self->getArg('verbose');
      next if $pdescription =~ /warning/i; ##alu sequence...
      if ($self->updateNameFromSP($obj,$qid,$aaSeqid,$pname,$pdescription,$ident,$ml,$sim_id,$aa_length,$sub_start,$sub_end,$num_matches)) {
        print MGI "DT.$id: $ext_db_rel_id|$accession, pVal=",$self->getPValue($mant,$exp),", ML=$ml, %=", $self->roundOff(($ident/$ml)*100),"\n" if $self->getArg('dots_mgi_file');
        $protStmt->finish();
        $haveName = 1;
        last;
      }
    }
    if ($haveName) { next; }

    $obj->set($self->getArg('attribute'),"No NR protein Similarities") unless $obj->get($self->getArg('attribute')) eq "No NR protein Similarities";
    $obj->submit() if $obj->hasChangedAttributes();
  }
  
  my ($updateRNA, $manualDescriptions);


  if ($self->getArg('copy_manual_descriptions')) {
      $manualDescriptions = $self->manualDescriptions($dbh);
  }

  if ($self->getArg('update_rna_descriptions')) {
      $updateRNA = $self->updateRNA($dbh);
  }

  $self->closeQueryHandle();
  ############################################################
  ###  put an informative summary in the results variable
  ############################################################
  my $results = "Processed $count Entries, updated $updateRNA RNA descriptions, copied $manualDescriptions manual descriptions to Assembly";

  print STDERR "\n$results\n";

  return "$results";
}

sub updateAssemblyFromRNA {
  my($self,$assem,$ext_db_rls_id,$gacc,$gname,$gdesc) = @_;
  print STDERR "updateNameFromRNA: ($assem,$gacc,$gname,$gdesc)\n" if $self->getArg('verbose');
  ##first update name and desc if different from above
##  $assem->set('name',$gname) if $gname ne $assem->get('name');
  my $desc = substr($gdesc,0,255);
  $assem->set('description',$desc) if $desc ne $assem->get('description');
  
  ##should do the evidence stuff regardless of whether is modified this time as don't have at all!
  ##external_db_id could be 78 or 22
  my $fact = GUS::Model::DoTS::ExternalNASequence->
    new({ 'external_database_release_id' => $ext_db_rls_id,
	  'source_id' => $gacc });
  
  ##if have retrieved fact tehn add evidence..
  if ($self->getArg('addEvidence')) {
    if($fact->retrieveFromDB()){
      $assem->addEvidence($fact,1,"description");
    }
  }
  
  ##now submit the RNA...will only be updated if sets were done...
  $assem->submit() if $assem->hasChangedAttributes();
}

sub updateNameFromSP {
  my($self,$obj,$qid,$aaSeqId,$name,$desc,$ident,$ml,$sim_id,$aa_length,$sub_start,$sub_end,$num_matches) = @_;
  return undef unless $desc;    ##don't want to do this on if no description

  my $pcid = $self->roundOff(($ident/$ml)*100);
  my $pcCov;
  
  ##if the matchlength > sub_end - sub_start (overlapping..) want to compute region covered..
  if ($ml > $sub_end - $sub_start + 1 || $self->getArg('addEvidence')) {
    print STDERR "\nEXTRACTING \% coverage from similarity...num_matches = $num_matches\n\n" if $self->getArg('verbose');
    my $simF = GUS::Model::DoTS::Similarity->new({'similarity_id' => $sim_id});
    $simF->retrieveFromDB();
    $obj->addEvidence($simF,1,$self->getArg('attribute')) if $self->getArg('addEvidence');
    if($ml > $sub_end - $sub_start + 1){
      print STDERR "MatchLength=$ml\n" if $self->getArg('verbose');
      $pcCov = $self->roundOff(($simF->getSubjectLengthCovered()/$aa_length)*100);
    }
  }
  $pcCov = $pcCov ? $pcCov : $self->roundOff(($ml/$aa_length)*100);
  $pcCov = $pcCov > 100 ? 100 : $pcCov;
  my $de = "$pcid\% identity to $pcCov\% of $desc";
  my $description = substr($de,0,255);
  print STDERR "\nDT.",$obj->getId()," DESC: $description\n" if $self->getArg('verbose');
  $obj->set($self->getArg('attribute'),$description) unless $description eq $obj->get($self->getArg('attribute'));
  $obj->submit() if $obj->hasChangedAttributes();
  return 1;
}

sub roundOff {
  my($self,$num) = @_;
  my $floor = int($num);
  if (($num - $floor) < 0.5) {
    return $floor;
  } else {
    return $floor + 1;
  }
}

sub getPValue {
  my($self,$mant,$exp) = @_;
  return $mant . (($exp != -999999 && $exp != 0) ? "e" . $exp : "");
}


sub updateRNA {
    my ($self,$dbh) = @_;
    my $rows = $dbh->prepareAndExecute ("update dots.RNA set description = substr ((select a.description from dots.rnafeature rf, dots.rnainstance rs, dots.assembly a where rs.rna_id = dots.rna.rna_id and rf.na_feature_id = rs.na_feature_id and a.na_sequence_id = rf.na_sequence_id ),0,255)where rna_id in (select rs1.rna_id from dots.rnainstance rs1, dots.rnafeature rf1, dots.assembly a1 where a1.taxon_id = " . $self->getArg('taxon_id') . " and a1.na_sequence_id = rf1.na_sequence_id and rs1.na_feature_id = rf1.na_feature_id) and (rna_id in ( select c.rna_id from dots.rnarnacategory c where c.rna_category_id != 17) or rna_id not in (select c.rna_id from dots.rnarnacategory c) or rna_id in (select c.rna_id from dots.rnarnacategory c where c.rna_category_id = 17 and (dots.rna.description like '%identity%' or dots.rna.description like 'No NR%')))"); 
    if ($self->getArg('commit')) { $dbh->commit;}
    return $rows;
}


sub manualDescriptions {
    my ($self,$dbh) = @_;
    my $rows = $dbh->prepareAndExecute ("update dots.assembly set description = 
     (select r.description from dots.RNA r, dots.rnainstance rs, dots.rnafeature rf
     where rf.na_sequence_id = dots.assembly.na_sequence_id
     and rs.na_feature_id = rf.na_feature_id and rs.rna_id = r.rna_id )
     where taxon_id = " . $self->getArg('taxon_id')
     . " and na_sequence_id in ( select rf1.na_sequence_id from dots.rnafeature rf1, dots.rnainstance rs1, dots.rna r1
     where r1.review_status_id = 1 and r1.description != 'No NR protein Similarities' and rs1.rna_id = r1.rna_id
     and rf1.na_feature_id = rs1.na_feature_id)");
    if ($self->getArg('commit')) { $dbh->commit;}
    return $rows;
}



1;

__END__

