############################################################
## Change Package name....
############################################################
package AssignSequenceDescription;

## Updates the RNA.name and RNA.description fields..

## makes evidence links to contained mRNA or Similarity...
## need to deal with deleting existing evidence....

##for now since updating all rows, just delete with sql

##  delete Evidence where target_table_id = 192

## Brian Brunk 7/21/00

use strict;
use DBI;
############################################################
# Add any specific objects (GUSdev::) here
############################################################
use Objects::GUSdev::RNA;
use Objects::GUSdev::Similarity;
use Objects::GUSdev::ExternalNASequence;



sub new {
	my $Class = shift;

	return bless {}, $Class;
}

sub Usage {
	my $M   = shift;
	return 'Assigns description to sequence based on similarity to nrdb neighbors';
}

############################################################
# put the options in this method....
############################################################
sub EasyCspOptions {
	my $M   = shift;
	{

 testnumber        => {
                       o => 'testnumber=i',
                       h => 'number of iterations for testing',
                      },

 restart           => {
                       o => 'restart=s',
                       h => 'sql string that returns primary_key list from --table to ignore',
                      },

 deleteEvidence    => {
                       o => 'deleteEvidence!',
                       h => 'delete evidence foreach RNA if update (more efficient to version/delete with sql if doing all entries)',
                      },
 doNotVersion      => {
                       o => 'doNotVersion!',
                       h => 'Sets globalDoNotVersion to 1',
                      },
 addEvidence       => {
                       o => 'addEvidence!',
                       h => 'evidence to support description',
                      },
 use_mrna          => {
                       o => 'use_mrna',
                       t => 'boolean',
                       h => 'use contained_mrna to assign description if Assembly',
                      },
 idSQL             => {
                       o => 'idSQL=s',
                       h => 'SQL statement:  should return table,query_table primary_key list from --table',
                      },
  ignoreAlgInv      => {
                        o => 'ignoreAlgInv',
                        t => 'string',
                        h => 'includes " and s.row_alg_invocation_id not in (ignoreAlgInv)" in similarity query',
                       },
  sim_mod_date     => {
                        o => 'sim_mod_date',
                        t => 'string',
                        h => 'includes " and s.modification_date >= \'sim_mod_date\'" in similarity query',
                       },
 query_table       => {
                       o => 'query_table=s',
                       h => 'query table which contains similarities',
                      },
 table             => {
                       o => 'table=s',
                       h => 'table to  update description',
                      },
 nrdb_ext_db_id    => {
                       o => 'nrdb_ext_db_id',
                       t => 'int',
                       h => 'table to  update description',
                       d => '4194',
                      },
 attribute         => {
                       o => 'attribute',
                       t => 'string',
                       h => 'name of attribute to udpate',
                       d => 'description',
                      },
 dots_mgi_file      => {
                       o => 'dots_mgi_file',
                       t => 'string',
                       h => 'name of a file to output iformation for Carol Bult...hack',
                      },
 update_rna_descriptions   => {
                       o => 'update_rna_descriptions!',
                       h => 'Update the RNA.description field',
		      },
 taxon_id            => {
                       o => 'taxon_id',
                       t => 'int',
                       h => 'taxon_id, required for update_rna_descriptions and copy_manual_descriptions',
		      },
 copy_manual_descriptions   => {
                       o => 'copy_manual_descriptions!',
                       h => 'scopy the  manually reviewed descriptions back to Assembly',
		      },
                    }
}

my $ctx;
my $debug = 0;
$| = 1;
my $simStmt;

sub Run {
  my $M   = shift;
  $ctx = shift;

  die "--table and --idSQL are required\n" unless $ctx->{cla}->{table} && $ctx->{cla}->{idSQL};
  my $query_table = $ctx->{cla}->{query_table} ?  $ctx->{cla}->{query_table} : $ctx->{cla}->{table};

  if ($ctx->{cla}->{update_rna_descriptions} || $ctx->{cla}->{copy_manual_descriptions}) {
      die "taxon_id required for updating RNA descriptions and copying manual descriptions\n" unless $ctx->{cla}->{taxon_id};
  }
  print $ctx->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $ctx->{'testnumber'}\n" if $ctx->{'testnumber'};
  #  print "Setting globalNoVersion(1)\n";
  #  $ctx->{self_inv}->setGlobalNoVersion(1);

  my $dbh = $ctx->{self_inv}->getQueryHandle();
  $dbh->{AutoCommit}=0;
  $ctx->{self_inv}->setGlobalNoVersion(1) if $ctx->{cla}->{doNotVersion};

  my $stmt;

  my %ignore;
  if ($ctx->{'restart'}) {
    $stmt = $dbh->prepare($ctx->{cla}->{restart});
    $stmt->execute();
    while (my($id) = $stmt->fetchrow_array()) {
      $ignore{$id} = 1;
    }
    print "Restarting:  ignoring ".scalar(keys%ignore). " rnas modified since $ctx->{'restart'}\n";
  }

  if($ctx->{cla}->{dots_mgi_file}){
    open(MGI,">>$ctx->{cla}->{dots_mgi_file}");
  }

  $stmt = $dbh->prepare($ctx->{cla}->{idSQL});
  $stmt->execute();
  my @todo;
  my $ctout = 0;
  while (my($id,$qid) = $stmt->fetchrow_array()) {
    next if exists $ignore{$id};
    $ctout++;
    if ($ctx->{cla}->{testnumber} && $ctout > $ctx->{cla}->{testnumber}) {
      $stmt->finish();
      last;
    }	
    print STDERR "getting ids: $ctout\n" if $ctout %10000 == 0;
    push(@todo,[$id,$qid]); 
  }
  my $totalToDo = scalar(@todo);
  print "Processing $totalToDo $ctx->{cla}->{table} objects\n";
  undef %ignore;                ##free memory.

  my $query_table_id = $ctx->{self_inv}->getTableIdFromTableName($query_table);
  my $table_pk = $ctx->{self_inv}->getTablePKFromTableId($ctx->{self_inv}->getTableIdFromTableName($ctx->{cla}->{table}));
  eval("require Objects::".$ctx->{self_inv}->getDatabase()->getDbName()."::$ctx->{cla}->{table}");

  my $protQuery = "select s.similarity_id,ea.aa_sequence_id,ea.external_db_id,ea.source_id,ea.name,ea.description,s.number_identical,s.total_match_length,s.pvalue_mant,s.pvalue_exp,ea.length,s.min_subject_start,s.max_subject_end,s.number_of_matches
      from Similarity s, ExternalAASequence ea
      where query_table_id = $query_table_id 
      and query_id = ?
      and s.subject_table_id = 83";

  $protQuery .= " and s.row_alg_invocation_id not in ($ctx->{cla}->{ignoreAlgInv})" if $ctx->{cla}->{ignoreAlgInv};
  $protQuery .= " and s.modification_date >= '$ctx->{cla}->{sim_mod_date}'" if $ctx->{cla}->{sim_mod_date};

  $protQuery .= " and s.subject_id = ea.aa_sequence_id
      and ea.external_db_id = $ctx->{cla}->{nrdb_ext_db_id}
      order by s.pvalue_exp,s.pvalue_mant,s.number_identical/s.total_match_length desc";

  print STDERR "debug: protQuery: $protQuery\n" if $debug || $ctx->{cla}->{verbose};

  my $protStmt = $dbh->prepare($protQuery);

  ##need a query to get Similarities if protein as there could be duplicate rows...
  my $simQuery = "select * from Similarity where
    query_table_id = $query_table_id 
    and query_id = ?
    and subject_table_id = 83
    and subject_id = ?
    order by similarity_id desc";  ##want the newest first...

      print STDERR "debug: simQuery: $simQuery\n" if $debug || $ctx->{cla}->{verbose};

  $simStmt = $dbh->prepare($simQuery);

  my $mRNAQuery = "select e.external_db_id,e.source_id,e.name,e.description 
  from AssemblySequence aseq, ExternalNASequence e 
  where aseq.assembly_na_sequence_id = ? 
  and e.na_sequence_id =  aseq.na_sequence_id 
  and e.sequence_type_id in (2,7) 
  order by e.length desc";

  print STDERR "debug: mRNAQuery: $mRNAQuery\n" if $debug || $ctx->{cla}->{verbose};
  my $mRNAStmt = $dbh->prepare($mRNAQuery);

  ##main loop
  my $count = 0;
  foreach my $r (@todo) {
    my($id,$qid) = @{$r};
    ##following must be in loop to allow garbage collection...
    $ctx->{'self_inv'}->undefPointerCache();
    $count++;

    print STDERR "debug: Processing $id: $qid\n" if $debug || $ctx->{cla}->{verbose};
		
    print "Processing $id: $count finished, " . ($totalToDo - $count) . " remaining\n" if $count % 10 == 0;

    ##first get the object...
    my $obj = $ctx->{cla}->{table}->new({ $table_pk => $id });
    if (!$obj->retrieveFromDB()) {
      print "ERROR: unable to retrieve $id from database\n";
      next;
    }

    ##need to get any current evidence and mark deleted...
    if ($ctx->{cla}->{deleteEvidence}) {
      $obj->retrieveAllEvidenceFromDB();
      foreach my $e ($obj->getAllEvidence()) {
        next unless $e->getAttributeName() eq $ctx->{cla}->{attribute};
        print STDERR "debug: marking evidence ".$e->getId()." deleted\n" if $debug || $ctx->{cla}->{verbose};
        $e->markDeleted();
      }
    }

    my $haveName = 0;

    ##if is an assembly then want to set name from contained mRNA if possible.
    if($ctx->{cla}->{use_mrna} && $ctx->{cla}->{table} eq "Assembly" && $obj->getContainsMrna()){
      print STDERR "Contains mRNA...getting name of longest one..\n" if $debug;
      $mRNAStmt->execute($id);
      while (my($ext_db_id,$accession,$gusName,$gusDesc) = $mRNAStmt->fetchrow_array()) {
        print STDERR "debug: mrnaQuery:($ext_db_id,$accession,$gusName,$gusDesc)\n" if $debug || $ctx->{cla}->{verbose};
        $gusDesc =~ s/^\s*(\S.*\S)\s*$/$1/;
        if($gusDesc){
          &updateAssemblyFromRNA($obj,$ext_db_id,$accession,$gusName,$gusDesc);
          $mRNAStmt->finish();
          $haveName = 1;
          print MGI "DT.$id: $ext_db_id|$accession - contained mRNA\n" if $ctx->{cla}->{dots_mgi_file};
          last;
        }
      }
      
    }
    if ($haveName) { next; }

    $protStmt->execute($qid);
    print STDERR "debug: Retrieving protein Similarities\n" if $debug || $ctx->{cla}->{verbose};
    while (my($sim_id,$aaSeqid,$ext_db_id,$accession,$pname,$pdescription,$ident,$ml,$mant,$exp,$aa_length,$sub_start,$sub_end,$num_matches) = $protStmt->fetchrow_array()) {

      print STDERR "debug:   protSim($aaSeqid,$ext_db_id,$accession,$pname,$pdescription,$ident,$ml)\n" if $debug || $ctx->{cla}->{verbose};
      next if $pdescription =~ /warning/i; ##alu sequence...
      if (&updateNameFromSP($obj,$qid,$aaSeqid,$pname,$pdescription,$ident,$ml,$sim_id,$aa_length,$sub_start,$sub_end,$num_matches)) {
        print MGI "DT.$id: $ext_db_id|$accession, pVal=",&getPValue($mant,$exp),", ML=$ml, %=", &roundOff(($ident/$ml)*100),"\n" if $ctx->{cla}->{dots_mgi_file};
        $protStmt->finish();
        $haveName = 1;
        last;
      }
    }
    if ($haveName) { next; }

    $obj->set($ctx->{cla}->{attribute},"No NR protein Similarities") unless $obj->get($ctx->{cla}->{attribute}) eq "No NR protein Similarities";
    $obj->submit() if $obj->hasChangedAttributes();
  }
  
  my ($updateRNA, $manualDescriptions);


  if ($ctx->{cla}->{copy_manual_descriptions}) {
      $manualDescriptions = &manualDescriptions($dbh);
  }

  if ($ctx->{cla}->{update_rna_descriptions}) {
      $updateRNA = &updateRNA($dbh);
  }

  $ctx->{self_inv}->closeQueryHandle();
  ############################################################
  ###  put an informative summary in the results variable
  ############################################################
  my $results = "Processed $count Entries, updated $updateRNA RNA descriptions, copied $manualDescriptions manual descriptions to Assembly";

  print STDERR "\n$results\n";

  return "$results";
}

sub updateAssemblyFromRNA {
  my($assem,$ext_db_id,$gacc,$gname,$gdesc) = @_;
  print STDERR "updateNameFromRNA: ($assem,$gacc,$gname,$gdesc)\n" if $ctx->{cla}->{verbose};
  ##first update name and desc if different from above
##  $assem->set('name',$gname) if $gname ne $assem->get('name');
  my $desc = substr($gdesc,0,255);
  $assem->set('description',$desc) if $desc ne $assem->get('description');
  
  ##should do the evidence stuff regardless of whether is modified this time as don't have at all!
  ##external_db_id could be 78 or 22
  my $fact = ExternalNASequence->new({ 'external_db_id' => $ext_db_id,
                                       'source_id' => $gacc });
  
  ##if have retrieved fact tehn add evidence..
  if ($ctx->{cla}->{addEvidence}) {
    if($fact->retrieveFromDB()){
      $assem->addEvidence($fact,1,"description");
    }
  }
  
  ##now submit the RNA...will only be updated if sets were done...
  $assem->submit() if $assem->hasChangedAttributes();
}

sub updateNameFromSP {
  my($obj,$qid,$aaSeqId,$name,$desc,$ident,$ml,$sim_id,$aa_length,$sub_start,$sub_end,$num_matches) = @_;
  return undef unless $desc;    ##don't want to do this on if no description

  my $pcid = &roundOff(($ident/$ml)*100);
  my $pcCov;
  
  ##if the matchlength > sub_end - sub_start (overlapping..) want to compute region covered..
  if ($ml > $sub_end - $sub_start + 1 || $ctx->{cla}->{addEvidence}) {
    print STDERR "\nEXTRACTING \% coverage from similarity...num_matches = $num_matches\n\n" if $debug;
    my $simF = Similarity->new({'similarity_id' => $sim_id});
    $simF->retrieveFromDB();
    $obj->addEvidence($simF,1,$ctx->{cla}->{attribute}) if $ctx->{cla}->{addEvidence};
    if($ml > $sub_end - $sub_start + 1){
      print STDERR "MatchLength=$ml\n" if $debug;
      $pcCov = &roundOff(($simF->getSubjectLengthCovered()/$aa_length)*100);
    }
  }
  $pcCov = $pcCov ? $pcCov : &roundOff(($ml/$aa_length)*100);
  $pcCov = $pcCov > 100 ? 100 : $pcCov;
  my $de = "$pcid\% identity to $pcCov\% of $desc";
  my $description = substr($de,0,255);
  print STDERR "\nDT.",$obj->getId()," DESC: $description\n" if $debug;
  $obj->set($ctx->{cla}->{attribute},$description) unless $description eq $obj->get($ctx->{cla}->{attribute});
  $obj->submit() if $obj->hasChangedAttributes();
  return 1;
}

sub roundOff {
  my($num) = @_;
  my $floor = int($num);
  if (($num - $floor) < 0.5) {
    return $floor;
  } else {
    return $floor + 1;
  }
}

sub getPValue {
  my($mant,$exp) = @_;
  return $mant . (($exp != -999999 && $exp != 0) ? "e" . $exp : "");
}

sub updateRNA {
    my ($dbh) = @_;
    my $rows = $dbh->prepareAndExecute ("update RNA set description = substr((
      select a.description from rnafeature rf, rnasequence rs, assembly a 
      where rs.rna_id = rna.rna_id and rf.na_feature_id = rs.na_feature_id
      and a.na_sequence_id = rf.na_sequence_id ),0,255)
      where rna_id in (
      select rs1.rna_id from rnasequence rs1, rnafeature rf1, assembly a1
      where a1.taxon_id = $ctx->{cla}->{taxon_id} and a1.na_sequence_id = rf1.na_sequence_id
      and rs1.na_feature_id = rf1.na_feature_id )
      and ( (rna_category_id != 17 or rna_category_id is null )
      or (rna_category_id = 17 and (description like '%identity%' or description like 'No NR%')))"); 
    if ($ctx->{'commit'}) { $dbh->commit;}
    return $rows;
}

sub manualDescriptions {
    my ($dbh) = @_;
    my $rows = $dbh->prepareAndExecute ("update assembly set description = 
     (select r.description from RNA r, rnasequence rs, rnafeature rf
     where rf.na_sequence_id = assembly.na_sequence_id
     and rs.na_feature_id = rf.na_feature_id and rs.rna_id = r.rna_id )
     where taxon_id = $ctx->{cla}->{taxon_id}
     and na_sequence_id in ( select rf1.na_sequence_id from rnafeature rf1, rnasequence rs1, rna r1
     where r1.manually_reviewed = 1 and rs1.rna_id = r1.rna_id
     and rf1.na_feature_id = rs1.na_feature_id)");
    if ($ctx->{'commit'}) { $dbh->commit;}
    return $rows;
}


1;

__END__

