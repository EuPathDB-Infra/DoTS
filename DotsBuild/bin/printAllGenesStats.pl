#!/usr/local/bin/perl

##Note that am presenting as the number of "genes" (for now represent as  clusters
## of RNAs)  Simplest to create a tmp table to hold the mapping btwn assemblies
## and genes to facilitate queries...

## Brian Brunk 11/28/00

use strict;
use DBI;
use Getopt::Long;
use Objects::dbiperl_utils::DbiDatabase;

$| = 1;

my ($login,$password,$ornl_ext_db_id,$verbose,$go_cvs_version,$createTmpTableOnly,$useExistingTmpTable);
my $dropTmpTable = 0;  #3drop the tmp table until learn how to keep it..
&GetOptions('login=s' => \$login,'dropTmpTable!' => \$dropTmpTable,'ornl_ext_db_id=i'=>\$ornl_ext_db_id,
            'verbose!' => \$verbose, 'go_cvs_version=s' => \$go_cvs_version,"createTmpTableOnly!" => \$createTmpTableOnly,
            'useExistingTmpTable' => \$useExistingTmpTable,
           );

#die "you must supply --ornl_ext_db_id <external_db_id or NASequences that have been annotated by ORNL> on the command line\n  usage: printAllGenesStats.pl --login [GUSdev] --ornl_ext_db_id 186 --dropTmpTable!" unless $ornl_ext_db_id;

$go_cvs_version = $go_cvs_version ? $go_cvs_version : '2.155';

#put all global vars here...
my($humAssSeqs,$musAssSeqs,$totalHumAssemblies,$totalMusAssemblies);
my $humNonSingletons = 0;
my $musNonSingletons = 0;
my $humDoTSGenes = 0;
my $musDoTSGenes = 0;
my $humORNLContigs = 0;
my $musORNLContigs = 0;
my $humGenScanPreds = 0;
my $musGenScanPreds = 0;
my $humGrailPreds = 0;
my $musGrailPreds = 0;
my $humGenesWithNR = 0;
my $musGenesWithNR = 0;
my $humAssWithNR = 0;
my $musAssWithNR = 0;
my $humDoTSWithNR = 0;
my $musDoTSWithNR = 0;
my $humGrailWithNR = 0;
my $musGrailWithNR = 0;
my $humGenScanWithNR = 0;
my $musGenScanWithNR = 0;
my $humGenesWithCellRoles = 0;
my $musGenesWithCellRoles = 0;
my $humDoTSWithCellRoles = 0;
my $musDoTSWithCellRoles = 0;
my $musAssWithCellRoles = 0;
my $humAssWithCellRoles = 0;
my $humGrailWithCellRoles = 0;
my $musGrailWithCellRoles = 0;
my $humGenScanWithCellRoles = 0;
my $musGenScanWithCellRoles = 0;
my $humGenesWithGOFunction = 0;
my $musGenesWithGOFunction = 0;
my $humDoTSWithGOFunction = 0;
my $musDoTSWithGOFunction = 0;
my $musAssWithGOFunction = 0;
my $humAssWithGOFunction = 0;
my $humGrailWithGOFunction = 0;
my $musGrailWithGOFunction = 0;
my $humGenScanWithGOFunction = 0;
my $musGenScanWithGOFunction = 0;
my %humDoTSRoles;
my %musDoTSRoles;
my $totalDoTShumRoles = 0;
my $totalDoTSmusRoles = 0;
my %humGrailRoles;
my %musGrailRoles;
my %humGenScanRoles;
my %musGenScanRoles;
my $totalGrailhumRoles = 0;
my $totalGenScanhumRoles = 0;
my $totalGrailmusRoles = 0;
my $totalGenScanmusRoles = 0;
my %humDoTSGOFun;
my %musDoTSGOFun;
my $totalDoTShumGOFun = 0;
my $totalDoTSmusGOFun = 0;
my %humGrailGOFun;
my %musGrailGOFun;
my %humGenScanGOFun;
my %musGenScanGOFun;
my $totalGrailhumGOFun = 0;
my $totalGenScanhumGOFun = 0;
my $totalGrailmusGOFun = 0;
my $totalGenScanmusGOFun = 0;
my $humConsAssemblies = 0;
my $musConsAssemblies = 0;
my $humConsDoTS = 0;
my $musConsDoTS = 0;
my $humDoTSConsGenes = 0;
my $musDoTSConsGenes = 0;
my %cellRoleIdMap;
my %goFunctionIdMap;

$login = 'brunkb' unless $login;
$password = 'cyhtn01' unless $password;
if(!$password){
  print STDERR "Please enter the password for login $login\n";
  $password = <STDIN>;
  chomp $password;
}

print STDERR "Establishing dbi login\n" if $verbose;
my $db = new DbiDatabase( undef, $login, $password, $verbose, 0, 1, 'GUSdev' ); ##note that am using GUSdev objects..

my $dbh = $db->getQueryHandle();
my $stmt;

# ##first create the tmp table that will be used to determine "genes"
print STDERR "creating tmp table\n";
##first check to see if it already exists....will drop it at end...
$stmt = $dbh->prepare("select table_name from all_tables where table_name = 'PRINTSTATSTMP'");
$stmt->execute();
my $tmpTableExists = 0;
while(my($tn) = $stmt->fetchrow_array()){
  $tmpTableExists = 1;
}

my $rows = 0;

if(! $useExistingTmpTable){
	##should drop the table first if it exists...
  print STDERR "\nDropping PrintStatsTmp: ",$dbh->do("drop table PrintStatsTmp"),"\n" if $tmpTableExists;
	print STDERR "\nCreating PrintStatsTmp\n";
	$rows = $dbh->do("create table PrintStatsTmp as select r.transcript_unit_id,a.na_sequence_id,a.number_of_contained_sequences, 1 as total_seqs,a.taxon_id,ts.aa_sequence_id from GUSdev.RNA r, GUSdev.RNASequence rs, GUSdev.NAFeature f, GUSdev.Assembly a, GUSdev.TranslatedAAFeature tf, GUSdev.TranslatedAASequence ts where r.rna_id = rs.rna_id and rs.na_feature_id = f.na_feature_id and a.na_sequence_id = f.na_sequence_id and a.taxon_id in (8,14) and tf.na_feature_id = f.na_feature_id and ts.aa_sequence_id = tf.aa_sequence_id ");
  $dbh->commit();
	print STDERR "tmp table created..entered $rows rows\n";
#	print STDERR "Creating stats_tmp table: ".$dbh->do("create global temporary table stats_tmp select transcript_unit_id from PrintStatsTmp group by transcript_unit_id having sum(number_of_contained_sequences) > 1")." non-singleton genes\n";
#	print STDERR "  Upating total_seqs: ".$dbh->do("update PrintStatsTmp set total_seqs = 2 where transcript_unit_id in (select transcript_unit_id from stats_tmp)")." rows\n";
	print STDERR "  Upating total_seqs: ".$dbh->do("update PrintStatsTmp set total_seqs = 2 where transcript_unit_id in (select transcript_unit_id from  PrintStatsTmp group by transcript_unit_id having sum(number_of_contained_sequences) > 1)")." rows\n";
  $dbh->commit();
  ##want to updat total_seqs to two where sequence is singleton but mRNA
  print STDERR "  Upating total_seqs: ".$dbh->do("update printstatstmp set total_seqs = 2 where na_sequence_id in ( select na_sequence_id from gusdev.assembly where contains_mrna = 1 ) and total_seqs = 1");
  $dbh->commit();
	print STDERR "  Updating number_of_contained_sequences: ".$dbh->do("update PrintStatsTmp set number_of_contained_sequences = 2 where number_of_contained_sequences > 1")." rows\n";
  $dbh->commit();
  ##create indexes to facillitate queries
  $dbh->do("create index tran_id_ind on printstatstmp (transcript_unit_id)");
  $dbh->do("create index na_id_ind on printstatstmp (na_sequence_id)");
  $dbh->do("create index aa_id_ind on printstatstmp (aa_sequence_id)");
  $dbh->commit();
  
}else{
	print STDERR "using existing PrintStatsTmp table and data\n";
}

if( $createTmpTableOnly){
  print STDERR "Creating TmpTable only....exiting\n";
  exit;
} 

##first get the number of Assembled Sequences
$stmt = $dbh->prepare("select a.taxon_id,sum(a.number_of_contained_sequences),count(*) from GUSdev.Assembly a group by a.taxon_id");
print STDERR "Retrieving number of AssemblySequences\n" ;
$stmt->execute();

while(my($taxon_id,$num,$count) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humAssSeqs = $num;
		$totalHumAssemblies = $count;
	}elsif($taxon_id == 14){
		$musAssSeqs = $num;
		$totalMusAssemblies = $count;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
	}
}
$stmt = $dbh->prepareAndExecute("select taxon_id,count(*) from GUSdev.Assembly where number_of_contained_sequences > 1 group by taxon_id");
while(my($taxon_id,$nonSing) = $stmt->fetchrow_array()){
  if($taxon_id == 8){
		$humNonSingletons = $nonSing;
	}elsif($taxon_id == 14){
		$musNonSingletons = $nonSing;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $nonSing nonSingletons\n";
	}
}
print STDERR "  Human: $humAssSeqs, $totalHumAssemblies, $humNonSingletons\n  Mouse: $musAssSeqs, $totalMusAssemblies, $musNonSingletons\n\n" ;

##now the number of DoTS genes....
$stmt = $dbh->prepare("select taxon_id,count(distinct transcript_unit_id) from PrintStatsTmp where total_seqs > 1 group by taxon_id");
print STDERR "Retrieving number of DoTS genes\n" ;
$stmt->execute();

while(my($taxon_id,$num) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humDoTSGenes = $num;
	}elsif($taxon_id == 14){
		$musDoTSGenes = $num;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
	}
}
print STDERR "  Human: $humDoTSGenes\n  Mouse: $musDoTSGenes\n\n" ;

##now want to get the number of consistent DoTS alignments on genomic sequence...only for human...
$stmt = $dbh->prepare("select taxon_id,count(distinct transcript_unit_id)
 from PrintStatsTmp t, GUSdev.ConsistentAlignment c
 where t.total_seqs > 1
 and c.transcript_na_sequence_id = t.na_sequence_id
 and c.is_consistent = 1
 group by taxon_id");
print STDERR "Retrieving number of Consistent DoTS genes\n";
$stmt->execute();

while(my($taxon_id,$num) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humDoTSConsGenes = $num;
	}elsif($taxon_id == 14){
		$musDoTSConsGenes = $num;
	}else{
#		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num AssemblySequences\n";
	}
}
print STDERR "  Human: $humDoTSConsGenes\n  Mouse: $musDoTSConsGenes\n\n" ;

##now the number of consistent DoTS and Assemblies...
$stmt = $dbh->prepare("select taxon_id,number_of_contained_sequences,count(*)
 from PrintStatsTmp t, GUSdev.ConsistentAlignment c
 where c.transcript_na_sequence_id = t.na_sequence_id
 and c.is_consistent = 1
 group by taxon_id,number_of_contained_sequences");
print STDERR "Retrieving number of Consistent DoTS sequences\n" ;
$stmt->execute();

while(my($taxon_id,$num_seqs,$count) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		if($num_seqs == 2){ $humConsAssemblies = $count;}
		$humConsDoTS += $count; ##total number...
	}elsif($taxon_id == 14){
		if($num_seqs == 2){ $musConsAssemblies = $count;}
		$musConsDoTS += $count; ##total number...
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $count consistent\n";
	}
}
print STDERR "  Human: $humConsDoTS, $humConsAssemblies \n  Mouse: $musConsDoTS, $musConsAssemblies\n\n" ;

##now the number of ORNL sequences assessed
# print STDERR "Retrieving number of ORNL contigs annotated with GeneFeatures\n" ;
# my $sql = "select s.taxon_id,f.external_db_id,count(*),count(distinct s.na_sequence_id) from NASequenceImp s, GeneFeature f where s.external_db_id = $ornl_ext_db_id and f.na_sequence_id = s.na_sequence_id group by s.taxon_id,f.external_db_id";
# #print STDERR "$sql\n";
# $stmt = $dbh->prepare($sql);
# $stmt->execute();

# while(my($taxon_id,$ext_db_id,$totalPreds,$totalSeqs) = $stmt->fetchrow_array()){
# #	print STDERR "($taxon_id,$ext_db_id,$totalPreds,$totalSeqs)\n";
# 	if($taxon_id == 8){
# 		$humORNLContigs = $totalSeqs < $humORNLContigs ? $humORNLContigs : $totalSeqs;
# 		if($ext_db_id == 183){$humGenScanPreds = $totalPreds;}
# 		elsif($ext_db_id == 184){$humGrailPreds = $totalPreds;}
# 		else{print STDERR "ERROR: retrieving GenePrediction numbers: unknown external_db_id '$ext_db_id'\n";}
# 	}elsif($taxon_id == 14){
# 		$musORNLContigs = $totalSeqs < $musORNLContigs ? $musORNLContigs : $totalSeqs;
# 		if($ext_db_id == 183){$musGenScanPreds = $totalPreds;}
# 		elsif($ext_db_id == 184){$musGrailPreds = $totalPreds;}
# 		else{print STDERR "ERROR: retrieving GenePrediction numbers: unknown external_db_id '$ext_db_id'\n";}
# 	}else{
# 		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totalSeqs ORNL contigs\n";
# 	}
# }
# print STDERR "  Human: $humORNLContigs, GenScan: $humGenScanPreds, Grail: $humGrailPreds\n  Mouse: $musORNLContigs, GenScan: $musGenScanPreds, Grail: $musGrailPreds\n\n" ;

## Now determine the "known" Genes ...
print STDERR "Determining number of DoTS genes with NR neigbors:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,count(distinct pst.transcript_unit_id) from PrintStatsTmp pst, GUSdev.Similarity s where s.query_table_id = 56 and s.query_id = pst.na_sequence_id and s.subject_table_id = 83 and pst.total_seqs > 1 group by pst.taxon_id");
$stmt->execute();

while(my($taxon_id,$totGenes) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humGenesWithNR = $totGenes;
	}elsif($taxon_id == 14){
		$musGenesWithNR = $totGenes;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totGenes Assembly neighbors\n";
	}
}
print STDERR "  Human: $humGenesWithNR\n  Mouse: $musGenesWithNR\n\n" ;

##now assemblies ...
print STDERR "Determining number of DoTS assemblies with NR neigbors:\n" ;
$stmt = $dbh->prepare("select a.taxon_id,a.number_of_contained_sequences,count(distinct a.na_sequence_id) from PrintStatsTmp a, GUSdev.Similarity s where s.query_table_id = 56 and s.query_id = a.na_sequence_id and s.subject_table_id = 83 group by a.taxon_id,a.number_of_contained_sequences");
$stmt->execute();

while(my($taxon_id,$num_seqs,$num) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humAssWithNR = $num unless $num_seqs == 1;
		$humDoTSWithNR += $num;
	}elsif($taxon_id == 14){
		$musAssWithNR = $num unless $num_seqs == 1;
		$musDoTSWithNR += $num;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num Assembly neighbors\n";
	}
}
print STDERR "  Human: Assem: $humAssWithNR, DoTS: $humDoTSWithNR\n  Mouse: Assem: $musAssWithNR, DoTS: $musDoTSWithNR\n\n" ;


##now ORNL things with know neighbors...
# print STDERR "Retrieving number of ORNL predictions with NR similarity\n" ;
# #print STDERR "NOTE: ASSUMING ALL SEQUENCES ARE HUMAN UNTIL CHANGED!!!\n";
# #$stmt = $dbh->prepare("select ts.external_db_id,count(distinct ts.aa_sequence_id) from TranslatedAASequence ts, Similarity s where s.query_table_id = 337 and s.query_id = ts.aa_sequence_id and s.subject_table_id = 83 group by ts.external_db_id");
# $stmt = $dbh->prepare("select nas.taxon_id,aas.external_db_id,count(distinct aas.aa_sequence_id) from Similarity s, TranslatedAASequence aas, TranslatedAAFeature tf, RNAFeature rf, GeneFeature gf, NASequenceImp nas where nas.external_db_id = $ornl_ext_db_id and gf.na_sequence_id = nas.na_sequence_id and rf.parent_id = gf.na_feature_id and tf.na_feature_id = rf.na_feature_id and aas.aa_sequence_id = tf.aa_sequence_id and s.query_id = aas.aa_sequence_id and s.query_table_id = 337 and s.subject_table_id = 83 group by nas.taxon_id,aas.external_db_id");
# $stmt->execute();

# while(my($taxon_id,$ext_db_id,$totalSeqs) = $stmt->fetchrow_array()){
# 	if($ext_db_id == 183){
# 		if($taxon_id == 8){$humGenScanWithNR = $totalSeqs;}
# 		elsif($taxon_id == 14){$musGenScanWithNR = $totalSeqs;}
# 		else{print STDERR "ERROR: retrieving known GenScan predictions..unknown taxon_id '$taxon_id'\n";}
# 	}elsif($ext_db_id == 184){
# 		if($taxon_id == 8){$humGrailWithNR = $totalSeqs;}
# 		elsif($taxon_id == 14){$musGrailWithNR = $totalSeqs;}
# 		else{print STDERR "ERROR: retrieving known grail predictions..unknown taxon_id '$taxon_id'\n";}
# 	}else{
# 		print STDERR "ERROR: unknown external_db_id '$ext_db_id' has $totalSeqs hits to NR\n";
# 	}
# }
# print STDERR "  Human: GrailVsNR: $humGrailWithNR, GenScanVsNR: $humGenScanWithNR\n  Mouse: GrailVsNR: $musGrailWithNR, GenScanVsNR: $musGenScanWithNR\n";

##now do the cell roles...
##DoTS Genes
##print STDERR "Determining number of DoTS genes with CellRoles:\n" ;
#$stmt = $dbh->prepare("select pst.taxon_id,count(distinct pst.transcript_unit_id) from PrintStatsTmp pst, AASequenceCellRole scr where pst.total_seqs > 1 and scr.aa_sequence_id = pst.aa_sequence_id group by pst.taxon_id");
#$stmt->execute();

#while(my($taxon_id,$totGenes) = $stmt->fetchrow_array()){
#	if($taxon_id == 8){
#		$humGenesWithCellRoles = $totGenes;
#	}elsif($taxon_id == 14){
#		$musGenesWithCellRoles = $totGenes;
#	}else{
#		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totGenes CellRoles\n";
#	}
#}
#print STDERR "  Human: $humGenesWithCellRoles Mouse: $musGenesWithCellRoles\n\n" ;

#DoTS assemblies and singletons
#print STDERR "Determining number of DoTS Assemblies with CellRoles:\n" ;
#$stmt = $dbh->prepare("select pst.taxon_id,pst.number_of_contained_sequences,count(distinct pst.na_sequence_id) from PrintStatsTmp pst, AASequenceCellRole scr where scr.aa_sequence_id = pst.aa_sequence_id group by pst.taxon_id,pst.number_of_contained_sequences");
#$stmt->execute();

#while(my($taxon_id,$num_seqs,$num) = $stmt->fetchrow_array()){
#	if($taxon_id == 8){
#		$humAssWithCellRoles = $num unless $num_seqs == 1;
#		$humDoTSWithCellRoles += $num;
#	}elsif($taxon_id == 14){
#		$musAssWithCellRoles = $num unless $num_seqs == 1;
#		$musDoTSWithCellRoles += $num;
#	}else{
#		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num CellRoles\n";
#	}
#}
#print STDERR "  Human: Assem: $humAssWithCellRoles, DoTS: $humDoTSWithCellRoles, Mouse: Assem: $musAssWithCellRoles, DoTS: $musDoTSWithCellRoles\n\n" ;
#
##ORNL cell roles
#print STDERR "Determining number of ORNL with CellRoles:\n" ;
#$stmt = $dbh->prepare("select nas.taxon_id,aas.external_db_id,count(distinct acr.aa_sequence_id) from AASequenceCellRole acr, TranslatedAASequence aas, TranslatedAAFeature tf, RNAFeature rf, GeneFeature gf, NASequenceImp nas where nas.external_db_id = $ornl_ext_db_id and gf.na_sequence_id = nas.na_sequence_id and rf.parent_id = gf.na_feature_id and tf.na_feature_id = rf.na_feature_id and aas.aa_sequence_id = tf.aa_sequence_id and acr.aa_sequence_id = aas.aa_sequence_id group by nas.taxon_id,aas.external_db_id");
#$stmt->execute();

#while(my($taxon_id,$ext_db_id,$totalSeqs) = $stmt->fetchrow_array()){
#	if($ext_db_id == 183){
#		if($taxon_id == 8){$humGenScanWithCellRoles = $totalSeqs;}
#		elsif($taxon_id == 14){$musGenScanWithCellRoles = $totalSeqs;}
#		else{print STDERR "ERROR: retrieving known GenScan predictions..unknown taxon_id '$taxon_id'\n";}
#	}elsif($ext_db_id == 184){
#		if($taxon_id == 8){$humGrailWithCellRoles = $totalSeqs;}
#		elsif($taxon_id == 14){$musGrailWithCellRoles = $totalSeqs;}
#		else{print STDERR "ERROR: retrieving known grail predictions..unknown taxon_id '$taxon_id'\n";}
#	}else{
#		print STDERR "ERROR: unknown external_db_id '$ext_db_id' has $totalSeqs hits to CellRoles\n";
#	}
#}
#print STDERR "  Human: Assem: $humAssWithCellRoles, DoTS: $humDoTSWithCellRoles, Mouse: Assem: $musAssWithCellRoles, DoTS: $musDoTSWithCellRoles\n\n" ;
##lastly do the GOFunctions
##DoTS
##DoTS Genes
print STDERR "Determining number of DoTS genes with GoFunctions:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,count(distinct pst.transcript_unit_id) from PrintStatsTmp pst, GUSdev.AASequenceGOFunction scr, GUSdev.GOFunction f where pst.total_seqs > 1 and scr.aa_sequence_id = pst.aa_sequence_id and f.go_function_id = scr.go_function_id and f.go_cvs_version = '$go_cvs_version' group by pst.taxon_id");
$stmt->execute();

while(my($taxon_id,$totGenes) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humGenesWithGOFunction = $totGenes;
	}elsif($taxon_id == 14){
		$musGenesWithGOFunction = $totGenes;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $totGenes GOFunction\n";
	}
}
print STDERR "  Human: $humGenesWithGOFunction Mouse: $musGenesWithGOFunction\n\n" ;

#DoTS assemblies and singletons
print STDERR "Determining number of DoTS Assemblies with GoFunctions:\n" ;
$stmt = $dbh->prepare("select pst.taxon_id,pst.number_of_contained_sequences,count(distinct pst.na_sequence_id) from PrintStatsTmp pst, GUSdev.AASequenceGOFunction scr, GUSdev.GOFunction f where scr.aa_sequence_id = pst.aa_sequence_id and f.go_function_id = scr.go_function_id and f.go_cvs_version = '$go_cvs_version' group by pst.taxon_id,pst.number_of_contained_sequences");

$stmt->execute();

while(my($taxon_id,$num_seqs,$num) = $stmt->fetchrow_array()){
	if($taxon_id == 8){
		$humAssWithGOFunction = $num unless $num_seqs == 1;
		$humDoTSWithGOFunction += $num;
	}elsif($taxon_id == 14){
		$musAssWithGOFunction = $num unless $num_seqs == 1;
		$musDoTSWithGOFunction += $num;
	}else{
		print STDERR "ERROR: unknown taxon_id '$taxon_id' has $num GOFunction\n";
	}
}
print STDERR "  Human: Assem: $humAssWithGOFunction, DoTS: $humDoTSWithGOFunction, Mouse: Assem: $musAssWithGOFunction, DoTS: $musDoTSWithGOFunction\n\n" ;

#ORNL GO functions 
# print STDERR "Determining number of ORNL with GOFunctions:\n" ;
# $stmt = $dbh->prepare("select nas.taxon_id,aas.external_db_id,count(distinct acr.aa_sequence_id) from GOFunction f, AASequenceGOFunction acr, TranslatedAASequence aas, TranslatedAAFeature tf, RNAFeature rf, GeneFeature gf, NASequenceImp nas where nas.external_db_id = $ornl_ext_db_id and gf.na_sequence_id = nas.na_sequence_id and rf.parent_id = gf.na_feature_id and tf.na_feature_id = rf.na_feature_id and aas.aa_sequence_id = tf.aa_sequence_id and acr.aa_sequence_id = aas.aa_sequence_id and f.go_function_id = acr.go_function_id and f.go_cvs_version = '$go_cvs_version' group by nas.taxon_id,aas.external_db_id");
# $stmt->execute();

# while(my($taxon_id,$ext_db_id,$totalSeqs) = $stmt->fetchrow_array()){
# 	if($ext_db_id == 183){
# 		if($taxon_id == 8){$humGenScanWithGOFunction = $totalSeqs;}
# 		elsif($taxon_id == 14){$musGenScanWithGOFunction = $totalSeqs;}
# 		else{print STDERR "ERROR: retrieving known GenScan predictions..unknown taxon_id '$taxon_id'\n";}
# 	}elsif($ext_db_id == 184){
# 		if($taxon_id == 8){$humGrailWithGOFunction = $totalSeqs;}
# 		elsif($taxon_id == 14){$musGrailWithGOFunction = $totalSeqs;}
# 		else{print STDERR "ERROR: retrieving known grail predictions..unknown taxon_id '$taxon_id'\n";}
# 	}else{
# 		print STDERR "ERROR: unknown external_db_id '$ext_db_id' has $totalSeqs hits to GOFunction\n";
# 	}
# }
# print STDERR "  Human: Grail: $humGrailWithGOFunction, GenScan: $humGenScanWithGOFunction, Mouse: Grail: $musGrailWithGOFunction, GenScan: $musGenScanWithGOFunction\n\n" ;

##now want to generate the data for the breakdowns of cellroles and perhaps GOFunctions..
##could present as bars in table rather than pie chart...
##Cell roles to use: (14,19,28,32,40,46,53)
#print STDERR "Determining breakdown of cellular roles for DoTS Genes:\n";
#$stmt = $dbh->prepare("select t.taxon_id,c.name,count(distinct t.transcript_unit_id) from AASequenceCellRole aacr, CellRole c, PrintStatsTmp t where t.total_seqs > 1 and t.aa_sequence_id = aacr.aa_sequence_id and c.cell_role_id in (14,19,28,32,40,46,53) and c.cell_role_id = aacr.cell_role_id group by t.taxon_id,c.name order by t.taxon_id,count(*) desc");
#$stmt->execute();
#
#while(my($taxon_id,$name,$count) = $stmt->fetchrow_array()){
#	if($taxon_id == 8){
#		$humDoTSRoles{$name} = $count;
#		$totalDoTShumRoles += $count;
#	}elsif($taxon_id == 14){
##		$musDoTSRoles{$name} = $count;
#		$totalDoTSmusRoles += $count;
#	}else{print STDERR "unknown taxon_id $taxon_id for $name, $count\n";}
#}
#foreach my $a (sort{$humDoTSRoles{$b} <=> $humDoTSRoles{$a}} keys %humDoTSRoles){
#	print STDERR "  $a: human=$humDoTSRoles{$a}, mus=$musDoTSRoles{$a}\n";
#}

#Cerllroles for ONRL
#print STDERR "Determining breakdown of cellular roles for ORNL predictions:\n";
#$stmt = $dbh->prepare("select nas.taxon_id,aas.external_db_id,cr.name,count(distinct acr.aa_sequence_id)
# from AASequenceCellRole acr, TranslatedAASequence aas, CellRole cr,
# TranslatedAAFeature tf, RNAFeature rf, GeneFeature gf, NASequenceImp nas
# where nas.external_db_id = $ornl_ext_db_id
# and gf.na_sequence_id = nas.na_sequence_id
#  and rf.parent_id = gf.na_feature_id
#  and tf.na_feature_id = rf.na_feature_id
#  and aas.aa_sequence_id = tf.aa_sequence_id
#  and acr.aa_sequence_id = aas.aa_sequence_id
#  and cr.cell_role_id = acr.cell_role_id
#  and cr.cell_role_id in (14,19,28,32,40,46,53)
#  group by nas.taxon_id,aas.external_db_id,cr.name");

# $stmt->execute();


# while(my($taxon_id,$ext_db_id,$name,$count) = $stmt->fetchrow_array()){
#  	if($ext_db_id == 183){
#  		if($taxon_id == 8){$humGenScanRoles{$name} = $count; $totalGenScanhumRoles += $count;}
#  		elsif($taxon_id == 14){$musGenScanRoles{$name} = $count; $totalGenScanmusRoles += $count;}
#  		else{print STDERR "ERROR: retrieving known GenScan predictions..unknown taxon_id '$taxon_id'\n";}
#  	}elsif($ext_db_id == 184){
#  		if($taxon_id == 8){$humGrailRoles{$name} = $count; $totalGrailhumRoles += $count;}
#  		elsif($taxon_id == 14){$musGrailRoles{$name} = $count; $totalGrailmusRoles += $count;}
#  		else{print STDERR "ERROR: retrieving known grail predictions..unknown taxon_id '$taxon_id'\n";}
#  	}else{
#  		print STDERR "ERROR: unknown external_db_id '$ext_db_id' has $count hits to GOFunction\n";
#  	}
# }

# foreach my $a (sort{$humGrailRoles{$b} <=> $humGrailRoles{$a}}keys%humGrailRoles){
# 	print STDERR "  $a: humGrail=$humGrailRoles{$a}, humGS=$humGenScanRoles{$a}, musGrail=$musGrailRoles{$a}, musGS=$musGenScanRoles{$a}\n\n";
# }

# ##now lets generate the GOFunction dist...
##first DoTS Genes..
print STDERR "Determining breakdown of GOFunctions DoTS Genes:\n";
$stmt = $dbh->prepare("select pst.taxon_id,f.name,count(distinct pst.transcript_unit_id)
      from PrintStatsTmp pst, GUSdev.AASequenceGOFunction scr, GUSdev.GOFunction f
      where pst.total_seqs > 1
      and scr.aa_sequence_id = pst.aa_sequence_id
      and f.go_function_id = scr.go_function_id
      and f.minimum_level = 1
  and f.go_cvs_version = '$go_cvs_version'
      group by pst.taxon_id,f.name");
$stmt->execute();

while(my($taxon_id,$name,$count) = $stmt->fetchrow_array()){
	next if $name eq 'obsolete' || $name =~ /unknown/;
	if($taxon_id == 8){
		$humDoTSGOFun{$name} = $count;
		$totalDoTShumGOFun += $count;
	}elsif($taxon_id == 14){
		$musDoTSGOFun{$name} = $count;
		$totalDoTSmusGOFun += $count;
	}else{print STDERR "unknown taxon_id $taxon_id for $name, $count\n";}
}
foreach my $a (keys%humDoTSGOFun){
	print STDERR "  $a: human=$humDoTSGOFun{$a}, mus=$musDoTSGOFun{$a}\n";
}

##now ORNL preds...
# print STDERR "Determining breakdown of GOFunctions for ORNL predictions:\n";
# $stmt = $dbh->prepare("select nas.taxon_id,aas.external_db_id,f.name,count(distinct agf.aa_sequence_id)
#   from AASequenceGOFunction agf, TranslatedAASequence aas, GOFunction f,
#   TranslatedAAFeature tf, RNAFeature rf, GeneFeature gf, NASequenceImp nas
#   where nas.external_db_id = $ornl_ext_db_id
#   and gf.na_sequence_id = nas.na_sequence_id
#   and rf.parent_id = gf.na_feature_id
#   and tf.na_feature_id = rf.na_feature_id
#   and aas.aa_sequence_id = tf.aa_sequence_id
#   and agf.aa_sequence_id = aas.aa_sequence_id
#   and f.go_function_id = agf.go_function_id
#  and f.minimum_level = 1
#   and f.go_cvs_version = '$go_cvs_version'
#    group by nas.taxon_id,aas.external_db_id,f.name");

# $stmt->execute();
# while(my($taxon_id,$ext_db_id,$name,$count) = $stmt->fetchrow_array()){
# 	next if $name eq 'obsolete' || $name =~ /unknown/;
#  	if($ext_db_id == 183){
#  		if($taxon_id == 8){$humGenScanGOFun{$name} = $count; $totalGenScanhumGOFun += $count;}
#  		elsif($taxon_id == 14){$musGenScanGOFun{$name} = $count; $totalGenScanmusGOFun += $count;}
#  		else{print STDERR "ERROR: retrieving known GenScan predictions..unknown taxon_id '$taxon_id'\n";}
#  	}elsif($ext_db_id == 184){
#  		if($taxon_id == 8){$humGrailGOFun{$name} = $count; $totalGrailhumGOFun += $count;}
#  		elsif($taxon_id == 14){$musGrailGOFun{$name} = $count; $totalGrailmusGOFun += $count;}
#  		else{print STDERR "ERROR: retrieving known grail predictions..unknown taxon_id '$taxon_id'\n";}
#  	}else{
#  		print STDERR "ERROR: unknown external_db_id '$ext_db_id' has $count hits to GOFunction\n";
#  	}
# }

# foreach my $a (sort{$humGrailGOFun{$b} <=> $humGrailGOFun{$a}}keys%humGrailGOFun){
# 	print STDERR "  $a: humGrail=$humGrailGOFun{$a}, humGS=$humGenScanGOFun{$a}, musGrail=$musGrailGOFun{$a}, musGS=$musGenScanGOFun{$a}\n";
# }

##get the current date...
#Wed Nov 29 10:00:01 EST 2000
my @date = split(/\s+/,`date`);
my $date = "$date[1] $date[2], $date[-1]";

##now need to get the hashes mapping cell roles and gofunctions to ids....
# $stmt = $dbh->prepare("select name,cell_role_id from CellRole where source = 'TIGR'");
# $stmt->execute();
# while (my($name,$id) = $stmt->fetchrow_array()){
# 	$cellRoleIdMap{$name} = $id;
# }

my $goQuery = "select name,go_function_id from GUSdev.GOFunction where minimum_level = 1 and go_cvs_version = '$go_cvs_version'";

print STDERR "\n$goQuery\n\n" if $verbose;

$stmt = $dbh->prepare($goQuery);
$stmt->execute();
while (my($name,$id) = $stmt->fetchrow_array()){
	$goFunctionIdMap{$name} = $id;
}

##print out the statistics...
&printStats();

if($dropTmpTable){
  $dbh->do("drop table PrintStatsTmp"); 
  $dbh->commit();
}

$stmt->finish();
$dbh->disconnect();

sub printStats {

print <<endOfStats;
<HTML>
<TABLE BGCOLOR="#ffffff" WIDTH="100%">
<TR><TD>
<FONT SIZE="+1"><I><B>1. Total input sequences</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="150" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff" WIDTH="100%">
<TR>
  <TH ROWSPAN="2">&nbsp;</TH>
  <TH BGCOLOR="#cc0000">
     <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS</FONT></TH>
</TR>
<TR>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">(ESTs and mRNAs)</FONT></TH>
</TR>
<TR><TD ALIGN="center">Human</TD><TD ALIGN="center">$humAssSeqs</TD></TR>
<TR><TD ALIGN="center">Mouse</TD><TD ALIGN="center">$musAssSeqs</TD></TR>
</TABLE>

</TD></TR></TABLE>
</TD></TR>
<TR><TD>

<BR CLEAR="both"><BR><BR>

<FONT SIZE="+1"><I><B>2. Predicted and "known" genes</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="200" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff">
<TR>
  <TH ROWSPAN="2">&nbsp;</TH>
  <TH BGCOLOR="#cc0000" COLSPAN="3">
    <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS</FONT></TH>
</TR>
<TR>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">Total</FONT></TH>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">Assemblies</FONT></TH>
  <TH><FONT FACE="helvetica,sans-serif" SIZE="+1">"Genes"</FONT></TH>

<TR>
<TD ALIGN="left">Human</TD><TD ALIGN="center">$totalHumAssemblies</TD><TD ALIGN="center">$humNonSingletons</TD><TD ALIGN="center">$humDoTSGenes</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  Consistent alignment</TD><TD ALIGN="center">$humConsDoTS</TD><TD ALIGN="center">$humConsAssemblies</TD><TD ALIGN="center">$humDoTSConsGenes</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  "Known"</TD><TD ALIGN="center">$humDoTSWithNR</TD><TD ALIGN="center">$humAssWithNR</TD><TD ALIGN="center">$humGenesWithNR</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  GO function assigned</TD><TD ALIGN="center">$humDoTSWithGOFunction</TD><TD ALIGN="center">$humAssWithGOFunction</TD><TD ALIGN="center">$humGenesWithGOFunction</TD>
</TR>

<TR>
<TD ALIGN="left">Mouse</TD><TD ALIGN="center">$totalMusAssemblies</TD><TD ALIGN="center">$musNonSingletons</TD><TD ALIGN="center">$musDoTSGenes</TD>
</TR>
endOfStats

if($musConsDoTS){  ##have consistent DoTS mus alignments..
	print <<endMusCons
<TR>
<TD ALIGN="left">&nbsp;&nbsp;  Consistent alignment</TD><TD ALIGN="center">$musConsDoTS</TD><TD ALIGN="center">$musConsAssemblies</TD><TD ALIGN="center">$musDoTSConsGenes</TD>
</TR>
endMusCons
}

print <<endOfStats;

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  "Known"</TD><TD ALIGN="center">$musDoTSWithNR</TD><TD ALIGN="center">$musAssWithNR</TD><TD ALIGN="center">$musGenesWithNR</TD>
</TR>

<TR>
<TD ALIGN="left">&nbsp;&nbsp;  GO function assigned</TD><TD ALIGN="center">$musDoTSWithGOFunction</TD><TD ALIGN="center">$musAssWithGOFunction</TD><TD ALIGN="center">$musGenesWithGOFunction</TD>
</TR>

</TD></TR></TABLE>
</TABLE>

</TD></TR>
<TR><TD>
<BR CLEAR="both"><BR>

<DIV ALIGN="left">
<B>DOTS Total:</B> all EST assemblies, including singletons<BR>
<B>DOTS Assemblies:</B> EST assemblies with >1 input sequence (i.e., excluding singletons)<BR>
<B>DOTS Genes:</B> Clusters (by similarity) of DoTS assemblies with >1 input sequence (i.e., excluding singletons)<BR>
<B>Consistent alignment:</B> DoTS sequences that align (SIM4) to UCSC Golden Path sequences in a manner consistent with being a processed transcript<BR>
<B>"Known" genes:</B> homology to a GenBank NR protein with BLAST p-value < 1E-5<BR>
<BR><BR>
</TD></TR>
<TR><TD>

endOfStats


##now make a table for the GOFunctions...
print <<gofuns;
<FONT SIZE="+1"><I><B>3. GO function distributions</B></I></FONT><IMG SRC="/images/dblue1.gif" WIDTH="350" HEIGHT="2"><BR><BR>
</TD></TR>
<TR><TD>
<TABLE CELLPADDING="1" CELLSPACING="0" BORDER="0" BGCOLOR="#000000" ALIGN="left">

<TR><TD>
<TABLE CELLPADDING="3" CELLSPACING="1" BORDER="0" BGCOLOR="#ffffff">
<TR>
  <TH ROWSPAN="3">&nbsp;</TH>
  <TH BGCOLOR="#cc0000" COLSPAN="4">
    <FONT FACE="helvetica,sans-serif" SIZE="+1" COLOR="#ffffff">DOTS "Genes"</FONT></TH>
</TR>
<TR>
  <TH COLSPAN="2"><FONT FACE="helvetica,sans-serif" SIZE="+1">Human</FONT></TH>
  <TH COLSPAN="2"><FONT FACE="helvetica,sans-serif" SIZE="+1">Mouse</FONT></TH>
</TR>
<TR>
<TD ALIGN="center">Total</TD><TD ALIGN="center">Percentage</TD>
<TD ALIGN="center">Total</TD><TD ALIGN="center">Percentage</TD>
</TR>
gofuns

foreach my $a (sort{$humDoTSGOFun{$b} <=> $humDoTSGOFun{$a}}keys%humDoTSGOFun){
	&printDistRow($a,$humDoTSGOFun{$a},$totalDoTShumGOFun,$musDoTSGOFun{$a},
								$totalDoTSmusGOFun,
								{'query' => 'goFunction','pm0' => $goFunctionIdMap{$a},
								 'method' => 'dots','organism' => 'both'});
}

print <<bottomOfStats;
</TD></TR></TABLE>
</TABLE>
</TD></TR>
<TR><TD>
<BR CLEAR="both">
<BR>
<FONT SIZE="+1" FACE="helvetica,sans-serif"><I><B>Statistics last updated $date</B></I></FONT>
</TD></TR>
</TABLE>

</HTML>

bottomOfStats

}

##note:  $link is hashref and must have the keys (query,method,organism,pm0 (value is id to be queried for))
sub printDistRow {
	my($name,$left,$lTotal,$right,$rTotal,$link) = @_;
	my $linkName = $name;
	if($link && scalar(keys%$link) == 4){
		$linkName = '<A href="http://www.allgenes.org/gc/servlet?page=query&rowsPerPage=50&run=submit';
		foreach my $key (keys %$link){
			if($key eq 'method'){
				$linkName .= '&pm1=' . ($link->{$key} =~ /dots/i ? '%3D+4' : '%3D+5');
			}elsif($key =~ /org/){
				$linkName .= '&pm2=' . ($link->{$key} =~ /both/i ? '%3E+7' : ($link->{$key} =~ /^h/i ? '%3D+8' : ' %3D+14'));
			}else{
				$linkName .= "\&$key=$link->{$key}";
			}
		}
		$linkName .= "\">$name</A>"
	}
	my ($dl,$ll) = &getSizes($left,$lTotal);
	my $ldkblue = "<IMG SRC=\"http://www.allgenes.org/images/blue1.gif\" HEIGHT=\"12\" WIDTH=\"$dl\">";
	my $lltblue = "<IMG SRC=\"http://www.allgenes.org/images/lblue1.gif\" HEIGHT=\"12\" WIDTH=\"$ll\">";
	my ($dr,$lr) = &getSizes($right,$rTotal);
	my $rdkblue = "<IMG SRC=\"http://www.allgenes.org/images/blue1.gif\" HEIGHT=\"12\" WIDTH=\"$dr\">";
	my $rltblue = "<IMG SRC=\"http://www.allgenes.org/images/lblue1.gif\" HEIGHT=\"12\" WIDTH=\"$lr\">";
	print "<TR><TD ALIGN=\"left\">$linkName</TD><TD ALIGN=\"center\">$left</TD><TD>$ldkblue$lltblue\&nbsp;", &getPercent($left,$lTotal),"\%</TD>\n";
	print "<TD ALIGN=\"center\">$right</TD><TD>$rdkblue$rltblue\&nbsp;", &getPercent($right,$rTotal),"\%</TD></TR>\n";
}

sub getPercent {
	my($a,$t) = @_;
	$t = $t ? $t : 1;
	my $p = $a/$t*100;
	$p =~ s/^(\d+\.\d).*/$1/;
	return $p;
}

sub getSizes {
	my($a,$t) = @_;
	$t = $t ? $t : 1;
	my $dark = int($a/$t*200) == 0 ? 1 : int($a/$t*200);
	return ($dark,200 - $dark);
}

1;

