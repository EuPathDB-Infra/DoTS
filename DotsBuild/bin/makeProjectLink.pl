#!/usr/bin/perl

## gets max (project_id) from project table, makes entry into project table and 
## inserts entries into ProjectLink table  

use strict;

use Getopt::Long;
use Objects::dbiperl_utils::DbiDatabase;


my ($login,$password,$verbose,$database,$allgenes_num,$restart,$commit,$taxon);
&GetOptions("verbose!"=> \$verbose,
            "login=s" => \$login,
	    "password=s" => \$password,
            "database=s" => \$database,
            "allgenes_num=s" => \$allgenes_num,
	    "restart!" => \$restart,
	    "commit!" => \$commit,
	    "taxon=i" => \$taxon);

if(!$login || !$password || $database || $allgenes_num){
	die "usage: makeProjectLink.pl --verbose --commit <commit inserts> --restart <optional, if restart required> --allgenes_num <version num of allgenes> --login [GUSrw] --database [GUSdev] --password <database password> --\n";
}


$database = $database ? $database : "GUSdev";
$login = $login ? $login : 'gusrw';

print STDERR "Establishing dbi login\n" if $verbose;
my $db = new DbiDatabase( undef, $login, $password, $verbose, 0, 1, $database );

my $dbh = $db->getQueryHandle();

my $project_id = &getProject($dbh);

my $idHash = &getHash($dbh,$taxon);

&idsDone($dbh,$idHash,$project_id) if ($restart);

&insertProjectLink($dbh,$idHash,$project_id);



sub getProject {
    
    my ($db) = @_; 

    my $sql = "select max(project_id) from project";
    
    my $stmt = $db->prepareAndExecute($sql);
    
    my ($project) = $stmt->fetchrow_array();
    
    my $allgenes = "Allgenes-$allgenes_num";
    
    my $description = "AllGenes release $allgenes_num containing public Assembly sequences";
    
    if (!$restart) {
	$project++;
	my $rows = $db->do("insert into Project Values ($project,$allgenes,$description,SYSDATE,1,1,1,1,1,0,12,0,0,0)");
	if ($rows) {
	    print STDERR ("Row inserted into Project table, project_id = $project\n");
	    $db->commit() if ($commit);
	    
	}
	else {
	    print STDERR ("Row insertion into Project table failed\n");
	    exit;
	}
    }
    return \$project;
}

sub getHash {
    
    my ($db, $taxon_id) = @_; 

    my $sql = "select na_sequence_id from assembly where taxon_id in ($taxon_id) 
            minus (select s.assembly_na_sequence_id from assemblysequence s, externalnasequence e 
            where s.na_sequence_id = e.na_sequence_id and e.taxon_id in ($taxon_id) and e.external_db_id = 3892 
            minus (select s.assembly_na_sequence_id from assemblysequence s, externalnasequence e 
            where s.na_sequence_id = e.na_sequence_id and e.taxon_id in ($taxon_id) and e.external_db_id != 3892))";

    my $stmt = $db->prepareAndExecute($sql);

    my %hash;

    while (my $id = $stmt->fetchrow_array()) {
	$hash{$id} = 1;
    }

    return \%hash;
}


sub idsDone {
    
    my ($db,$hash,$project) = @_;

    my $sql = "select id from projectlink where project_id = $project and table_id = 56";

    my $stmt = $db->prepareAndExecute($sql);

    while (my $id = $stmt->fetchrow_array()) {
	if (exists $hash->{$id}) {
	    delete $hash->{$id};
	}
	else {
	    next;
	}
    }
}

sub insertProjectLink {

    my ($db,$hash,$project) = @_;

    my $count = 0;

    my $num = scalar (keys %{$hash});

    print STDERR ("$num assembly ids to be entered into ProjectLink\n");

    my $sql = "insert into ProjectLink (ProjectLink_sq.nextval,$project,56,?,null,SYSDATE,1
,1,1,1,1,0,12,0,$project,0)";

    my $stmt = $db->prepare($sql);

    foreach my $id (keys %{$hash}) {

	$stmt->execute($id);

	$db->commit if ($commit);

	$count++;

	print STDERR ("$count\n") if (($count % 1000) == 0);

    }

    if ($num-$count == 0) {
	print STDERR ("ProjectLink entry completed - $count assembly ids entered\n");
    }
    else {
	print STDERR ("ProjectLink entry NOT completed - only $count out of $num assembly ids entered\n");
    }
}

	

	

	

 
