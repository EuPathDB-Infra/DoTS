#!@perl@

## gets max (project_id) from project table, makes entry into project table and 
## inserts entries into ProjectLink table  

use strict;

use Getopt::Long;
use lib "$ENV{GUS_HOME}/lib/perl";
use GUS::ObjRelP::DbiDatabase;
use GUS::Common::GusConfig;


my ($gusConfigFile,$verbose,$allgenes_num,$restart,$commit,$taxon,$project_id);
&GetOptions("verbose!"=> \$verbose,
	    "gusConfigFile=s" => \$gusConfigFile,
            "allgenes_num=s" => \$allgenes_num,
	    "project_id=i" => \$project_id,
	    "restart!" => \$restart,
	    "commit!" => \$commit,
	    "taxon=i" => \$taxon);

if(! $allgenes_num || ! $taxon){
	die "usage: makeProjectLink.pl --verbose --commit <commit inserts> --restart <optional, if restart required> --allgenes_num <version num of allgenes> --gusConfigFile [\$GUS_CONFIG_FILE] --\n";
}


print STDERR "Establishing dbi login\n" if $verbose;
my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

my $db = GUS::ObjRelP::DbiDatabase->new($gusconfig->getDbiDsn(),
					$gusconfig->getReadOnlyDatabaseLogin(),
					$gusconfig->getReadOnlyDatabasePassword,
					$verbose,0,1,
					$gusconfig->getCoreSchemaName());


my $dbh = $db->getQueryHandle();

my $project_id = $project_id ? $project_id : &getProject($dbh);

my $idHash = &getHash($dbh,$taxon);

&idsDone($dbh,$idHash,$project_id) if ($restart);

&insertProjectLink($dbh,$idHash,$project_id);



sub getProject {
    
    my ($db) = @_;

    my $sql = "select max(project_id) from core.projectinfo";
    
    my $stmt = $db->prepareAndExecute($sql);
    
    my ($project) = $stmt->fetchrow_array();
    
    my $allgenes = "Allgenes-$allgenes_num";
    
    my $description = "AllGenes release $allgenes_num containing public Assembly sequences";
    
    if (!$restart) {
	$project++;
	my $rows = $db->do("insert into core.ProjectInfo Values ($project,'$allgenes','$description',SYSDATE,1,1,1,1,1,0,12,0,0,0)");
	if ($rows) {
	    print STDERR ("Row inserted into core.ProjectInfo table, project_id = $project\n");
	    $db->commit() if ($commit);
	    
	}
	else {
	    print STDERR ("Row insertion into ProjectInfo table failed\n");
	    exit;
	}
    }
    return $project;
}

sub getHash {
    
    my ($db, $taxon_id) = @_; 

    # subtract away assemblies that are exclusively imclone ESTs.
    my $sql = "select na_sequence_id from dots.assembly where taxon_id = $taxon_id";

    my $stmt = $db->prepareAndExecute($sql);

    my %hash;

    while (my $id = $stmt->fetchrow_array()) {
	$hash{$id} = 1;
    }

    return \%hash;
}


sub idsDone {
    
    my ($db,$hash,$project) = @_;

    my $sql = "select id from core.projectlink where project_id = $project and table_id = 56";

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

    my $nextvalSql = $db->getDbPlatform()->getNextValSql("dots.ProjectLink");

    my $sql = "insert into dots.ProjectLink Values ($nextvalSql,$project,56,?,null,SYSDATE,1,1,1,1,1,0,12,0,$project,0)";
    print STDERR ("$sql\n");

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

