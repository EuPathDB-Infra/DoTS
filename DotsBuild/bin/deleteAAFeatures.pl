#!@perl@

##takes in an SQL query and deletes entries one at a time from an AAFeature table as well as corresponding entries in AALocation...
## query must return primary key of feature table..

use strict;
use lib "$ENV{GUS_HOME}/lib/perl";
use Getopt::Long;
use GUS::ObjRelP::DbiDatabase;
use GUS::Model::Core::TableInfo;
use GUS::Common::GusConfig;

$| = 1;

my ($verbose,$idSQL,$table,$batch_size, $gusConfigFile);
&GetOptions("verbose!"=> \$verbose,
	    "idSQL=s" => \$idSQL, 
	    "gusConfigFile=s" => \$gusConfigFile, 
	    "table=s" => \$table,
	    "batch_size=i" => \$batch_size);

die "usage: deleteEntries.pl --idSQL 'sql query returns primary keys of table' --table <tablename in schema.table format> --verbose --batch_size [100]\n" unless $idSQL && $table;

$batch_size = $batch_size ? $batch_size : 1000;

print STDERR "Establishing dbi login\n" if $verbose;
my $gusconfig = GUS::Common::GusConfig->new($gusConfigFile);

my $db = GUS::ObjRelP::DbiDatabase->new($gusconfig->getDbiDsn(),
					$gusconfig->getReadOnlyDatabaseLogin(),
					$gusconfig->getReadOnlyDatabasePassword,
					$verbose,0,1,
					$gusconfig->getCoreSchemaName());

my $dbh = $db->getQueryHandle();

my $o = GUS::Model::Core::TableInfo->new();

my $pk = $o->getTablePKFromTableId($o->getTableIdFromTableName($table));

print STDERR "Query: $idSQL\n" if $verbose;

my $stmt = $dbh->prepareAndExecute($idSQL);
my @ids;
my $ct = 0;
while(my($id) = $stmt->fetchrow_array()){
  $ct++;
  push(@ids,$id);
  print STDERR "Retrieving ids: $ct\n" if $ct % 1000 == 0;
}
print STDERR "deleting ",scalar(@ids)," ids from $table\n";
$table =~ s/\:\:/\./;
$ct = 0;
for(my $i=0;$i<scalar(@ids);$i += $batch_size){
  $dbh->do("delete from dots.aalocation where aa_feature_id in (".join(', ',@ids[$i..($i + $batch_size - 1 < scalar(@ids) ? $i + $batch_size - 1 : scalar(@ids) - 1)]).")");
  $ct += $dbh->do("delete from $table where $pk in (".join(', ',@ids[$i..($i + $batch_size - 1 < scalar(@ids) ? $i + $batch_size - 1 : scalar(@ids) - 1)]).")");
  $dbh->commit();
  print STDERR "$ct deleted ",`date`;
}
$dbh->disconnect();
$db->logout();
  
