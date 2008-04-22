package DoTS::DotsBuild::Plugin::LoadGeneIdInfo;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::PluginMgr::Plugin;
use GUS::Model::SRes::DbRef;
use FileHandle;

$| = 1;

sub new {
    my ($class) = @_;

    my $self = {};
    bless($self,$class);

    my $purpose = <<PURPOSE;
Insert information from downloaded file into sres.dbref pertaining to NCBI Gene DB entries 
PURPOSE

my $purposeBrief = <<PURPOSE_BRIEF;
Insert information pertaining to NCBI Gene DB entries
PURPOSE_BRIEF

my $notes = <<NOTES;
tab delimited input file with
tax_id,GeneId,Symbol,LocusTag,Synonyms,dbXrefs,chromosome,map location,description,type of gene,symbol from nomenclature authority,full name for nomenclature authority,nomenclature status,otherdesignation 
NOTES

my $tablesAffected = <<TABLES_AFFECTED;
    [ ['SRes::DbRef', '']
    ];
TABLES_AFFECTED

my $tablesDependedOn = <<TABLES_DEPENDED_ON;
  [ ['SRes::DbRef','']];
TABLES_DEPENDED_ON

my $howToRestart = <<RESTART;
RESTART

my $failureCases = <<FAIL_CASES;
FAIL_CASES

my $documentation = { purpose          => $purpose,
		      purposeBrief     => $purposeBrief,
		      notes            => $notes,
		      tablesAffected   => $tablesAffected,
		      tablesDependedOn => $tablesDependedOn,
		      howToRestart     => $howToRestart,
		      failureCases     => $failureCases };





    my $argsDeclaration =
	[
	 integerArg({name  => 'testnumber',
            descr => 'number of iterations for testing',
            reqd  => 0,
            constraintFunc=> undef,
            isList=> 0
            }),
	 integerArg({name  => 'externalDbRel',
            descr => 'sres.externaldatabaserelease.id for NCBI Gene DB',
            reqd  => 1,
            constraintFunc=> undef,
            isList=> 0
            }),
	 stringArg({name => 'infoFile',
            descr => 'file downloaded from NCBI Gene containing additional information for rows in DbRef',
            reqd  => 1,
            constraintFunc=> undef,
            isList=> 0
            }),
	 integerArg({name  => 'ncbiTaxId',
            descr => 'ncbi tax_id, e.g. mouse-10090, human-9606',
            reqd  => 1,
	    constraintFunc=> undef,
		     isList=> 0
	    })
	];

    $self->initialize({requiredDbVersion => 3.5,
		       cvsRevision => '$Revision$', # cvs fills this in!
		       name => ref($self),
		       argsDeclaration   => $argsDeclaration,
		       documentation     => $documentation});


  return $self;
}


sub run {

  my $self  = shift;

  my $msg = $self->updateDbRef();

  return $msg;
}


sub updateDbRef {
  my ($self) = @_;

  my $infoFile = $self->getArg('infoFile');

  my $ncbiTaxId = $self->getArg('ncbiTaxId');

  my $externalDbRel = $self->getArg('externalDbRel');

  my $testnumber = $self->getArg('testnumber') if $self->getArg('testnumber');

  open (FILE, $infoFile);

  my $num = 0;

  while (<FILE>) {
    my @line = split(/\t/,$_);

    next if $line[0] != $ncbiTaxId;

    my $newDbRef = GUS::Model::SRes::DbRef->new({'primary_identifier'=>$line[1],'external_database_release_id'=>$externalDbRel});

    $newDbRef->retrieveFromDB;

    $line[6] =  substr($line[6],0,2);

    if ($line[6] ne $newDbRef->get('chromosome')) {

      $newDbRef->setChromosome($line[6]);
    }
    if ($line[2] ne $newDbRef->get('gene_symbol')) {

      $newDbRef->setGeneSymbol($line[2]);
    }
    if ($line[8] ne $newDbRef->get('remark')) {

      $newDbRef->setRemark($line[8]);
    }

    $num += $newDbRef->submit();

    $newDbRef->undefPointerCache();

    exit if ($testnumber && $testnumber >= $num);
  }

  my $msg = "$num DbRef rows updated\n";

  $self->log ($msg);
  close (FILE);

  return $msg;
}



1;
__END__
=pod
=head1 Description
B<LoadGeneIdInfo> - a plug_in that udates rows in DbRef with information from GeneId for C<ga> (GUS application) package.

=head1 Purpose
B<LoadGeneIdInfo> plug_in that updates informtion in DbRef rows.
=cut
