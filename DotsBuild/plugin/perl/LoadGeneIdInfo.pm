package DoTS::DotsBuild::Plugin::LoadGeneIdInfo;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::Model::SRes::DbRef;
use FileHandle;

$| = 1;

sub new {
    my ($class) = @_;

    my $self = {};
    bless($self,$class);
    my $usage = 'A package to update sres.DBRef from a tab delimited gene_info file.';

    my $easycsp =
	[{o => 'testnumber',
	  t => 'int',
	  h => 'number of iterations for testing',
         },
	 {o => 'externalDbRel',
	  t => 'int',
	  h => 'sres.externaldatabaserelease.external_database_release_id of GeneId',
	 },
	 {o => 'infoFile',
	  t => 'string',
	  h => 'file downloaded from NCBI Gene containing additional information for rows in DbRef',
         },
	 {o => 'ncbiTaxId',
	  t => 'int',
	  h => 'ncbi tax_id, e.g. mouse-10090, human-9606',
         }
	 ];
    
    $self->initialize({requiredDbVersion => {},
		       cvsRevision => '$Revision$', # cvs fills this in!
		       cvsTag => '$Name$', # cvs fills this in!
		       name => ref($self),
		       revisionNotes => 'make consistent with GUS 3.0',
		       easyCspOptions => $easycsp,
		       usage => $usage
		       });
    
    return $self;
}


sub run {

  my $self  = shift;

  $self->log ($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");

  $self->log ("Testing on". $self->getArgs()->{'testnumber'}."insertions\n") if $self->getArgs()->{'testnumber'};

  if (!$self->getArgs()->{'infoFile'} || !$self->getArgs()->{'ncbiTaxId'} || !$self->getArgs()->{'externalDbRel'}) {

    die "--infoFile --ncbiTaxId --externalDbRel must be supplied\n";

  }

  $self->updateDbRef();
}


sub updateDbRef {
  my ($self) = @_;

  my $infoFile = $self->getArgs()->{'infoFile'};

  my $ncbiTaxId = $self->getArgs()->{'ncbiTaxId'};

  my $externalDbRel = $self->getArgs()->{'externalDbRel'};

  my $testnumber = $self->getArgs()->{'testnumber'} if $self->getArgs()->{'testnumber'};

  my $fh = new FileHandle($infoFile) || die "Could not open". $self->getArgs()->{'infoFile'};

  my $num = 0;

  while (<$fh>) {
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

  $self->log ("$num DbRef rows updated\n");
}



1;
__END__
=pod
=head1 Description
B<LoadGeneIdInfo> - a plug_in that udates rows in DbRef with information from GeneId for C<ga> (GUS application) package.

=head1 Purpose
B<LoadGeneIdInfo> plug_in that updates informtion in DbRef rows.
=cut
