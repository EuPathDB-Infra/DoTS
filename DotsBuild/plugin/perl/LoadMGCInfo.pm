package DoTS::DotsBuild::Plugin::LoadMGCInfo;

@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::Model::SRes::DbRef;
use FileHandle;

$| = 1;

sub new {
    my $class = shift;

    my $self = {};
    bless($self,$class);
    my $usage = 'A package to update sres.DBRef from a tab delimited MGC file.';

    my $easycsp =
	[{o => 'testnumber',
	  t => 'int',
	  h => 'number of iterations for testing',
         },
	 {o => 'externalDbRel',
	  t => 'int',
	  h => 'sres.externaldatabaserelease.external_database_release_id of MGC',
	 },
	 {o => 'infoFile',
	  t => 'string',
	  h => 'file downloaded from MGC containing additional information for rows in DbRef',
         },
	 ];
    
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

  my $self  = shift;

  $self->log ($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");

  $self->log ("Testing on". $self->getArgs()->{'testnumber'}."insertions\n") if $self->getArgs()->{'testnumber'};

  if (!self->getArgs()->{'infoFile'} || !$self->getArgs()->{'externalDbRel'}) {

    die "--infoFile --ncbiTaxId --externalDbRel must be supplied\n";

  }

  $self->updateDbRef();
}


sub updateDbRef {
  my ($self) = @_;

  my $infoFile = $self->getArgs()->{'infoFile'};

  my $externalDbRel = $self->getArgs()->{'externalDbRel'};

  my $testnumber = $self->getArgs()->{'testnumber'} if $self->getArgs()->{'testnumber'};

  my $fh = new FileHandle($infoFile) || die "Could not open". $self->getArgs()->{'infoFile'};

  my $num = 0;

  while (<$fh>) {
    chomp;

    my @line = split(/\t/,$_);

    my $primaryId;

    if ($line[1] =~ /(MGC:\d+)/) {
      $primaryId = $1;
    }
    else {
      next;
    }
    my $newDbRef = GUS::Model::SRes::DbRef->new({'lowercase_primary_identifier'=>lc($primaryId),'external_database_release_id'=>$externalDbRel});

    $newDbRef->retrieveFromDB;

    if ($primaryId ne $newDbRef->get('primary_identifier')) {

      $newDbRef->setPrimaryIdentifier($primaryId);
    }
    if ($line[4] ne $newDbRef->get('secondary_identifier')) {

      $newDbRef->setSecondaryIdentifier($line[4]);
    }
    if (lc ($line[4]) ne $newDbRef->get('lowercase_secondary_identifier')) {

      $newDbRef->setLowercaseSecondaryIdentifier(lc($line[4]));
    }
    if ($line[0] ne $newDbRef->get('gene_symbol')) {

      $newDbRef->setGeneSymbol($line[0]);
    }
    if ($line[1] ne $newDbRef->get('remark')) {

      $newDbRef->setRemark($line[1]);
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
B<LoadLocusLinkInfo> - a plug_in that udates rows in DbRef with information from MGC for C<ga> (GUS application) package.

=head1 Purpose
B<LoadLocusLinkInfo> plug_in that updates informtion in DbRef rows.

=cut
