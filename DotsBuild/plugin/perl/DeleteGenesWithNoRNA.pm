package DoTS::DotsBuild::Plugin::DeleteGenesWithNoRNA;



@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::ObjRelP::DbiDatabase;
use GUS::Model::DoTS::Gene;

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);
  
  my $usage = 'Deletes Genes and children that have no RNA';
  
  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'user_id',
      t => 'int',
      h => 'row_user_id of dots.gene rows that will be deleted',
     }
    ];
  
  $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$ $',  # cvs fills this in!
		     cvsTag => '$Name$', # cvs fills this in!
		     name => ref($self),
		     revisionNotes => ' ',
		     easyCspOptions => $easycsp,
		     usage => $usage
		    });
  
  return $self;
}


$| = 1;

sub run {
  my $self   = shift;
  
  $self->log ($self->getArgs()->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");

  $self->log ("Testing on ".$self->getArgs()->{'testnumber'}."\n") if $self->getArgs()->{'testnumber'};

  my $userId = $self->getArgs()->{'user_id'} if $self->getArgs()->{'user_id'};

  my $rnaHash = $self->getRNA();

  my $geneHash;

  if ($self->getArgs()->{'user_id'})  {

    $geneHash = $self->getGene($rnaHash,$userId);
  }

  else {

    $geneHash = $self->getGene($rnaHash);
  }

  my $ret = $self->deleteGene($geneHash);

  return $ret;

}

sub getRNA {
  my ($self) = @_;

  my %rnaHash;

  my $dbh = $self->getQueryHandle();

  my $query = "select gene_id from dots.rna";

  my $stmt = $dbh->prepareAndExecute($query);

  while(my ($id) = $stmt->fetchrow_array()){
    $rnaHash{$id} = 1;
  }

  return \%rnaHash;

}

sub getGene {
  my ($self,$rnaHash,$userId) = @_;

  my %geneHash;

  my $dbh = $self->getQueryHandle();

  my $query = "select gene_id from dots.gene"; 

  $query .= " where row_user_id = $userId" if $userId;

  my $stmt = $dbh->prepareAndExecute($query);

  my $num;

  while(my ($id) = $stmt->fetchrow_array()){
    $geneHash{$id} = 1 if ($rnaHash->{$id} != 1);
    $num = scalar (keys %geneHash);
    last if $self->getArgs()->{'testnumber'} && $num >= $self->getArgs()->{testnumber};
  }

  $self->log ("$num gene_ids to be deleted\n");

  return \%geneHash;

}

sub deleteGene {

  my ($self,$geneHash) = @_;

  my $algoInvo = $self->getAlgInvocation;

  $algoInvo->setGlobalDeleteEvidenceOnDelete(0);

  my $ctDel = 0;

  foreach my $gene_id (keys %$geneHash){
    $self->log ("Deleting Gene $gene_id\n");
    my $gene = GUS::Model::DoTS::Gene->new({'gene_id' => $gene_id});
    $gene->markDeleted(1);
    if($gene->retrieveFromDB()){
      $gene->retrieveAllChildrenFromDB(1);
      $ctDel += $gene->submit();
      $gene->undefPointerCache();
    }else{
      $self->log ("Unable to retrieve entry for Gene $gene_id from DB\n");
    }
  }

  my $ret = "Deleted $ctDel Gene and \n";

  $self->log ("$ret");

  return "$ret";
}




1;

