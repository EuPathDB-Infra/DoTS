package DoTS::DotsBuild::Plugin::MakeGenesForRNA;



@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use GUS::ObjRelP::DbiDatabase;
use GUS::Model::DoTS::Gene;
use GUS::Model::DoTS::RNA;

sub new {
  my ($class) = @_;
  my $self = {};
  bless($self,$class);
  
  my $usage = 'Makes gene entries for dots.RNA rows with null gene_id';
  
  my $easycsp =
    [
     {o => 'taxon',
      t => 'int',
      h => 'taxon_id',
     }
    ];
  
  $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$',  # cvs fills this in!
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

  my $rnaArray = $self->getRNA();

  $self->addGene($rnaArray);

}

sub getRNA {
  my ($self) = @_;

  my @rnaArray;

  my $dbh = $self->getQueryHandle();

  my $taxon = $self->getArgs()->{'taxon'} || die "must supply --taxon\n"; 

  my $query = "select r.rna_id from dots.rna r, dots.rnainstance i,dots.rnafeature f, dots.assembly a where r.gene_id is null and r.rna_id = i.rna_id and i.na_feature_id = f.na_feature_id and f.na_sequence_id = a.na_sequence_id and a.taxon_id = $taxon";

  $self->log ("$query\n");

  my $stmt = $dbh->prepareAndExecute($query);

  while(my ($rna_id) = $stmt->fetchrow_array()){
    push (@rnaArray, $rna_id);
  }

  my $num = @rnaArray;
  $self->log ("$num dots.rna rows have null gene_id\n");
  return \@rnaArray;

}

sub addGene {
  my ($self,$rnaArray) = @_;
  my $count;
  foreach my $rna_id (@$rnaArray) {
    my $rna = GUS::Model::DoTS::RNA->new({'rna_id' => $rna_id});
    $rna->retrieveFromDB();
    my $gene = $rna->getParent('DoTS::Gene',1) ? $rna->getParent('DoTS::Gene') : $self->makeGene($rna);
    $rna->addToSubmitList($gene);
    $count += $rna->submit();
    $rna->undefPointerCache();
  }
  $self->log ("$count\n");  
}



sub makeGene {
  my ($self,$rna) = @_;
  my $review_status_id = 0;  
  my $gene = GUS::Model::DoTS::Gene->new({'review_status_id'=>$review_status_id});
  $gene->addChild($rna);
  return $gene;
}






1;

