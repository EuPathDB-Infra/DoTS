package DoTS::DotsBuild::Plugin::MakeRNAClustersForAssemblies;

@ISA = qw(GUS::PluginMgr::Plugin);
use strict;

use GUS::Model::DoTS::MergeSplit;
use GUS::Model::DoTS::Assembly;
use GUS::Model::DoTS::RNA;
use GUS::Model::DoTS::RNAInstance;
use GUS::Model::DoTS::Gene;
use GUS::Model::DoTS::Protein;
use GUS::Model::DoTS::RNARNACategory;

sub new {
  my ($class) = @_;
    
  my $self = {};
  bless($self,$class);
    
  my $usage = 'Creates RNA and Gene entries given a .cluster file...';
    
  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'sort_desc',
      t => 'boolean',
      h => 'indicates that user has sorted cluster file by number of sequences DESCENDING!!',
     },
     
     {o => 'clusterfile',
      t => 'string',
      h => 'name of cluster file for input',
     },
     {o => 'logfile',
      t => 'string',
      h => 'name of log file...used for  restarting',
      d => 'makeClusters.log',
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


my $debug = 0;
$| = 1;
my $ms_id =  0;

sub run {
  my $self   = shift;
  my $args = $self->getArgs;

    
  $self->log ($args->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n");
  $self->log ("Testing on $args->{'testnumber'}\n") if $args->{'testnumber'};
    
  die "\nYou must sort the cluster file by descending number of sequences and include --sort_desc on command line\n\n" unless $args->{sort_desc};
    
  open(F,"$args->{'clusterfile'}") || die "clusterfile $args->{'clusterfile'} not found\n";
    
  my $algoInvo = $self->getAlgInvocation;

  ##don't delete evidence
  $algoInvo->setGlobalDeleteEvidenceOnDelete(0);
    
  ##get entries already processed
  my %done;
  if (-e $args->{logfile}) {
    open(L, "$args->{logfile}");
    while (<L>) {
      if (/^(Cluster_\S+)/) {
        $done{$1} = 1;
      }
    }
    $self->log ("restarting: completed ",scalar(keys%done)," clusters\n");
    close L;
  }
  open(L, ">>$args->{logfile}");
  select L; $| = 1; select STDOUT;
    
  my $ct = 0;
  while (<F>) {
    if (/^(Cluster_\S+)\s.*:\s\((.*)\)/) {
      $ms_id++;
      my $cluster = $1;
      $ct++;
      next if $done{$1};
      last if $args->{testnumber} && $ct > $args->{testnumber};
      $self->log ("Processing $ct: $1\n") if $ct % 100 == 0;
      my @ids = split(', ',$2); ##array of related assembly.na_sequence_ids
      my %rnas;
      my %genes;
      foreach my $id (@ids) {
        my $ass = GUS::Model::DoTS::Assembly->new({'na_sequence_id' => $id});
        ##first get all RNA's if exist each of which should have a gene...
        my $rna = $ass->getRNA(1);
        if ($rna) {
          $self->log ("assembly.$id has RNA ",$rna->getId(),"\n") if $debug;
          my $g = $rna->getParent('DoTS::Gene',1);
		    
          if ($g) {
            $genes{"$g"} = $g; 
          }
        } else {
          $self->log ("assembly.id does not have RNA creating new one\n") if $debug;
          $rna = $self->createNewRNAAndSequence();
          $rna->getChild('DoTS::RNAInstance')->setParent($ass->getChild('DoTS::RNAFeature',1)); ##connect feature to rnasequence
          #          my $g = $rna->getParent('DoTS::Gene',1);
          #          $genes{"$g"} = $g; 
        }
        $rnas{"$rna"} = $rna;
      }
      ##now do something..
      ##if more than one gene...get one that is manually_reviewed if possible to use and mark rest
      ##deleted...
      my @genes = values %genes;
      my $gene;
	    
      ##NOTE: assuming that every manually reviewed Gene will have an RNA that 
      ##  is a reference sequence (rna_category_id = 17)....if there is a gene
      ##  that has been manually reviewed and does not have reference RNA, that gene
      ##  may be lost with all associated annotation.  The annotated gene will thus always
      ##  stay associated with the RNA.reference sequence if it exists.
	    
      my $otherReference = 0;
      my %genesWithReference;
      my %manRevGenes;
      foreach my $g (@genes) {
        $manRevGenes{$g} = $g if $g->getReviewStatusId == 1;
        $g->removeAllChildPointers(1);
        $self->log ("Retrieving RNAs for gene ",$g->getId(),"\n") if $debug;
        foreach my $rna ($g->getRNAs(1)) {
          $self->log ("setting gene id to null for RNA.",$rna->getId(),"\n") if $debug;
          
          ##if  is  a reference sequence..keep associated with  gene unless is one of rnas working with..
          if ($rna->getChild('DoTS::RNARNACategory')->getRnaCategoryId() == 17) { ##is a reference rna 
            if (exists $rnas{"$rna"}) { #3and in this cluster.. 
              if (!$gene) {	## don't yet have a gene..
                $gene = $g;
              }
            } else {
              $otherReference++;
              $genesWithReference{"$g"} = 1;
              ##remove from manRevGenes as will NOT use this for this gene..
              delete $manRevGenes{$g};
              next;		##this is a reference RNA but not in this cluster....leave with gene...
            } 
          }      
          #set the gene_id to null
          $rna->removeParent($g);
          ##don't do next if is rna I want as may cause update when unnecessary
          next if $rnas{"$rna"};
          $rna->setGeneId('NULL');
          ##add to list of rna's to submit
          $algoInvo->addChild($rna);
          
        }                       ##end of foreach $rna... 
        $algoInvo->addChild($self->createMergeSplit($g,$gene,1)) if "$g" ne "$gene";
      } 
      ##now get a gene if  don't  have one...
      ##want a gene that is manually reviewed if it does not have other RNA with reference..
      if (!$gene && scalar(@genes) > 0) {
        foreach my $v (values%manRevGenes) {
          $gene = $v;
          last;
        }
        if (!$gene) {
          foreach my $g (@genes) {
            next if $genesWithReference{$g};
            $gene = $g;
            last;
          }
        }
      }
	
      if (!$gene) {
        $gene = $self->createNewGene();
      }				##don't have any so create a new one...
	
      ##now foreach  rna in  cluster...add to the gene
      foreach my $r (values %rnas) {
        $gene->addChild($r);
      }
	
      ##now mark extra genes deleted...
      foreach my $g (@genes) {
        if ("$g" eq "$gene") {
          $self->log ("This gene $g eq my gene $gene\n") if $debug;
          next;
        }
        ##may still have reference rna associated...don't want  to  alter gene...
        next if $genesWithReference{$g};
          
        $g->retrieveAllChildrenFromDB(); ##don't do it recursively
        $g->removeChildren($g->getChildren('RNA'));
        #        foreach my $c ($gene->getChildren('RNA')){
        #          $c->removeParent($g);
        #        }
        $g->markDeleted(1);
        $gene->addToSubmitList($g);
      }
      ##and submit...
      $self->submit($algoInvo,$gene);
      $algoInvo->removeAllChildren();
      $algoInvo->undefPointerCache();
      print L "$cluster complete\n";
    }
      
    close F;
      
    ############################################################
    # return status
    # replace word "done" with meaningful return value/summary
    ############################################################
    my $res = "Processed $ct clusters";
    $self->log ("\n$res\n");
    print L "\n$res\n";
    close L;
    return "$res";
  }
}

  ##need to  submit the rna children first before the genes..
sub submit {
  
  my $self   = shift;
  my($algoInvo,$gene) = @_;
  $algoInvo->manageTransaction(undef,'begin');
  ##need to submit the RNAs first before genes..
  foreach my $rna ($algoInvo->getAllChildren()) {
    $self->log ("Submitting RNA....should have null transcript_unit_id...\n") if $debug;
    $rna->submit(undef,1);
  }
  $gene->submit(undef,1);
  #  foreach my $g ($algoInvo->getChildren('DoTS::Gene')){
  #    $g->submit(undef,1);
  #  }
  return $algoInvo->manageTransaction(undef,'commit');
}

sub createNewRNAAndSequence {
  
  my $self   = shift;
  my $rs = GUS::Model::DoTS::RNAInstance->new({'review_status_id' => 0,
					       'is_reference' => 0,
					       'rna_instance_category_id' => 0 });
  
  ##last the RNA
  my $rna = GUS::Model::DoTS::RNA->new({'review_status_id' => 0});
  $rna->addChild($rs);
  
  #    my $gene = GUS::Model::DoTS::Gene->new({'review_status_id' => 0});
  #    $gene->addChild($rna);
  
  ##need to also create a protein entry so can do GOFunction predictions..
  my $prot = GUS::Model::DoTS::Protein->new({ 'review_status_id' => 0});
  $rna->addChild($prot);
  
  return $rna;
}

#creates a new gene with a TU child...I (DP-3/26/03)changed $gene->addChild(GUS::Model::DoTS::RNA->new({'review_status_id' => 0})); was new TU-however, RNAs are created above so this should be commented out?
sub createNewGene {
  
  my $self   = shift;
  
  my $gene = GUS::Model::DoTS::Gene->new({'is_reference' => 0,
					  'review_status_id' => 0 });
  ##add a transcriptUnit
  #$gene->addChild(GUS::Model::DoTS::RNA->new({'review_status_id' => 0}));
  return $gene;
}

sub createMergeSplit {
  
  my $self   = shift;
  my($o,$n,$is_merge) = @_;
  $self->log ("Creating MergeSplit: Old:".$o->getId().", New:".$n->getId().", is_merge:'$is_merge'\n") if $debug;
  my $ms = GUS::Model::DoTS::MergeSplit->new({'old_id' => $o->getId(),
					      'new_id' => $n->getId(),
					      'table_id' => $o->getTableIdFromTableName($o->getClassName()),
					      'merge_split_group_id' => $ms_id,
					      'is_merge' => $is_merge });
  return $ms;
}


1;

__END__

=pod
  =head1 Description
  B<Template> - a template plug-in for C<ga> (GUS application) package.
  
  =head1 Purpose
  B<Template> is a minimal 'plug-in' GUS application.
  
  =cut
  



