############################################################
## Change Package name....
############################################################
package MakeRNAClustersForAssemblies;

use strict;
use DBI;

############################################################
# Add any specific objects (Objects::GUSdev::Objectname) here
############################################################
use Objects::GUSdev::MergeSplit;
use Objects::GUSdev::Assembly;
use Objects::GUSdev::RNA;
use Objects::GUSdev::RNASequence;
use Objects::GUSdev::Gene;
use Objects::GUSdev::TranscriptUnit;
use Objects::GUSdev::Protein;


my $Cfg;  ##global configuration object....passed into constructor as second arg


sub new {
	my $Class = shift;
	$Cfg = shift;  ##configuration object...

	return bless {}, $Class;
}

sub Usage {
	my $M   = shift;
	return 'Creates RNA and Gene entries given a .cluster file...';
}

############################################################
# put the options in this method....
############################################################
sub EasyCspOptions {
	my $M   = shift;
	{

#		test_opt1 => {
#									o => 'opt1',  ##variable name in app = $ctx->{cla}->{opt1}
#		              t => 'int',   ## type of this argument ('int','float','boolean','string')
#									h => 'option 1 for test application',  ##help text
#									d => 4,  ##default value
#									l => 1,	 ##is list valued if true
#                ld => ':',  ##list delimiter (default is ",")
#									e => [ qw( 1 2 3 4 ) ], ##allowed values
#								 },

  testnumber        => {
                        o => 'testnumber',
                        t => 'int',
                        h => 'number of iterations for testing',
                       },
  sort_desc         => {
                        o => 'sort_desc',
                        t => 'boolean',
                        h => 'indicates that user has sorted cluster file by number of sequences DESCENDING!!',
                       },

  clusterfile       => {
                        o => 'clusterfile=s',
                        h => 'name of cluster file for input',
                       },
  logfile           => {
                        o => 'logfile=s',
                        h => 'name of log file...used for  restarting',
                        d => 'makeClusters.log',
                       },


	}
}

my $ctx;
my $debug = 0;
$| = 1;
my $ms_id =  0;

sub Run {
  my $M   = shift;
  $ctx = shift;
  
  print $ctx->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $ctx->{'testnumber'}\n" if $ctx->{'testnumber'};
  
  ############################################################
  # Put loop here...remember to undefPointerCache()!
  ############################################################
  
  die "\nYou must sort the cluster file by descending number of sequences and include --sort_desc on command line\n\n" unless $ctx->{cla}->{sort_desc};
  
  open(F,"$ctx->{'clusterfile'}") || die "clusterfile $ctx->{cla}->{'clusterfile'} not found\n";

  my $algoInvo = $ctx->{self_inv};

  ##don't delete evidence
  $algoInvo->setGlobalDeleteEvidenceOnDelete(0);

  ##get entries already processed
  my %done;
  if(-e $ctx->{cla}->{logfile}){
    open(L, "$ctx->{cla}->{logfile}");
    while(<L>){
      if(/^(Cluster_\S+)/){
        $done{$1} = 1;
      }
    }
    print STDERR "restarting: completed ",scalar(keys%done)," clusters\n";
    close L;
  }
  open(L, ">>$ctx->{cla}->{logfile}");
  select L; $| = 1; select STDOUT;
  
  my $ct = 0;
  while(<F>){
    if(/^(Cluster_\S+)\s.*:\s\((.*)\)/){
      $ms_id++;
      my $cluster = $1;
      $ct++;
      next if $done{$1};
      last if $ctx->{cla}->{testnumber} && $ct > $ctx->{cla}->{testnumber};
      print STDERR "Processing $ct: $1\n" if $ct % 100 == 0;
      my @ids = split(', ',$2);  ##array of related assembly.na_sequence_ids
      my %rnas;
      my %genes;
      foreach my $id (@ids){
        my $ass = Assembly->new({'na_sequence_id' => $id});
        ##first get all RNA's if exist each of which should have a gene...
        my $rna = $ass->getRNA(1);
        if($rna){
          print STDERR "assembly.$id has RNA ",$rna->getId(),"\n" if $debug;
          my $tu = $rna->getParent('TranscriptUnit',1);
          my $g;
          if($tu){ $g = $tu->getParent('Gene',1) ;}
          if($g){
            $genes{"$g"} = $g; 
            $rna->removeParent($tu);  ##remove the parent to gene can retrieve all rnas..
          }
        }else{
          print STDERR "assembly.id does not have RNA creating new one\n" if $debug;
          $rna = &createNewRNAAndSequence();
          $rna->getChild('RNASequence')->setParent($ass->getChild('RNAFeature',1));  ##connect feature to rnasequence
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
      foreach my $g (@genes){
        $manRevGenes{$g} = $g if $g->getManuallyReviewed();
        $g->removeAllChildPointers(1);
        print STDERR "Retrieving RNAs for gene ",$g->getId(),"\n" if $debug;
        foreach my $rna ($g->getRNAs(1)){
          print STDERR "setting transcript unit id to null for RNA.",$rna->getId(),"\n" if $debug;

          ##if  is  a reference sequence..keep associated with  gene unless is one of rnas working with..
          if($rna->getRnaCategoryId() == 17){##is a reference rna 
            if(exists $rnas{"$rna"}){#3and in this cluster.. 
              if(!$gene){ ## don't yet have a gene..
                $gene = $g;
              }else{ ##should make mergesplit here...
                $algoInvo->addChild(&createMergeSplit($g,$gene,1));
              }
            }else{
              $otherReference++;
              $genesWithReference{"$g"} = 1;
              ##remove from manRevGenes as will NOT use this for this gene..
              delete $manRevGenes{$g};
              next;  ##this is a reference RNA but not in this cluster....leave with gene...
            }
          }
          #set the transcriptunitid to null
          $rna->removeParent($g->getChild('TranscriptUnit'));
          ##don't do next if is rna I want as may cause update when unnecessary
          next if $rnas{"$rna"};
          $rna->setTranscriptUnitId('NULL');
          ##add to list of rna's to submit
          $algoInvo->addChild($rna);
        }
      } 
      ##now get a gene if  don't  have one...
      ##want a gene that is manually reviewed if it does not have other RNA with reference..
      if(!$gene && scalar(@genes) > 0){
        foreach my $v (values%manRevGenes){
          $gene = $v;
          last;
        }
        if(!$gene){
          foreach my $g (@genes){
            next if $genesWithReference{$g};
            $gene = $g;
            last;
          }
        }
      }
      
      if(!$gene){$gene = &createNewGene();} ##don't have any so create a new one...

      ##now foreach  rna in  cluster...add to the gene
      foreach my $r (values %rnas){
        $gene->getChild('TranscriptUnit')->addChild($r);
      }

      ##now mark extra genes deleted...
      foreach my $g (@genes){
        if("$g" eq "$gene"){
          print STDERR "This gene $g eq my gene $gene\n" if $debug;
          next;
        }
        ##may still have reference rna associated...don't want  to  alter gene...
        next if $genesWithReference{$g};
          
        $g->retrieveAllChildrenFromDB(); ##don't do it recursively
        $g->markDeleted(1);
        $gene->addToSubmitList($g);
      }
      ##and submit...
      &submit($algoInvo,$gene);
      $algoInvo->removeAllChildren();
      $algoInvo->undefPointerCache();
      print L "$cluster complete\n";
    }
  }
  close F;
    
  ############################################################
  # return status
  # replace word "done" with meaningful return value/summary
  ############################################################
  my $res = "Processed $ct clusters";
  print STDERR "\n$res\n";
  print L "\n$res\n";
  close L;
  return "$res";
}

##need to  submit the rna children first before the genes..
sub submit {
  my($algoInvo,$gene) = @_;
  $algoInvo->manageTransaction(undef,'begin');
  ##need to submit the RNAs first before genes..
  foreach my $rna ($algoInvo->getAllChildren()){
    print STDERR "Submitting RNA....should have null transcript_unit_id...\n" if $debug;
    $rna->submit(undef,1);
  }
  $gene->submit(undef,1);
#  foreach my $g ($algoInvo->getChildren('Gene')){
#    $g->submit(undef,1);
#  }
  return $algoInvo->manageTransaction(undef,'commit');
}

sub createNewRNAAndSequence {
 my $rs = RNASequence->new({'manually_reviewed' => 0,
                             'is_reference' => 0,
                             'rna_sequence_type_id' => 0 });

  ##last the RNA
  my $rna = RNA->new({'manually_reviewed' => 0,
                      'is_reference' => 0,
                      'name' => 'DOTS Assembly',
                     });
  $rna->addChild($rs);

  ##need to also create a protein entry so can do GOFunction predictions..
  my $prot = Protein->new({ 'manually_reviewed' => 0,
                            'is_reference' => 0 });
  $rna->addChild($prot);

  return $rna;
}

#creates a new gene with a TU child...
sub createNewGene {
  my $gene = Gene->new({'is_reference' => 0,
                        'manually_reviewed' => 0 });
  ##add a transcriptUnit
  $gene->addChild(TranscriptUnit->new({'is_reference' => 0,
                                       'manually_reviewed' => 0 }));
  return $gene;
}

sub createMergeSplit {
  my($o,$n,$is_merge) = @_;
  print STDERR "Creating MergeSplit: Old:".$o->getId().", New:".$n->getId().", is_merge:'$is_merge'\n" if $debug;
  my $ms = MergeSplit->new({'old_id' => $o->getId(),
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
