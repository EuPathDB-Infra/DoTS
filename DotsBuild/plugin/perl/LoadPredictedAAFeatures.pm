#!/usr/local/bin/perl

# ----------------------------------------------------------
# ImportPlasmoDBAAFeatures.pm
#
# Description : ga-plugin used to import predicted AA features
#               from a tab-delimited file.
#
#              Note the following:
#
# Modified  By               Description
# _________________________________________________________
#
# 5/24/00   Sharon Diskin    Created
#11/30/00   Martin Fraunholz Modified for PlasmoDB
#12/02/02   Jonathan Schug   Modified for DoTS.
#----------------------------------------------------------
package DoTS::DotsBuild::Plugin::LoadPredictedAAFeatures;
@ISA = qw(GUS::PluginMgr::Plugin);

use strict;
use FileHandle;
use GUS::Model::DoTS::AALocation;
use GUS::Model::DoTS::TranslatedAASequence;
use GUS::Model::DoTS::PredictedAAFeature;
use GUS::Model::DoTS::PfamEntry;
use GUS::Model::DoTS::SignalPeptideFeature;
use Disp;



#----------------------------------------------------------
# Global Variables:
#----------------------------------------------------------

my $seq_source_id = "";;
my $seq_description = "";
my %finished;
my $countInserts = 0;
my %counter =();
my $verbose = 0;
my %alg_id=(
						'TOPPRED2'    => '306',
						'SignalP'     => '305',
						'PATMATMOTIFS'=> '304',
						'HMMPFAM'     => '303',
						'TMAP'        => '302',
						'TMPRED'      => '301',
						'TMHMM2'      => '7392'
           );
#-----------------------------------------------------------
# GUSApplication
#-----------------------------------------------------------

sub new {
	 my $class = shift;
	 my $self = {};
	 bless($self, $class);
	 my $usage = 'A package to insert predicted AA features from a tab delimited file.';

	 my $easycsp = 

	   [
	    {
	     o => 'testnumber',
	     t => 'int',
	     h => 'Number of iterations for testing',
	    },

	    {
	     o => 'restart',     
	     t => 'string',
	     h => 'For restarting script...takes list of row_alg_invocation_ids to exclude',
	    },

	    {
	     o => 'filename',
	     t => 'string',
	     h => 'Name of file containing the predicted aa features',
	    },

	    {
	     o => 'reload',
	     t => 'boolean',
	     h => 'For loading new Ids and scores for HmmPfam entries',
	    },
	    {
	     o => 'project_id',
	     t => 'int',
	     h => 'applicable project_id from core.projectinfo table',
	    }
	   ];

	    $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$', # keyword filled in by cvs
		     cvsTag => '$Name$ ',             # keyword filled in by cvs
		     name => ref($self),
		     revisionNotes => 'first pass of plugin for GUS 3.0',
		     easyCspOptions => $easycsp,
		     usage => $usage
		    }); 
  return $self;
}

sub run {
	 my $self   = shift;
	 $self->log("Testing on " . $self->getArgs()->{'testnumber'}."\n") if $self->getArgs()->{'testnumber'};
	 #GusApplication::Log('INFO', 'RAIID',  $ctx->{self_inv}->{algorithm_invocation_id});
	 $self->getArgs()->{'commit'} ? $self->log("***COMMIT ON***\n") : $self->log("**COMMIT TURNED OFF**\n");

	 my $aa_sequence_id;

	 if (!-e $self->getArgs()->{'filename'}) {
			die "You must provide a valid tab delimited file on the command line\n";
	 }
	 if ($self->getArgs()->{'debug'}) {
			$verbose=1;
	 }
	 if ($self->getArgs()->{'verbose'}) {
			$verbose=1;
	 }

	 my $f  = $self->getArgs()->{'filename'} =~ /\.Z$/ ? "zcat".$self->getArgs()->{'filename'}."|" : $self->getArgs()->{'filename'};
	 my $fh = FileHandle->new($f);
	 print STDERR "Could not open". $self->getArgs()->{'filename'}.": $!\n\n" unless $fh;

	 my %aa_sequence_ids_seen = ();
	 my $rows_n = 0;
	 my $seq;

	 while (<$fh>) {
			chomp;
			#
			# skip comments and blank lines, shouldn't be in there anyway
			#
			next if ( (/^\s*\#/) || (/^\s+$/) );
			print STDERR $_ ."\n" if $self->getArgs()->{'debug'};
			my @tmp = split(/\t/,$_);
			#GusApplication::Log('INPUT', @tmp);
			$rows_n++;
			$aa_sequence_ids_seen{$tmp[0]}++;

			if (my $seq = $self->get_object(@tmp)) {
				 $self->process($seq,@tmp);
				 Disp::Display(\%counter) if $rows_n % 1000 == 0;
				 $seq->undefPointerCache;
			}
			if ($self->getArgs()->{'testnumber'} && $rows_n > $self->getArgs()->{'testnumber'}) {
			  last();
			}
	 } # end while LOOP
	 $fh->close;

	 $counter{sequence} = scalar keys %aa_sequence_ids_seen;

	 #
	 # generate summary
	 #

	 my $results= "Run finished: Processed $counter{sequence} peptides\n";
	 foreach my $key (keys %counter) {
			$results .= $counter{$key} . " $key features\n";
	 }

	 print STDOUT "\n" . $results . "\n\n"; 
	 return $results;
}






############################################################
# Subroutines
############################################################
sub get_object{
	 my ($self,@tmp)=@_;
	 my $seq;
	 my $aa_sequence_id;
	 my $project_id=$self->getArgs()->{'project_id'};

	 if ($tmp[0]=~m/^(\d+)/) {
			my $id=$1; 
			$seq = GUS::Model::DoTS::TranslatedAASequence->new({'aa_sequence_id' => $id});
			$self->log ("Found UID $id !\n") if $self->getArgs()->{'verbose'};
	 } else {
			$self->log ("Illformed AaSequenceId $tmp[0]\n");
	 }

	 return $seq;
}

###########################################################
sub error{
	 my $errcode=shift;
	 my $errtxt=shift;
	 print STDERR "###########################################\n";
	 print STDERR $errtxt ."\n\n";
	 print STDERR "\n\ngoing on to next sequence\n\n";
	 print STDERR "and deleting PointerCache\n";
	 print STDERR "###########################################\n";
}

############################################################
sub process {
  my ($self,$seq,@tmp) = @_;

  my $aa_sequence_id = $tmp[0];

  if (!$seq) {
    die ("ERROR: sequence object does not exist!\n");
  }

  #
  # Set the fields that we need.
  #
  my $source_id = $tmp[0];
  $source_id=~s/\.pep$//;			#in case it is still there
  my $alg_name = $tmp[1];

  my $feature;

  if ($self->getArgs()->{'verbose'} && $source_id) {
    print STDERR "source_id is $source_id\n";
    print STDERR "alg_name is $alg_name\n\n";
  }


  ## WILL HAVE TO CHANGE THAT TO A GFF FORMAT INPUT, 
  ## SO I CAN GET RID OF THE ELSIF LOOPS
  ## On the other hand: GFF sucks !!


  if ($alg_name eq "SignalP") {
    $feature = $self->createNewSignalPeptideFeature(@tmp);
    if ($feature) {
      $counter{'signal'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  }														# end SIGNALP

  # ........................................
  if ($alg_name=~m/SignalPHMM/i) {
    $feature = $self->UpdateSignalPeptideHMMFeature(@tmp);

    if ($feature) {
      $counter{'signalhmm'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  }														# end SIGNALP

  # ........................................
  elsif ($alg_name eq "TMPRED") {

    my ($file,$algorithm,$featuretype,$score,$start,$stop,$center,$orientation,$helix_min_len,$helix_max_len)=@tmp;
    $feature = $self->createNewPredictedAAFeature($source_id, 'transmembrane region',
						  $alg_name, 'TMhelix' , $score, $start, $stop);
    if ($feature) {
      $counter{'tmpred'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }

  }														# END TMPRED


  # ........................................
  elsif ($alg_name=~/TMHMM2/) {

    my ($file,$algorithm,$featuretype,$start,$stop,undef,undef,undef)=@tmp;
    $algorithm=~s/\.\d+$//;		#remove trailing .0

    $feature = $self->createNewPredictedAAFeature($file, 'transmembrane region',
						  $algorithm, 'TMhelix', '',
						  $start, $stop);
    if ($feature) {
      $counter{'tmhmm'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  } # END TMHMM2


  elsif ($alg_name eq "TOPPRED2") {

    my ($file,$algorithm,$featuretype,$helix_num,$helix_begin,$helix_end,$helix_score,$helix_certainity,$hydrophob_file,$cyt_ext_file,$org_type,$charge_pair,$full_wdw,$core_wdw,$num_res,$crit_len,$upcandcutoff,$lowcandcutoff)=@tmp;
    $feature = $self->createNewPredictedAAFeature($source_id, 'transmembrane region',
						  $alg_name, 'TMhelix', $helix_score, $helix_begin, $helix_end);

    if ($feature) {
      $counter{'toppred'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  }														# END TOPPRED

  elsif ($alg_name eq "TMAP") {
    my ($file,$algorithm,$featuretype,$start,$stop,$tmap_helix_len)=@tmp;
    $feature = $self->createNewPredictedAAFeature($source_id,'transmembrane region',
						  $alg_name, 'TMhelix', '', $start, $stop); ########### no score available
    if ($feature) {
      $counter{'tmap'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++; 
    }
  }														# END TMAP

  elsif ($alg_name eq "PFSCAN") {

    my ($file,$algorithm,$featuretype,$database,$score1,$score2,$start,$stop,$donknow1,$donknow2,$PSID,$PSABBR,$PS_TITLE)=@tmp;
    $feature = $self->createNewPredictedAAFeature($source_id, "PROSITE motif $PSID, accession $PS_TITLE",
						  $alg_name, 'PROSITE motif', $score1, $start, $stop);

    if ($feature) {
      $counter{'pfscan'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  }														#END PFSCAN

  elsif ($alg_name=~/PATMATMOTIF/) { ######### WORKAROUND DUE TO TYPO -> patmatmotifs

    my ($file,$algorithm,$featuretype,$database,$start,$stop,$PSID,$PSACC)=@tmp;
    $feature = $self->createNewPredictedAAFeature($source_id, "PROSITE motif $PSID, accession $PSACC",
						  "PATMATMOTIFS", 'PROSITE motif', '0', $start, $stop);
    if ($feature) {
      $counter{'patmat'}++;
    } else {
      print STDERR "ERROR: Unable to create PredictedAAFeature for $source_id $alg_name\n\n";
      $counter{'unprocessed'}++;
    }
  }														#END PATMATMOTIFS

  elsif ($alg_name eq "HMMPFAM") {

    my ($file,$algorithm,$featuretype,$database,$accession,$model,$info,$domain,$start,$stop,$hmmf,$hmmt,$score,$evalue)=@tmp;

    # return if there is no matching entry in PFAM
    my $pf = GUS::Model::DoTS::PfamEntry->new({'release' => '7.5',
			     'accession' => $accession});
    $pf->retrieveFromDB() || print STDERR "Could not retrieve PFAMEntry for SOURCE $source_id :  ACC: $accession\n\n" && return undef;

    my $motif_id=$pf->getPfamEntryId;

    my ($feature,$aafid)=$self->existenceHmmPfamFeature($aa_sequence_id, $alg_name, 'PFAM motif', 
							$score, $start, $stop);
    $self->log("Got HmmPfamFeature with ID $aafid....\n\n") if $self->getArgs()->{'verbose'};

    ## if the feature exists, we might want to update it (e.g. is values are missing)....
    if ($feature) {
      
      ## ... but only if 'reload' is set !!! 
      unless ($self->getArgs->{'reload'}) {
	
	## else return to parser !
	print STDERR "\n\nCannot create the new HmmPFAM Feature. Returning....\n\n";
	$counter{'existsPfam'}++;
	return;
      }
      
      print STDERR "\n\nSetting MotifID $motif_id for FeatureId $aafid....\n\n";
      $feature->set('score', $score);
      $feature->set('motif_id', $motif_id) ; 
      $feature->setParent($pf);
      
      $counter{'updatePfam'}++;
      
    } else {
      ## if the feature does not exist, we create it
      
      
      my $newFeature = $self->createNewPredictedAAFeature( $source_id, $info,
							   $alg_name, 'PFAM motif',
							   $score, $start, $stop, $motif_id);
      
      if ($newFeature) {
	print STDERR "\n\nNo Entry yet!!  Creating new HMMPFAM Feature for $source_id...\n\n";
	$seq->addChild($newFeature);
	$newFeature->setParent($pf);
	$counter{'newPfam'}++;
      } else {	
	print STDERR "\n\nCannot create the new HmmPFAM Feature.\n\n";
	$counter{'missedPfam'}++;
      }
    }
    
  }
  
  $feature->submit if $feature;
  
}																

######################################################################

sub existenceAAFeature{
	 my ($self,$source_id, $alg_name, $feature_name, $start, $stop) = @_;

	 my $dbh = $self->getQueryHandle();
	 my $sql = <<Sql;

  select paf.name, paf.algorithm_name, paf.description
       , paf.score, aal.start_min, aal.end_max
	  from dots.PredictedAAFeature paf
       , dots.AALocation         aal
	 where paf.aa_sequence_id = $source_id
     and paf.algorithm_name = '$alg_name'
	   and paf.name           = '$feature_name'
     and aal.aa_feature_id  = paf.aa_feature_id
	   and aal.start_min      = $start
	   and aal.end_max        = $stop

Sql

	 #print STDERR "$sql\n" if $ctx->{'debug'};

	 my $stmt = $dbh->prepare( $sql );

	 $stmt->execute();

	 my $countRows = 0;
	 while (my $row = $stmt->fetchrow_hashref('NAME_lc')) {
      $countRows++;
	 }

	 if ($countRows > 0) {
			print STDERR ">>>>>>>>>>Skipping feature $feature_name for $source_id (from $start to $stop)....\n\n\n";
			print STDERR "$countRows rows returned for :\n$sql\n\n";
			return $countRows;
	 } else {
			return undef;
	 }

}																#end Sub

#####################################################################

sub existenceHmmPfamFeature{
	 my ($self,$aa_sequence_id, $alg_name, $feature_name,$score, $start, $stop) = @_;

	 my $project_id=$self->getArgs()->{'project_id'};
	 my $dbh = $self->getQueryHandle();
	 my $sql = "select paf.*
	        from dots.TranslatedAAFeature taf, dots.PredictedAAFeature paf, dots.AALocation aal
	        where taf.aa_sequence_id=$aa_sequence_id
	        and taf.aa_sequence_id = paf.aa_sequence_id
	        and paf.aa_feature_id = aal.aa_feature_id
                and paf.algorithm_name = \'$alg_name\'
                and paf.score = $score
	        and aal.start_min=$start
	        and aal.end_max=$stop
	        and paf.name = \'$feature_name\'";

	 my $stmt = $dbh->prepare( $sql );

	 #print STDERR "$sql\n";
	 $stmt->execute();
	 my $paf;
	 my $aafid;
	 my $countRows = 0;
	 while (my $row = $stmt->fetchrow_hashref('NAME_lc')) {
      $paf=GUS::Model::DoTS::PredictedAAFeature->new($row);
      $aafid=$$row{'aa_feature_id'};
      $countRows++;
	 }

	 if ($countRows==1) {
      print STDERR "\n>>>>> SQL returned  1 row  <<<<<<\n\n";
      return ($paf,$aafid);
	 } elsif ($countRows>1) {
      print STDERR "#### $countRows rows returned, featurequery is not stringent enough \n\n\n"; 
      return undef;
	 } else {
      print STDERR "No rows returned for query\n\n";
			return undef;
	 }

}																#end Sub

## END INTRODUCED SECTION
#################################################################


#
# Create the new PredictedAAFeature
#
sub createNewPredictedAAFeature {
  my ($self,$source_id, $description, $alg_name, $feature_name, $score, $start, $stop, $motif_id) = @_;


  # check, if feature already exists
  # rather do this with the method call retrieveFromDB <- see below 
  # since the constraints of retrieveFromDB() do not check for aa_sequence_id or AAlocation

  my $exists=$self->existenceAAFeature($source_id, $alg_name, $feature_name, $start, $stop);
  if ($exists) {
    print STDERR "Feature exists...\n\n";
    return undef;
  }

  my $newFeature = GUS::Model::DoTS::PredictedAAFeature->new({
					    'aa_sequence_id'         => $source_id,
					    'is_predicted'           => 1,
					    'review_status_id'       => 0,
					    'algorithm_name'         => $alg_name,
					    'prediction_algorithm_id'=> $alg_id{$alg_name},
					    'description'            => $description,
					    'name'                   => $feature_name,
					    'score'                  => $score
					   });

  #
  # Create the AALocation object and add it as a child of
  # the created PredictedAAFeature.
  #

  my $aa_location = $self->createNewAALocation($start, $stop);
  $newFeature->addChild($aa_location) if $aa_location;

  return $newFeature;
}


##############################################
sub createNewAALocation {
  my($self,$start,$end) = @_;
  return undef if (!$start || !$end);
  my $aa_loc = GUS::Model::DoTS::AALocation->new({'start_min' => $start,
						  'start_max' => $start,
						  'end_min'   => $end,
						  'end_max'   => $end});
  return $aa_loc;
}
###############################################

sub createNewSignalPeptideFeature {

  ### The feature will be created using the NN data from SignalP
  ### UpdateSignalPeptideFeature will add additional information:
  ### The signal_probability and the anchor_probability
  ### Ideally we would have to have two tables, storing the different sets od data,
  ### so view this as a workaround

  my($self,$source_id,$algorithm,$featuretype,$maxC_position,$maxC_value,$maxC_cutoff,$maxC_conclusion,$maxY_position,$maxY_value,$maxY_cutoff,$maxY_conclusion,$maxS_position,$maxS_value,$maxS_cutoff,$maxS_conclusion,$meanS_start,$meanS_stop,$meanS_value,$meanS_cutoff,$meanS_conclusion,$quality,$signal) = @_;

  # check, if feature already exists. Do this with specific SQL
  # rather than with the retrieveFromDB() call, 
  # since the constraints of retrieveFromDB() do not check for aa_sequence_id or AAlocation

  my ($exists,$aafid) = $self->existenceSPFeature($source_id, $algorithm, $featuretype,
						  $meanS_start, $meanS_stop);
  if ($exists) {
    $self-log ("SignalPeptideFeature with aa_feature_id $aafid exists. Skipping ...\n\n") if $self->getArgs()->{'verbose'};
    return undef;
  }


  #
  # Create the new SignalPeptideFeature if it doe not exist yet.
  #
  my $newSignalPeptide = GUS::Model::DoTS::SignalPeptideFeature
    ->new({
	   'aa_sequence_id'    => $source_id,
	   'is_predicted'      => 1,
	   'review_status_id' => 0,
	   'algorithm_name'    => $algorithm,
	   'prediction_algorithm_id' => '305',
	   'description'       => 'Signal Peptide',
	   'name'              => $featuretype,
	   'maxc_score'        => $maxC_value,
	   'maxc_conclusion'   => $self->conclude($maxC_conclusion),
	   'maxy_score'        => $maxY_value,
	   'maxy_conclusion'   => &conclude( $maxY_conclusion),
	   'maxs_score'        => $maxS_value,
	   'maxs_conclusion'   => $self->conclude( $maxS_conclusion),
	   'means_score'       => $meanS_value,
	   'means_conclusion'  => $self->conclude( $meanS_conclusion),
	   'num_positives'     => $quality,
	  });
  
  #
  # Create the AALocation object and add it as a child of
  # the created PredictedAAFeature.
  #
  my $aa_location = $self->createNewAALocation($meanS_start, $meanS_stop);
  $newSignalPeptide->addChild($aa_location) if $aa_location;
  
  $self->log ($newSignalPeptide->toString()."\n") if $self->getArgs()->{'debug'};
  
  
  ### DOES NOT WORK FOR THAT PURPOSE    
  #    #if ($newSignalPeptide->retrieveFromDB()){
  #	print STDERR "$algorithm feature exists for $source_id !!!\n\n" if $ctx->{'verbose'};
  #	return undef;
  #    }
  ##################################
  
  return $newSignalPeptide;
}

####################################

sub UpdateSignalPeptideHMMFeature {

  my ($self,$source_id, $algorithm, $featuretype,
      $prediction, $SPP, $SAP, $CSP, $start, $signal
     ) = @_;

  #    if ($prediction=~/Non-/){
  #	# do not load or attempt to load entries for non secretory proteins
  #	print STDERR "No SP (HMM) in dataset. Returning \'undef\' ... \n\n" if $ctx->{'verbose'};
  #	return undef;
  #    }else{
  my ($SignalPeptideFeature, $aafid) =
    $self->existenceSPFeature($source_id, $algorithm, $featuretype, 1, $start);

  if ($SignalPeptideFeature) {
    $self->log("UPDATE $source_id\n") if $self->getArgs()->{'verbose'};
    $SignalPeptideFeature->set('anchor_probability', $SAP);
    $SignalPeptideFeature->set('signal_probability', $SPP);
    return $SignalPeptideFeature;
  } else {
    $self->log("NOSIGNAL $source_id\n");
    return undef;
  }
}

######################

sub conclude{
  my($self) = @_;
  my $value=shift;
  my $answer=($value=~m/yes/i) ? 1 : 0;
  return $answer;
}

######################

sub existenceSPFeature{
  my ($self,$source_id, $alg_name, $feature_name, $start, $stop) = @_;
  
  my $project_id = $self->getArgs()->{'project_id'};
  my $dbh = $self->getQueryHandle();
  
  #  my $sql = "select spf.*
  #	        from dots.GeneFeature gf, dots.RNAFeature rf,
  #	        dots.TranslatedAAFeature taf,
  #	        dots.SignalPeptideFeature spf, dots.AALocation aal, dots.ProjectLink pl
  #	        where gf.source_id = \'$source_id\'
  #	        and rf.parent_id = gf.na_feature_id
  #	        and rf.na_feature_id = taf.na_feature_id
  #	        and taf.aa_sequence_id = spf.aa_sequence_id
  #	        and spf.aa_feature_id = aal.aa_feature_id
  #               and spf.algorithm_name = \'$alg_name\'
  #	        and aal.start_min=$start
  #	        and aal.end_max=$stop
  #	        and spf.name = 'SIGNAL'
  #	        and gf.na_feature_id = pl.id
  #	        and pl.table_id = 108
  #	        and pl.project_id = $project_id";
  
  #	 my $sql = "select spf.*
  #	        from dots.GeneFeature gf, dots.RNAFeature rf,
  #	        dots.TranslatedAAFeature taf,
  #	        dots.SignalPeptideFeature spf, dots.ProjectLink pl
  #	        where gf.source_id = \'$source_id\'
  #	        and rf.parent_id = gf.na_feature_id
  #	        and rf.na_feature_id = taf.na_feature_id
  #	        and taf.aa_sequence_id = spf.aa_sequence_id
  #                and spf.algorithm_name = 'SignalP'
  #	        and spf.name = 'SIGNAL'
  #	        and gf.na_feature_id = pl.id
  #	        and pl.table_id = 108
  #	        and pl.project_id = $project_id";

  my $sql = "select spf.*
   from dots.SignalPeptideFeature spf
      , dots.AALocation           aal
  where spf.aa_sequence_id = $source_id
    and aal.aa_feature_id  = spf.aa_feature_id
    and aal.start_min      = $start
    and aal.end_max        = $stop
    and spf.name           = 'SIGNAL'
    and spf.algorithm_name = 'SignalP'";

  my $stmt = $dbh->prepareAndExecute( $sql );

  my $countRows = 0;
  my $aafeatureid;
  my $spf;
  while (my $row = $stmt->fetchrow_hashref('NAME_lc')) {
    $spf = GUS::Model::DoTS::SignalPeptideFeature->new($row);
    $aafeatureid = $row->{'aa_feature_id'};
    $countRows++;
  }

  if ($countRows>0) {
    $self->log("EXISTS $source_id, $feature_name, $start, $stop, $countRows\n" );
    return ($spf,$aafeatureid);
  } else {
    return undef;
  }

}



1;




__END__

