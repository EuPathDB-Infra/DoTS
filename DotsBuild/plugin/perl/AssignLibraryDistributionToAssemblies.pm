package DoTS::DotsBuild::Plugin::AssignLibraryDistributionToAssemblies;

@ISA = qw(GUS::PluginMgr::Plugin);
use strict;
use GUS::Model::DoTS::AssemblyAnatomyPercent;

# ----------------------------------------------------------------------
# create and initalize new plugin instance.

sub new {
  my $Class = shift;

  my $m = bless {}, $Class;

  my $usage = 'assign library distribution to assemblies';

  my $easycsp =
    [
     {o => 'testnumber',
      t => 'int',
      h => 'number of iterations for testing',
     },
     {o => 'taxon_id',
      t => 'int',
      h => 'taxon_id',
     },
     {o => 'restart',
      t => 'string',
      h => 'restarts: ignores assembies in the AssemblyAnatomyPercent table >= this date',
     },
     {o => 'idSQL',
      t => 'string',
      h => 'SQL query that returns Assembly na_sequence_ids, deletes these from AssemblyAnatomyPercent unless >= restart date and then creates new entries for each na_sequence_id',
     },
     {o => 'update',
      t => 'boolean',
      h => 'Deletes existing rows from AssemblyAnatomyPercent for these na_sequence_ids',
     }
   ];

  $m->initialize({requiredDbVersion => {},
		  cvsRevision => '$Revision$', # cvs fills this in!
		  cvsTag => '$Name$', # cvs fills this in!
		  name => ref($m),
		  revisionNotes => 'make consistent with GUS 3.0',
		  easyCspOptions => $easycsp,
		  usage => $usage
		 });

  return $m;
}

my $debug = 0;

$| = 1;

my $ctx;
my $dbh;
my %multAss;
my %anat;
my %hier;

sub run {
  my $M   = shift;
  $ctx = shift;

  my $stmt;

  print $ctx->{cla}->{'commit'} ? "***COMMIT ON***\n" : "***COMMIT TURNED OFF***\n";
  print "Testing on $ctx->{cla}->{'testnumber'}\n" if $ctx->{cla}->{'testnumber'};

  die "--idSQL is required\n" unless $ctx->{cla}->{idSQL};

  die "--taxon_id is required\n" unless $ctx->{cla}->{taxon_id};

  my $taxon_id = $ctx->{cla}->{taxon_id}; 

  $dbh = $ctx->{self_inv}->getQueryHandle();

  my %ignore;
  if ($ctx->{cla}->{'restart'}) {
    $stmt = $dbh->prepare("select distinct na_sequence_id from DoTS.AssemblyAnatomyPercent where modification_date >= '$ctx->{cla}->{restart}'");
    $stmt->execute();
    while ( my($id) = $stmt->fetchrow_array()) {
      $ignore{$id} = 1;
    }
    print "Restarting: ignoring ".scalar(keys%ignore). " Assemblies\n";
  }

  ##now need to step through the assemblies and get the source_library for each EST
  ##calculating the % contribution...

  $stmt = $dbh->prepare("select anatomy_id,parent_id,hier_level from SRes.Anatomy");
  $stmt->execute();
  while (my ($anatomy_id,$parent_id,$hier_level) = $stmt->fetchrow_array()) {
    print STDERR "$anatomy_id: parent=$parent_id, level=$hier_level\n" if $debug;
    $anat{$anatomy_id} = $parent_id;
    $hier{$anatomy_id} = $hier_level; ##put this into an array!!
  }
  print STDERR "Anatomy_ids: ",scalar(keys%anat),"\n"; # if $debug;

  ##could get those things that are assigned to multiple anat_ids...track number
  ##and when processing take into account somehow...
  $stmt = $dbh->prepare("select dbest_library_id,anatomy_id from DoTS.AnatomyLibrary order by dbest_library_id");
	
  $stmt->execute();

  my $prevLib = 0;
  my @tmp;
  while (my($lib,$anat) = $stmt->fetchrow_array()) {
    if ($lib != $prevLib ) {    ##have array of anat_ids for each multAss lib_id
      @{$multAss{$prevLib}} = @tmp if scalar(@tmp) > 1;
      undef @tmp;
    }
    $prevLib = $lib;
    push (@tmp,$anat);
  }

  print STDERR "have ".scalar(keys%multAss)." multiply assigned libraries\n"; # if $debug;

  ##now get the suckers towork on!!
  $stmt = $dbh->prepare($ctx->{cla}->{idSQL});

  $stmt->execute();

  my $countEntries = 0;
  my $countAssign = 0;
  my @todo;
  print STDERR "Getting Assembly entries to process\n"; # if $debug;ssemb
  my $cte = 0;
  while (my($nas) = $stmt->fetchrow_array()) {
    next if exists $ignore{$nas};
    $cte++;
    print STDERR "Retrieving entry $cte\n" if($cte % 10000 == 0);
    push(@todo,$nas);
    last if ($ctx->{cla}->{'testnumber'} && $cte >= $ctx->{cla}->{'testnumber'}); ##for testing...
  }
  my $totalToDo = scalar(@todo);
  print "$totalToDo entries selected to process for distribution\n";

  undef %ignore;                ##free this memory...

  ##prepare the library query statment to reuse in following loop
  my $libQuery = " select /*+ RULE*/ e.accession,al.dbest_library_id,al.anatomy_id
        from DoTS.AssemblySequence ass,
        DoTS.AnatomyLibrary al, DoTS.EST e, DoTS.Library l
        where ass.assembly_na_sequence_id = ?
        and e.na_sequence_id = ass.na_sequence_id
        and l.library_id = e.library_id
        and al.dbest_library_id = l.dbest_id";

  $stmt = $dbh->prepare($libQuery);

  my $updateStmt = $dbh->prepare('select * from DoTS.AssemblyAnatomyPercent where na_sequence_id = ?');

  foreach my $na_sequence_id (@todo) {
    my %library;
    my %totalLibs;
    my %mult;
    $ctx->{'self_inv'}->undefPointerCache();

    #	print "processing $na_sequence_id entry number $countEntries\n";

    $countEntries++;
		
    ##first delete existing entries if updating...NOTE: this is in separate transaction...
    if ($ctx->{cla}->{update}) {
      $updateStmt->execute($na_sequence_id);
      my @del;
      while (my $row = $updateStmt->fetchrow_hashref()) {
        push(@del,GUS::Model::DoTS::AssemblyAnatomyPercent->new($row));
      }
      if (@del) {
        $ctx->{'self_inv'}->manageTransaction(undef,'begin'); ##starts a transaction...
        foreach my $l (@del) {
          $l->markDeleted();
          $l->submit(undef,1);  ##submit without starting a new transaction...
        }
        $ctx->{'self_inv'}->manageTransaction(undef,'commit'); ##commits a transaction...
      }
    }

    $stmt->execute($na_sequence_id);

    while (my($accession,$library_id,$anatomy_id) = $stmt->fetchrow_array()) {
      $library{$anatomy_id}++;
      $totalLibs{$accession}++;	#for keeping the total number of unique accessions
      if ($totalLibs{$accession} == 2) { ##only record this the first time!!
        $mult{$library_id}++;   ##keeping the number from these multiply assigned libs.
      }
    }
    #3now analyze the libraries.
    my $est_count = scalar(keys%totalLibs);
    print STDERR "IDS: ",join(', ',keys%library),"\n" if $debug;
    if ($est_count > 0) {
      #      last;  ##for debugging the multiple problem!!
      $countAssign += &analyzeLib($na_sequence_id,$taxon_id,$est_count,\%mult,%library);
    }
    print "Processing $na_sequence_id: $countEntries, $countAssign assigned, ",($totalToDo - $countEntries)," remaining\n" if $countEntries % 100 == 0; 

  }

  ##results here.....
  my $results = "run finished, processed $countEntries and assigned distribution to $countAssign DOTS RNAs";

  $ctx->{self_inv}->closeQueryHandle();

  return $results;
}





1;

sub analyzeLib{
  my($na_sequence_id,$taxon_id,$total,$mult,%libs) = @_;
  my %totLibs;
  my %multAnat;
  my %track;
  my %haveSub;                  ##tracks if have subtracted
  if ($debug) { 
    print STDERR "Analyzing $na_sequence_id: $total ESTs\n";
    print STDERR "Multiple lib_ids = ",join(', ',keys%{$mult}),"\n";
  }
	
  ##create hash of multiple assigned things!!
  foreach my $key (keys %{$mult}) {
    foreach my $id (@{$multAss{$key}}) {
      $multAnat{$id}->{$key} = $mult->{$key}; ##datastructure is hash with anat_id keys
      ##containing hash with lib_id keys value is that combo....
      print STDERR "\(Lib_id $key, Anat_id $id\) = $mult->{$key}\n" if $debug;
    }
  }
  
  ##need to sort by reverse hierarchy level.
  ##need to add the children to the parent even if not directly the parent.
  foreach my $key (sort{$hier{$b} <=> $hier{$a}}keys%libs) {
    print STDERR "$key Level $hier{$key}\n" if $debug;
    my $parID = $anat{$key};
    $totLibs{$key} += $libs{$key}; ##starts the percolation up...
    until ($parID == 0) {
      ##deal with multiples.....
      if (exists $multAnat{$key}) {
        foreach my $id (keys%{$multAnat{$key}}) {
          next if exists $haveSub{$id}->{$key}; ##have already subtracted this one...
          ##need to track each one and if seen then substract the value from $libs{$key}
          if (exists $track{$id}->{$parID}) {
            $libs{$key} -= $multAnat{$key}->{$id};
            $haveSub{$id}->{$key} = 1; ##have subtracted this one....don't do again!!
          }
          $track{$id}->{$parID} = 1;
        }
      }
      $totLibs{$parID} += $libs{$key};
      print STDERR "Adding $libs{$key} from $key to $parID = $totLibs{$parID}\n" if $debug;
      $parID = $anat{$parID};   ##gets the parent...
    }
  }
  my $totalPercent = 0;
  my @libdist;
  foreach my $key (keys%totLibs) { ##this should do it...
    my $percent = ($totLibs{$key}/$total)*100;
    $totalPercent += $percent;
    push(@libdist,&makeLibDist($na_sequence_id,$taxon_id,$key,$percent,$totLibs{$key},$total)) if $percent > 0;
    print STDERR "\tAnat_id: $key = $percent percent\n" if $debug;
  }
  print STDERR "\tTotalPercent: $totalPercent\n" if $debug;
	
  if (scalar(@libdist) > 0) {
    ##submit the distribution here!!
    $ctx->{'self_inv'}->manageTransaction(undef,'begin'); ##starts a transaction...
    foreach my $l (@libdist) {
      $l->submit(undef,1);      ##submit without starting a new transaction...
    }
    $ctx->{'self_inv'}->manageTransaction(undef,'commit'); ##commits a transaction...
    return 1;                   ##returns 1 if it does inserts...
  }
}

sub makeLibDist {
  my($na_sequence_id,$taxon_id,$anat_id,$percent,$anatomy_ests,$est_count) = @_;
  my$libH = GUS::Model::DoTS::AssemblyAnatomyPercent->
    new({'na_sequence_id' => $na_sequence_id,
	 'anatomy_id' => $anat_id,
	 'est_count' => $est_count,
	 'anatomy_ests' => $anatomy_ests,
	 'percent' => $percent,
	 'taxon_id' => $taxon_id});
  return $libH;
}

__END__

=pod
=head1 Description
B<Template> - a template plug-in for C<ga> (GUS application) package.

=head1 Purpose
B<Template> is a minimal 'plug-in' GUS application.

=cut
