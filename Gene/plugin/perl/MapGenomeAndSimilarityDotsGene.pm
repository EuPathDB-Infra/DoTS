package DoTS::Gene::Plugin::MapGenomeAndSimilarityDotsGene;
@ISA = qw( GUS::PluginMgr::Plugin);

use strict 'vars';

use CBIL::Util::PropertySet;

use GUS::PluginMgr::Plugin;

sub new {
    my ($class) = @_;
    my $self = {};
    bless($self,$class);

    my $purposeBrief = 'Map genome dots genes to similarity dots genes';
  
    my $purpose = $purposeBrief;

    my $tablesDependedOn = [];
    my $tablesAffected = [];

    my $howToRestart = <<RESTART; 
Run in specific modes to rerun failed steps.
RESTART
;
    my $failureCases = <<FAILURE_CASES;
Expensive queries may upset DBMS.
FAILURE_CASES
;
    my $notes = &_getNotes;

    my $documentation = {purpose=>$purpose, purposeBrief=>$purposeBrief,tablesAffected=>$tablesAffected,tablesDependedOn=>$tablesDependedOn,howToRestart=>$howToRestart,failureCases=>$failureCases,notes=>$notes};
;

    my $mode_sub =  sub { my $mode = $_[1];
			  die "mode is $mode. if set, it must be prep|init|oneone|onemore|transfer"
			      if $mode && $mode !~ /^(prep|init|oneone|onemore|transfer)$/i ;
		      };
    my $argsDeclaration  =
	[
	 stringArg({name => 'mode',
		    descr => 'what step to do (all steps will be done if not set)',
		    constraintFunc=> $mode_sub,
		    reqd  => 0,
		    isList => 0, 
		}),

     stringArg({name => 'temp_login',
		descr => 'login for temp table usage',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0, 
	    }),

    stringArg({name => 'temp_map_table',
		descr => 'temp table for intermediate mapping info cache',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'GenomeAndSimilarityDotsGeneMap',
		isList => 0 
		}),

    stringArg({name => 'genome_dots_gene_cache',
		descr => 'table that caches info about genome dots genes',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'GenomeDotsGene',
		isList => 0 
		}),

    stringArg({name => 'genome_dots_transcript_cache',
		descr => 'table that caches info about genome dots transcripts',
		constraintFunc=> undef,
		reqd  => 0,
		default => 'GenomeDotsTranscript',
		isList => 0 
		}),

    integerArg({name => 'taxon_id',
		descr => 'taxon id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),

    integerArg({name => 'genome_db_rls_id',
		descr => 'genome external database release id',
		constraintFunc=> undef,
		reqd  => 1,
		isList => 0 
		}),
    ];

    $self->initialize({requiredDbVersion => {},
		     cvsRevision => '$Revision$',
		     cvsTag => '$Name$',
		     name => ref($self),
		     revisionNotes => '',
		     argsDeclaration => $argsDeclaration,
		     documentation => $documentation
		    });
    return $self;
}


sub run {
    my $self = shift;
    
    $self->logAlgInvocationId();
    $self->logArgs();

    my $dbh = $self->getQueryHandle();
    my $tmpTab = $self->getArg('temp_login') . '.' . $self->getArg('temp_map_table');
    my $gdgTab = $self->getArg('temp_login') . '.' . $self->getArg('genome_dots_gene_cache');
    my $gdtTab = $self->getArg('temp_login') . '.' . $self->getArg('genome_dots_transcript_cache');
    my $taxonId = $self->getArg('taxon_id');
    my $genomeId = $self->getArg('genome_db_rls_id');

    my $mode = $self->getArg('mode');
    $self->prepareForMapping($dbh, $gdtTab, $taxonId, $genomeId) if !$mode || $mode =~ /^prep$/i;
    $self->initMapping($dbh, $tmpTab, $gdtTab, $taxonId, $genomeId) if !$mode || $mode =~ /^init$/i;
    $self->uniqueMapOneToOne($dbh, $tmpTab, $taxonId, $genomeId) if !$mode || $mode =~ /^oneone$/i;
    $self->uniqueMapOneToMore($dbh, $tmpTab, $gdgTab, $gdtTab, $taxonId, $genomeId)
	if !$mode || $mode =~ /^onemore$/i;
    $self->transferGeneId($dbh, $tmpTab, $gdgTab, $taxonId, $genomeId)
	if !$mode || $mode =~ /^transfer$/i;

    return "done";
}

############### subroutines ###########################

sub _getNotes {
    my $notes =<<NOTES;
To fill in later.
NOTES
;
    $notes;
}

sub prepareForMapping {
    my ($self, $dbh, $gdtTab, $taxonId, $genomeId) = @_;

    $dbh->do("alter table $gdtTab add similarity_dots_gene_id NUMBER(10)") or print "";
    my ($sql, $sth);

    # get an autocommiting query handle
    $dbh = $self->getDb()->getQueryHandle(1);
=pod    
    $sql = "select na_sequence_id from $gdtTab "
	. "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $self->log("get dt ids in $gdtTab for which to find similarity dots gene id: sql=$sql");
    $sth = $dbh->prepareAndExecute($sql);
    my %dt_ids;
    while (my ($id) = $sth->fetchrow_array) { $dt_ids{$id} = undef; }
    $self->log("found " . scalar(keys %dt_ids) . " unique dt ids");
=cut

    my $tmpTab = "DtDgTmp";
    $sql = "create table $tmpTab as select rf.na_sequence_id, r.gene_id "
         . "from DoTS.RnaFeature rf, DoTS.RnaInstance ri, DoTS.RNA r "
         . "where rf.na_feature_id = ri.na_feature_id and ri.rna_id = r.rna_id";
    $self->log("caching dt id vs sDG id mapping: sql=$sql");
    $dbh->sqlexec($sql); 
    $sql = "create index $tmpTab" . "_ind01 on $tmpTab" . "(na_sequence_id)";
    $self->log("creating index: sql=$sql");
    $dbh->sqlexec($sql);
    $dbh->sqlexec("analyze table $tmpTab compute statistics");

    $sql = "update $gdtTab gdt set gdt.similarity_dots_gene_id = "
         . "(select tmp.gene_id from $tmpTab tmp where gdt.na_sequence_id = tmp.na_sequence_id) "
         . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";
    $self->log("updating sDG id in $gdtTab: sql=$sql");
    $dbh->sqlexec("$sql");

    $self->log("deleting temp table $tmpTab");
    $dbh->do("drop table $tmpTab");
=pod
    my $c = 0;
    foreach my $dt_id (keys %dt_ids) {
	$sql = "select r.gene_id from DoTS.RnaFeature rf, DoTS.RnaInstance ri, DoTS.RNA r "
	    . "where rf.na_sequence_id = $dt_id "
	    . "and rf.na_feature_id = ri.na_feature_id and ri.rna_id = r.rna_id";
	$sth = $dbh->prepareAndExecute($sql);
	my $gene_id = $sth->fetchrow_array;
	if (!$gene_id) {
	    $gene_id = 'NULL';
	    $self->log("WARNING: not similarity gene id found for DT.$dt_id");
	}

	$sql = "update $gdtTab set similarity_dots_gene_id = $gene_id "
	    . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	    . "and na_sequence_id = $dt_id";
	$dbh->sqlexec("$sql");

	unless (++$c % 500) {
	    $self->log("set sDG id for $c entries");
	    $dbh->commit();
	}
    }
=cut

}

sub initMapping {
    my ($self, $dbh, $tmpTab, $gdtTab, $taxonId, $genomeId) = @_;

    $dbh->do("drop table $tmpTab") or print "";
    
    my $selSql = "select distinct taxon_id, genome_external_db_release_id, "
	. "genome_dots_gene_id, similarity_dots_gene_id, 0 as selected from $gdtTab "
        . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId";

    my $sql = "create table $tmpTab as ($selSql)";
    $self->log("creating $tmpTab with initial id mapping: sql = $sql");
    $dbh->sqlexec($sql);
    $self->log("done");
}

# CASE 1: set the one to one mappings first
sub uniqueMapOneToOne {
    my ($self, $dbh, $tmpTab, $taxonId, $genomeId) = @_;
    
    my $sql = "select genome_dots_gene_id from $tmpTab "
	. "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	. "group by genome_dots_gene_id having count(distinct similarity_dots_gene_id) = 1 "
	. " intersect  "
	. "select max(genome_dots_gene_id) as genome_dots_gene_id from $tmpTab "
	. "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	. "group by similarity_dots_gene_id having count(distinct genome_dots_gene_id) = 1"; 
    $self->log("Finding the one to one mappings, sql = $sql");
    my $sth = $dbh->prepareAndExecute($sql);

    $self->log("setting 1 to 1 mappings in db, 500 at a time...");
    my ($gdgids, $count);
    while (my ($id) = $sth->fetchrow_array) {
	$gdgids .=  $id . ',';
	if (++$count % 500 == 0) {
	    &setSelected($dbh, $tmpTab, $taxonId, $genomeId, $gdgids);
	    $gdgids = undef;
	    $dbh->commit();
	    $self->log("set $count");
	}
    }
    &setSelected($dbh, $tmpTab, $taxonId, $genomeId, $gdgids);
    $dbh->commit();
    $self->log("set $count");
}

# CASE 2: for aligned_gene_id's mapping to more than one gene_id's
sub uniqueMapOneToMore {
    my ($self, $dbh, $tmpTab, $gdgTab, $gdtTab, $taxonId, $genomeId) = @_;

    my $sql = "select m.genome_dots_gene_id, m.similarity_dots_gene_id "
	. "from $tmpTab m, $gdtTab gdt, $gdgTab gdg, DoTS.Assembly a "
	. "where m.taxon_id = $taxonId and m.genome_external_db_release_id = $genomeId "
	. "and m.selected != 1 "
	. "and m.genome_dots_gene_id = gdt.genome_dots_gene_id "
	. "and m.similarity_dots_gene_id = gdt.similarity_dots_gene_id "
	. "and gdt.na_sequence_id = a.na_sequence_id and "
	. "m.genome_dots_gene_id = gdg.genome_dots_gene_id "
	. "group by m.genome_dots_gene_id, decode(greatest(gdg.max_intron, 14), 14, 0, 1), "
	. " gdg.gene_size, m.similarity_dots_gene_id "
	. "order by decode(greatest(gdg.max_intron, 14), 14, 0, 1) desc, gdg.gene_size desc, "
	. " max(a.contains_mrna) desc, sum(a.number_of_contained_sequences) desc";
    $self->log("finding one to many mappings, sql = $sql");
    my $sth = $dbh->prepareAndExecute($sql);

    # get an array of ordered mappings
    #     "better" aligned genes first: spliced first, larger gene size first
    #     "better" genes first: contain mrna first, contain more sequences first 
    my @ag_ords;
    my %ag_ids;
    my %id_maps;
    my $prev_ag;
    while (my ($agid, $gid) = $sth->fetchrow_array) {
	push @ag_ords, $agid if $agid != $prev_ag;
	$ag_ids{$agid} = undef;
	$prev_ag = $agid;
	
	if (!exists $id_maps{$agid}) { $id_maps{$agid} = []; }
	push @{ $id_maps{$agid} }, $gid;
    }

    $self->log("found ", scalar(keys %ag_ids), " genome_dots_gene_ids in this case");

    my %used_gids;
    my $count;
    foreach my $agid (@ag_ords) {
	my $gids = $id_maps{$agid};
	my $mapped = 0;
	foreach my $gid (@$gids) {
	    # exit inner loop if this agid is already mapped
	    if ($mapped) {
		delete $id_maps{$agid};
		last;
	    }

	    # next if this gene id is already taken
	    next if exists $used_gids{$gid};

	    # o.w. update table to indicate this mapping
	    $sql = "update $tmpTab set selected = 1 "
		. "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
		. "and genome_dots_gene_id = $agid and similarity_dots_gene_id = $gid";
	    $dbh->sqlexec($sql);
	    if (++$count % 500 == 0) { 
		$dbh->commit();
		$self->log("set $count 1 to 2+ entries");
	    }

	    # take the gene_id off from further use, and flag mapped
	    $used_gids{$gid} = "";
	    $mapped = 1;
	}
    }
    $self->log("Mapped $count in the 1:2+ case, " . scalar(keys %ag_ids) - $count . " remain unmapped");
}

sub setSelected {
    my ($dbh, $tmpTab, $taxonId, $genomeId, $agids) = @_;
    $agids =~ s/,$//;
    my $sql = "update $tmpTab set selected = 1 " .
	      "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	      . "and genome_dots_gene_id in ($agids)";
    $dbh->sqlexec($sql);
}

sub transferGeneId {
    my ($self, $dbh, $tmpTab, $gdgTab, $taxonId, $genomeId) = @_;

    my $sql = "select genome_dots_gene_id, similarity_dots_gene_id from $tmpTab "
            . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId and selected = 1";
    $self->log("select id for transfer from $tmpTab to $gdgTab: sql=$sql");
    my $sth = $dbh->prepareAndExecute($sql);

    my @id_pairs;
    while (my @a = $sth->fetchrow_array) { push @id_pairs, \@a; }
    $self->log("found " . scalar(@id_pairs) . " to transfer");

    my $c = 0;    
    foreach (@id_pairs) {
	my ($gdg_id, $sdg_id) = @$_;
        my $sql = "update $gdgTab set gene_id = $sdg_id "
            . "where taxon_id = $taxonId and genome_external_db_release_id = $genomeId "
	    . "and genome_dots_gene_id = $gdg_id";
        $dbh->sqlexec($sql);
	
	unless (++$c % 500) {
	    $dbh->commit();
	    $self->log("transfered $c gene ids");
	}
    }
    $dbh->commit();
    $self->log("total gene id transfered: $c");
}


1;
