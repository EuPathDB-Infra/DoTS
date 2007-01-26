package DoTS::DotsBuild::Plugin::UpdateAssemblySourceId;

@ISA = qw(GUS::PluginMgr::Plugin);
use strict;


use GUS::PluginMgr::Plugin;
use GUS::Model::DoTS::Assembly;


sub new {
  my ($class) = @_;

  my $self = {};
  bless($self,$class);

my $purpose = <<PURPOSE;
Update dots.assembly.source_id from a mapping file.
PURPOSE

my $purposeBrief = <<PURPOSE_BRIEF;
Updating dots.assembly.source_id from file.
PURPOSE_BRIEF

my $notes = <<NOTES;
NOTES

my $tablesAffected = <<TABLES_AFFECTED;
    [ ['DoTS::Assembly', ''] ];
TABLES_AFFECTED

my $tablesDependedOn = <<TABLES_DEPENDED_ON;
TABLES_DEPENDED_ON

my $howToRestart = <<RESTART;
No additional restart procedure.
RESTART

my $failureCases = <<FAIL_CASES;
Regex incorrect or file incomplete or inconsistent.
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
   fileArg({name => 'mappingFile',
            descr => 'file with assembly na_sequence_id to source_id mapping',
            constraintFunc=> undef,
            reqd  => 1,
            isList => 0,
            mustExist => 1,
            format => '',}),
   stringArg({name => 'assemblyRegex',
	      descr => 'regex to obtain na_sequence_id from line of mapping file',
	      constraintFunc=> undef,
	      reqd  => 1,
	      isList => 0
	     }),
   stringArg({name  => 'sourceIdRegex',
	      descr => 'regex to obtain source_id from line of mapping file',
	      constraintFunc=> undef,
	      reqd => 1,
              isList => 0
              }),
   stringArg({name  => 'unassigned',
	      descr => 'value in place of source_id when none is mapped',
              constraintFunc=> undef,
	      reqd => 0,
	      isList => 0
	      }),
   stringArg({name  => 'prefix',
	      descr => 'string prepended to the source ids in the mapping file, e.g. DT.',
	      constraintFunc=> undef,
	      reqd => 0,
	      isList => 0
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

    my $self = shift;

    my $count = $self->processMappingFile();

    return "$count dots.assembly.source_ids updated";

}

sub processMappingFile {
  my $self   = shift;
  $self->log("processing mapping file");

  my $count;

  my $file = $self->getArg('mappingFile');

  open (FILE, $file) || die "Can't open file $file for reading\n";

  my $sourceIdRegex = $self->getArg('sourceIdRegex');

  my $assemblyRegex = $self->getArg('assemblyRegex');

  my $sourceId;
  my $naSequenceId;

  while(<FILE>){
    if (/$sourceIdRegex/ && $1) {
      $sourceId = $1;
    }
    else {
      my $forgotParens = ($sourceIdRegex !~ /\(/) ? "(Forgot parens?)" : "";
      $self->userError("Unable to parse source_id from $_ using regex '$sourceIdRegex' $forgotParens");
    }

    if (/$assemblyRegex/ && $1) {
      $naSequenceId = $1;
    } else {
      my $forgotParens = ($assemblyRegex !~ /\(/)? "(Forgot parens?)" : "";
      $self->userError("Unable to parse source_id from $_ using regex '$assemblyRegex' $forgotParens");
    }

    $count += $self->updateAssembly($sourceId,$naSequenceId);

    $self->log("$count assembly rows updated") if $count % 10000 == 0;
  }

  return $count;
}



sub updateAssembly {
  my ($self,$sourceId,$naSequenceId) = @_;

  my $assemblyRow  = GUS::Model::DoTS::Assembly->new({'na_sequence_id'=>$naSequenceId});

  $assemblyRow->retrieveFromDB();

  my $unassigned = $self->getArg('unassigned') ? $self->getArg('unassigned') : "";

  $sourceId = $naSequenceId if $sourceId eq $unassigned;

  my $prefix = $self->getArg('prefix') ? $self->getArg('prefix') : "";

  $sourceId = "${prefix}$sourceId";

  if ($assemblyRow->get('source_id') ne $sourceId) {
    $assemblyRow->set('source_id',$sourceId);
  }

  my $submit = $assemblyRow->submit();

  if ($submit == 1) {
    $self->log("$naSequenceId updated to source_id = $sourceId");
  }
  else {
    $self->log("$naSequenceId FAILED to update to source_id = $sourceId");
  }

  $self->undefPointerCache();

  return $submit;
}

1;

