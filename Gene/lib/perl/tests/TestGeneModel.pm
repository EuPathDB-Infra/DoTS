package DoTS::Gene::tests::TestGeneModel;

use strict;

use base qw(Test::Unit::TestCase);
use DoTS::Gene::GeneModel;

use constant DEBUG => 0;

sub new {
    my $self = shift()->SUPER::new(@_);
    # no need to do anything
    return $self;
}

sub set_up {
    my $self = shift;
    my (@t_inputs, @t_answers);
    push @t_inputs, { id => 999,
		      strand => '-',
		      coords => [[25,35], [10,15], [5,20], [40,50], [30,45]],
                      debug => DEBUG,
		    };
    push @t_answers, { model_string => '5-20 .. 25-50',
		       gene_size => 40,
		       genomic_start => 5,
		       genomic_end => 50,
		       min_exon_size => 15,
		       max_exon_size => 25,
		       min_intron_size => 5,
		       max_intron_size => 5,
		     };

    my @t_genemodels;
    foreach (@t_inputs) {
      my $genemodel =  GeneModel->new($_->{id}, $_->{strand}, $_->{coords}, $_->{debug});
      push @t_genemodels, $genemodel;
    }

    $self->{t_inputs} = \@t_inputs;
    $self->{t_answers} = \@t_answers;
    $self->{t_genemodels} = \@t_genemodels;
}

sub tear_down {
    my $self = shift;
    if (DEBUG) {
	print " no need to do anything ";
    }
}

# test subs

sub test_get_model_string {
    my $self = shift;

    my @t_genemodels = @{ $self->{t_genemodels} };
    my @t_answers = @{ $self->{t_answers} };

    foreach (@t_genemodels) {
      my $actual = $_->getModelString;
      my $ans = pop @t_answers;
      my $expect = $ans->{model_string};
      $self->assert($actual eq $expect, "actual=$actual, expect=$expect");
    }
}

sub test_get_gene_size {
    my $self = shift;

    my @t_genemodels = @{ $self->{t_genemodels} };
    my @t_answers = @{ $self->{t_answers} };

    foreach (@t_genemodels) {
      my $actual = $_->getGeneSize;
      my $ans = pop @t_answers;
      my $expect = $ans->{gene_size};
      $self->assert($actual == $expect, "actual=$actual, expect=$expect");
    }
}

1;
