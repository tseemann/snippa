package Snippa::Command::report;

use parent Snippa::Command;

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;
  
  my $sc = $ap->add_parser('report', help=>'Generate a report from BAM + VCF', parents=>[$shared]);
  
  $sc->add_arg('--vcf', '-v', type=>'Scalar', required=>1,
    help => "VCF input",
  );
  $sc->add_arg('--bam', '-b', type=>'Scalar',
    help => "FASTQ input",
  );
  $sc->add_arg('--outdir', '-o', type=>'Scalar', required=>1,
    help => "Output folder",
  );
  
  return $sc;
}

#----------------------------------------------------------------------

1;

