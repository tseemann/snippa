package Snippa::Command::snpeff;

use parent Snippa::Command;

#----------------------------------------------------------------------

use Snippa::Logger ':all';
use Snippa::Util 'run_cmd';
use File::Path qw(mkpath);

#----------------------------------------------------------------------

my %SNPEFF = (
  'snpeff' => 'snpEff --version',
   # http://samtools.github.io/bcftools/howtos/csq-calling.html
  'bcftools'  => 'bcftools csq -f {REF} -g {GFF} {VCF_IN} -p R -Ov -o {VCF_OUT}',
);

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;
  
  my $sc = $ap->add_parser('snpeff', help=>'Determine consequences of variants', parents=>[$shared]);
  
  $sc->add_arg('--method', '-m', choices_i => [ sort keys %SNPEFF ], reset=>1,
    help => "Variant caller to use", default => 'bcftools',
  );
  $sc->add_arg('--refdir', '-r', type=>'Scalar', required=>1,
    help => "Reference genome folder",
  );
  $sc->add_arg('--vcf', '-v', type=>'Scalar', required=>1,
    help => "VCF file to use",
  );
  $sc->add_arg('--gff', '-g', type=>'Scalar', required=>1,
    help => "Annotations of reference in GFF3 format",
  );
  $sc->add_arg('--outdir', '-o', type=>'Scalar', required=>1,
    help => "Output folder",
  );
  
  return $sc;
}

#----------------------------------------------------------------------
sub run {
  my($self, $arg) = @_;
  my $ref = $arg->refdir.'/ref.fa';
  my $dir = $arg->outdir;
  my $fname = 'snpeff';
  my $vcf = $arg->vcf;
  my $gff = $arg->gff;

  my $cmd = $SNPEFF{$arg->method};
  $cmd =~ s/{REF}/$ref/g;
  $cmd =~ s/{GFF}/$gff/g;
  $cmd =~ s/{VCF_IN}/$vcf/g;
  $cmd =~ s/{VCF_OUT}/$dir\/$fname.vcf/g;
  msg("CMD = $cmd");
  mkpath $dir;
  run_cmd($cmd);
}

#----------------------------------------------------------------------


#----------------------------------------------------------------------

1;

