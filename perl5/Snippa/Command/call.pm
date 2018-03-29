package Snippa::Command::call;

use parent Snippa::Command;

#----------------------------------------------------------------------

use Snippa::Logger ':all';
use Snippa::Util 'run_cmd';
use File::Path qw(mkpath);

#----------------------------------------------------------------------

my %CALLER = (
  'freebayes' => 'freebayes -f {REF} -b {BAM} --vcf {VCF} -P 0.5 --no-partial-observations --min-repeat-entropy -0',
  'bcftools'  => 'bcftools mpileup -f {REF} {BAM} | bcftools call -mv -Ov -o {VCF}',
#  'varscan'   => 'samtools mpileup {BAM} | varscan mpileup2snp /dev/stdin --variants --output-vcf > {VCF}',
  'pilon'     => 'pilon --threads {CPUS} --genome {REF} --bam {BAM} --vcf --outdir {VCF_DIR} --output {VCF_NAME} --variant --diploid --minmq 30 --mindepth 0.25',
);

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;
  
  my $sc = $ap->add_parser('call', help=>'Call SNPs from BAM', parents=>[$shared]);
  
  $sc->add_arg('--method', '-m', choices_i => [ sort keys %CALLER ], reset=>1,
    help => "Variant caller to use", default => 'freebayes',
  );
  $sc->add_arg('--refdir', '-r', type=>'Scalar', required=>1,
    help => "Reference genome folder",
  );
  $sc->add_arg('--bam', '-b', type=>'Scalar',
    help => "BAM file of aligned reads",
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
  my $bam = $arg->bam;
  my $dir = $arg->outdir;
  my $fname = 'snps';

  my $cmd = $CALLER{$arg->method};
  $cmd =~ s/{REF}/$ref/g;
  $cmd =~ s/{BAM}/$bam/g;
  $cmd =~ s/{VCF}/$dir\/$fname.vcf/g;
  $cmd =~ s/{VCF_DIR}/$dir/g;
  $cmd =~ s/{VCF_NAME}/$fname/g;
  msg("CMD = $cmd");

  mkpath $dir;
  run_cmd($cmd);
}

#----------------------------------------------------------------------

1;

