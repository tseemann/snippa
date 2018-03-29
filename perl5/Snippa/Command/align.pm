package Snippa::Command::align;

use parent Snippa::Command;

use Snippa::Logger ':all';
use Snippa::Util 'run_cmd';
use File::Path qw(mkpath);

#----------------------------------------------------------------------

my %ALIGNER = (
  'minimap2' => 'minimap2 -a -x sr -t {CPUS} {REF} {R1} {R2}',
  'bwamem'   => 'bwa mem -Y -t {CPUS} {REF} {R1} {R2}',
  'bowtie2'  => 'bowtie2 -X 1000 --end-to-end -p {CPUS} -x {REF} -1 {R1} -2 {R2}', 
);

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;

  my $sc = $ap->add_parser('align', help=>'Align reads to reference to make BAM', parents=>[$shared]);

  $sc->add_arg('--method', '-m', choices_i => [ sort keys %ALIGNER ], reset=>1,
    help => "Align reads using this", default=>'minimap2',
  );
  $sc->add_arg('--refdir', '-r', type=>'Scalar', required=>1,
    help => "Reference genome folder",
  );
  $sc->add_arg('--R1', '-1', type=>'Scalar', required=>1,
    help => "FASTQ for R1",
  );
  $sc->add_arg('--R2', '-2', type=>'Scalar', required=>1,
    help => "FASTQ for R2",
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
  my $cpus = $arg->threads;
  my $R1 = $arg->R1;
  my $R2 = $arg->R2;
  my $dir = $arg->outdir;

  my $cmd = $ALIGNER{$arg->method};
  $cmd =~ s/{REF}/$ref/g;
  $cmd =~ s/{R1}/$R1/g;
  $cmd =~ s/{R2}/$R2/g;
  $cmd =~ s/{CPUS}/$cpus/g;
  msg("CMD = $cmd");

  mkpath $dir;
  run_cmd(
    $cmd
    . " | samtools sort -n"
    . " | samtools fixmate -m - -"
    . " | samtools sort"
    . " | samtools markdup -r -s - -"
    . " > \Q$dir/reads.bam\E"
    . " && samtools index \Q$dir/reads.bam\E"
  );
    
}

#----------------------------------------------------------------------

1;
