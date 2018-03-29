package Snippa::Command::makeref;

use parent Snippa::Command;

use Snippa::Logger ':all';
use Snippa::Util ':all';
use Data::Dumper;
use File::Path qw(mkpath);
use File::Copy;

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;

  my $sc = $ap->add_parser('makeref', help=>'Make a shared reference folder', parents=>[$shared]);

  $sc->add_arg('--reference', '-r', type=>'Scalar', required=>1,
    help => "Reference genome",
  );
  $sc->add_arg('--outdir', '-o', type=>'Scalar', required=>1,
    help => "Output folder",
  );
  return $sc;
}

#----------------------------------------------------------------------

sub run {
  my($self, $arg) = @_;
  msg("running 'makeref'");
  my $ref = $arg->reference;
  my $dir = $arg->outdir;
  my $fmt = guess_format($ref);
  msg("Formatting '$fmt' reference: $ref");
  $fmt eq 'fasta' or err("Only 'fasta' format is supported");
  mkpath $dir;
  my $fasta = "$dir/ref.fa";
  copy($ref, $fasta);
  run_cmd("samtools faidx \Q$fasta\E");
  run_cmd("bwa index \Q$fasta\E");
#  run_cmd("bowtie2-build \Q$fasta\E \Q$fasta\E");
}

#----------------------------------------------------------------------

sub guess_format {
  my($fname) = @_;
  open my $fh, '<', $fname or err("Could not read '$fname'");
  my $line = <$fh>;
  return 'genbank' if $line =~ m/^LOCUS/;
  return 'fasta' if $line =~ m/^>\S/;
  return 'gff' if $line =~ m/^##gff/;
  return 'unknown';  
}

#----------------------------------------------------------------------

1;
