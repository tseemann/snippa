#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

#----------------------------------------------------------------------
# globals

my $EXE = basename($0);
my $VERSION = "0.1.0";
my $AUTHOR = "Torsten Seemann";

#----------------------------------------------------------------------
# command line parameters

my $max = 5;
my $ref = '';
my $invert = 0;
my $debug = 0;

#----------------------------------------------------------------------
# getopts

GetOptions(
  "help"      => \&usage,
  "V|version" => \&version,
  "ref=s"     => \$ref,
  "max=i"     => \$max,
  "v|invert"  => \$invert,
  "ends"      => \$invert,
  "debug"     => \$debug,
) or usage(1);
             
$ref or err("Please supply reference genome with --ref");
$max >= 0 or err("Please supply --max >= 0");  
$ref .= ".fai" unless $ref =~ m/\.fai$/;
-r $ref or err("Can't see '$ref' index. Run 'samtools faidx $ref' ?"); 
!@ARGV and -t STDIN and err("Please provide or pipe a SAM file to $EXE");

#----------------------------------------------------------------------
# main
msg("$EXE $VERSION by $AUTHOR");

# get a hash of { seqname => length }
my $len = fai_to_dict($ref);
msg(Dumper($len)) if $debug;
msg("Found", scalar keys %$len, "sequences in $ref");

# read SAM
my $total=0;
my $removed=0;
my $kept=0;

while (my $line = <ARGV>) {
  $total++;
  my @sam = split m/\t/, $line;
  if ($sam[5] =~ /\d[SH]/) {
    my($HL, $SL, undef, $SR, $HR) 
      = ($sam[5] =~ m/ ^ (?:(\d+)H)? (?:(\d+)S)? (.*?) (?:(\d+)S)? (?:(\d+)H)? $/x);
    $HL ||= 0; $SL ||= 0; $SR ||= 0; $HR ||= 0;
    msg("## CHROM=$sam[2] POS=$sam[3] CIGAR=$sam[5] HL=$HL SL=$SL SR=$SR HR=$HR (max=$max)") if $debug;
    if ($HL+$SL > $max or $HR+$SR > $max) {
      $removed++;
      next;
    }
    print $line if $invert;
    $kept++;
  }
  print $line unless $invert;
}
msg("Read $total, removed $removed, allowed $kept, passed", $total-$removed);

#----------------------------------------------------------------------
sub fai_to_dict {
  my($fname) = @_;
  open my $fai, '<', $fname or err("Can't read FAI '$fname'");
  my $len = {};
  while (<$fai>) {
    my($name, $bp) = split m/\t/;
    $len->{$name} = $bp;
  }
  close $fai;
  return $len;
}

#----------------------------------------------------------------------
sub usage {
  my $errcode = @_;
  my $fh = $errcode ? \*STDERR : \*STDOUT;
  print $fh "Usage: $EXE --ref ref.fasta [--max=$max] [-v/--invert] < in.sam > out.sam\n";
  exit($errcode);
}

#----------------------------------------------------------------------
sub version {
  print "$EXE $VERSION\n";
  exit(0);
}

#----------------------------------------------------------------------
sub msg {
  print STDERR "[$EXE] @_\n";
}

#----------------------------------------------------------------------
sub err {
  msg("ERROR:", @_);
  exit(1);
}
