#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../perl5";
use Getopt::ArgParse;
use Data::Dumper;
use Module::Load;
use Snippa::Logger ':all';

our $TOOLNAME = 'snippa';
our $VERSION = '0.1';
my @COMMAND = (qw(makeref align call report core version));

my $ap = Getopt::ArgParse->new_parser(
  prog        => $TOOLNAME,
  description => 'Modular variant calling for bacterial sequences',
  epilog      => "Written by Torsten Seemann\nAvailable from http://github.com/tseemann/snippa",
);
$ap->add_arg('--version', '-V', type=>'Bool',
  help => "Print version and exit",
);

my $shared = Getopt::ArgParse->new_parser();
$shared->add_arg('--quiet', '-q', type=>'Bool',
  help => "Don't output anything to screen",
);
$shared->add_arg('--threads', '-t', type=>'Scalar', default=>1,
  help => "Number of CPU threads to use",
);
#$shared->add_arg('--outdir', '-o', type=>'Scalar', required=>0,
#  help => "Output folder for results",
#);

$ap->add_subparsers(title => 'subcommands');

my %parser;
for my $cmd (@COMMAND) {
  my $mod = "Snippa::Command::$cmd";
  load $mod;
  $parser{$cmd} = $mod->add_command($ap, $shared);
}
#print Dumper($ap);

my $arg = $ap->parse_args(); # uses @ARGV
print Dumper($arg);

my $cc = $arg->current_command || 'help';
wrn("current_command = $cc");

my $mod = "Snippa::Command::$cc";
$mod->run($arg);


