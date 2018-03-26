#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../perl5";
use Snippa::Logger ':all';
use Getopt::ArgParse;

msg("Hello", $ENV{USER});
wrn("It is ".localtime);
err("Nothing to do!");

