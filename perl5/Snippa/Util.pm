package Snippa::Util;

use base Exporter;
@EXPORT_OK = qw(run_cmd);
%EXPORT_TAGS = ( 'all' =>  [ @EXPORT_OK ] );

use strict;
use Snippa::Logger ':all';

#----------------------------------------------------------------------

sub run_cmd {
  msg("Running: @_");
  system(@_)==0 or err("Error: $!");
}

#----------------------------------------------------------------------

1;

