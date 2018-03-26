package Snippa::Logger;

use base Exporter;
@EXPORT_OK = qw(msg wrn err);
%EXPORT_TAGS = ( 'all' =>  [ @EXPORT_OK ] );

#----------------------------------------------------------------------

our $quiet = 0;

#----------------------------------------------------------------------

sub quiet {
  my($self, $value) = @_;
  $quiet = $value;
  return $quiet;
}

#----------------------------------------------------------------------

sub msg {
  return if $quiet;
#  my($pkg,undef,$line) = caller(1);
#  print STDERR "[$pkg:$line] @_\n";
#  my $t = localtime;
  print STDERR "@_\n";
}
      
#----------------------------------------------------------------------

sub wrn {
  msg("WARNING", @_);
}

#----------------------------------------------------------------------

sub err {
  msg("ERROR", @_);
  exit(1);
}

#----------------------------------------------------------------------

1;

