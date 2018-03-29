package Snippa::Command::version;

use parent Snippa::Command;

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;
  
  my $sc = $ap->add_parser('version', help=>'Print version and exit');
  
  return $sc;
}

#----------------------------------------------------------------------

sub run {
  my($self, $arg) = @_;
  print "$main::TOOLNAME $main::VERSION\n";
  exit(0);
}

#----------------------------------------------------------------------

1;

