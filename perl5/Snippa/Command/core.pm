package Snippa::Command::core;

use parent Snippa::Command;

#----------------------------------------------------------------------

sub add_command {
  my($self, $ap, $shared) = @_;
  
  my $sc = $ap->add_parser('core', help=>'Build core gene alignments from VCF', parents=>[$shared]);
  
  $sc->add_arg('--ambiguous', '-a', type=>'Bool',
    help => "Allow ambiguous bases",
  );
  $sc->add_arg('--gaps', '-g', type=>'Bool',
    help => "Allow missing bases / gaps",
  );
  
  return $sc;
}

#----------------------------------------------------------------------

1;

