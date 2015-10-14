#!/usr/bin/env perl 
use Carp;
use LWP::Simple;
use File::Slurp;
use ConciseDescriptions;
use strict;
use warnings;
#
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions
#
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my @args;
my $status;
#
my $species = "c_elegans";

my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $tissue_dir = $home . "release/$PRODUCTION_RELEASE/$species/tissue_expression/input_files/";
if (-e $tissue_dir){
   print "$tissue_dir exists\n"; 
} else {
   mkdir $tissue_dir or die "could not create $tissue_dir";
}
my $local_file  = $tissue_dir . "anatomy_association.$RELEASE.wb";

if (-e $local_file){
   @args = ("/bin/rm", "-f", $local_file);
   system(@args) == 0 or die("could not delete file $local_file\n");
}
#
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/ONTOLOGY/anatomy_association.$RELEASE.wb";
#
getstore($url, $local_file);
#
exit 0;
