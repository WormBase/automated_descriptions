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
my $html = ConciseDescriptions::get_html_dir();
#
my $species = "c_elegans";
my $home = $html . "concise_descriptions/release/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $tissue_dir = $home . "$PRODUCTION_RELEASE/$species/tissue_expression/input_files/";
my $local_file = $tissue_dir . "WBbt.obo";
if (-e $local_file){
   @args = ("/bin/rm", "-f", $local_file);
   system(@args) == 0 or die("could not delete file $local_file\n");
}
my $unzipped_file;
#
my $output_file = $tissue_dir . "WBbt.obo";
if (-e $output_file){
   @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
#
#my $url = "https://raw.githubusercontent.com/raymond91125/Wao/master/WBbt.obo";
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/ONTOLOGY/anatomy_ontology.$RELEASE.obo";
#
getstore($url, $local_file);
#
exit 0;
