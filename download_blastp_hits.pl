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
my $species = $ARGV[0];
my $project = $ARGV[1];
my $home = $html . "concise_descriptions/release/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $orthology_dir = $home . "$PRODUCTION_RELEASE/$species/orthology/input_files/";
my $local_file = $orthology_dir . $species . ".$RELEASE.best_blastp_hits.txt.gz";
if (-e $local_file){
   @args = ("/bin/rm", "-f", $local_file);
   system(@args) == 0 or die("could not delete file $local_file\n");
}
my $unzipped_file;
#
my $output_file = $orthology_dir . $species . ".best_blastp_hits.txt";
if (-e $output_file){
   @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
#
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/species/$species/$project/$species.$project.$RELEASE.best_blastp_hits.txt.gz";
#
getstore($url, $local_file);
if (-e $local_file){
    $unzipped_file =$local_file;
    $unzipped_file =~s/\.gz//g;
}
if (-e $unzipped_file){
   @args = ("/bin/rm", "-f", $unzipped_file);
   system(@args) == 0 or die("could not delete file $unzipped_file\n");
}
 @args = ("/bin/gunzip", $local_file);
$status = system(@args);
#
if (-e $unzipped_file){
 @args = ("/bin/ln","-s",$unzipped_file,$output_file);
 $status = system(@args);
}
#
exit 0;
