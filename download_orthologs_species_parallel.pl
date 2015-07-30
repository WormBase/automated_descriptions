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
my $home = $html . "concise_descriptions/";
#
my $RELEASE = ConciseDescriptions::get_release(); 
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $AND = "AND";
my $species_project = $ARGV[0];
 chomp($species_project);
my ($species, $project) = split(/$AND/, $species_project);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
#
my $orthology_dir = $home . "release/$PRODUCTION_RELEASE/$species/orthology/input_files/";
#
my $local_file = $orthology_dir . "$species.$project.$RELEASE.orthologs.txt.gz";
if (-e $local_file){
   @args = ("/bin/rm", "-f", $local_file);
   system(@args) == 0 or die("could not delete file $local_file\n");
}
my $unzipped_file;
#
my $output_file = $orthology_dir . "$species.orthologs.txt";
if (-e $output_file){
   @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
#
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/species/$species/$project/annotation/$species.$project.$RELEASE.orthologs.txt.gz";
my $url2 = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/species/$species/$project/annotation/$species.$project.$RELEASE.orthologs.txt";
#
if (head($url)){
 print "$url exists\n";
 LWP::Simple::getstore($url, $local_file);
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
} elsif (head($url2)){
  print "$url2 exists\n";
  LWP::Simple::getstore($url2, $output_file);
  if (-e $output_file){
   print "file copied to $output_file\n";
  } else {
   print "$output_file not copied\n";
  }
}
#
exit 0;
