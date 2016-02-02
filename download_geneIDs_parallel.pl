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
# c_elegans	PRJNA13758
#my $species = "c_elegans";
#my $project = "PRJNA13758";
my $AND = "AND";
my $species_project = $ARGV[0];
 chomp($species_project);
my ($species, $project) = split(/$AND/, $species_project);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
#
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $gene_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene_lists/";

my $local_file  = $gene_dir . "geneIDs.txt.gz";

if (-e $local_file){
   my $file = $local_file;
      $file =~s/\.gz//gi;
   @args = ("/bin/rm", "-f", $file);
   system(@args) == 0 or die("could not delete file $file\n");
}
#
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/species/$species/$project/annotation/$species.$project.$RELEASE.geneIDs.txt.gz";
#
getstore($url, $local_file);
#
if (-e $local_file) {
 @args = ("/bin/gunzip", $local_file);
 system(@args) == 0 or die("could not unzip file $local_file\n");
}
#
exit 0;
