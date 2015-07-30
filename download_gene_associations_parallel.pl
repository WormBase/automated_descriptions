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
my $download_elegans_script = "./download_gene_associations_elegans.pl";
my $AND = "AND";
my @args;
my $status;
#
my $species_project = $ARGV[0];
my ($species, $project) = split(/$AND/, $species_project);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;

if ($species !~/elegans/){

my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();

#
my $go_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene_ontology/input_files/";
if (-e $go_dir){
   print "$go_dir exists\n"; 
} else {
   mkdir $go_dir or die "could not create $go_dir";
}
my $local_file  = $go_dir . "gene_association.$RELEASE.wb.$species";
if (-e $local_file){
   @args = ("/bin/rm", "-f", $local_file);
   system(@args) == 0 or die("could not delete file $local_file\n");
}
#
my $url = "ftp://ftp.sanger.ac.uk/pub/wormbase/releases/$RELEASE/ONTOLOGY/gene_association.$RELEASE.wb.$species";
#
getstore($url, $local_file);
#
} else{
 @args=($download_elegans_script, $species, $project);
 $status=system(@args);
}
exit 0;
