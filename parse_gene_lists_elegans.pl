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
my $species = "c_elegans";
my $project = "PRJNA13758";
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();

#
my $gene_dir = $home . "release/$PRODUCTION_RELEASE/$species/gene_lists/";
my $local_file  = $gene_dir . "geneIDs.txt";

my $output_file = $gene_dir . "wormbase_gene_id_name.list";
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}

my @lines = read_file($local_file);


foreach my $line (@lines){

my ($id, $wb_gene_id, $gene_name, $alt_gene_name, $live) = split(/\,/, $line);

if ($live =~/Live/){
 my $name;
 if ($gene_name){
    $name = $gene_name;
 } elsif ($alt_gene_name) {
    $name = $alt_gene_name;
 }
 if (($name) and ($wb_gene_id)){
  my $term = "$wb_gene_id\t$name\n";
  write_file($output_file, {append => 1 }, $term);
 }
}

}
#
exit 0;
