#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;
use ConciseDescriptions;
#
my $species_release = $ARGV[0];

my $AND = "AND";
my ($species, $release) = split(/$AND/, $species_release);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $release =~s/^\s+//;
 $release =~s/\s+$//;

my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/";

my $release_directory = $home. $release . "\/";
my $species_directory =  $release_directory . $species . "\/";
my $descriptions_directory = $species_directory . "descriptions\/";
my $individual_descriptions_directory = $species_directory . "descriptions\/individual\_gene_descriptions\/";
my $gene_lists_directory = $species_directory . "gene\_lists\/";
#
my $GO_directory = $species_directory . "gene\_ontology\/";
my $GO_input_directory = $species_directory . "gene\_ontology\/input\_files\/";
my $GO_output_directory = $species_directory . "gene\_ontology\/output\_files\/";
my $GO_individual_output_directory = $species_directory . "gene\_ontology\/output\_files\/individual\_gene\_sentences\/";
#
my $orthology_directory = $species_directory . "orthology\/";
my $orthology_input_directory = $species_directory  . "orthology\/input\_files\/";
my $orthology_output_directory = $species_directory . "orthology\/output\_files\/";
my $orthology_individual_output_directory = $species_directory . "orthology\/output\_files\/individual\_gene\_sentences\/";
my $orthology_individual_output_directory_from_go_elegans = $species_directory . "orthology\/output\_files\/individual\_gene\_sentences\_from\_GO\_elegans\/";
my $tissue_expression_elegans = $species_directory . "tissue\_expression";
my $tissue_expression_elegans_input = $species_directory . "tissue\_expression\/input\_files";
my $tissue_expression_elegans_output = $species_directory . "tissue\_expression\/output\_files";
my $tissue_expression_elegans_output_directory = $species_directory . "tissue\_expression\/output\_files\/individual\_gene\_sentences\/";

if (($species) and ($release)) {

 make_dir($release_directory);
 make_dir($species_directory);
 make_dir($descriptions_directory);
 make_dir($individual_descriptions_directory);
 make_dir($gene_lists_directory);
 make_dir($GO_directory);
 make_dir($GO_input_directory);
 make_dir($GO_output_directory);
 make_dir($GO_individual_output_directory);
 make_dir($orthology_directory);
 make_dir($orthology_input_directory);
 make_dir($orthology_output_directory);
 make_dir($orthology_individual_output_directory);

 if ($species !~/elegans/){
   make_dir($orthology_individual_output_directory_from_go_elegans);
 } elsif  ($species =~/elegans/){
   make_dir($tissue_expression_elegans);
   make_dir($tissue_expression_elegans_input);
   make_dir($tissue_expression_elegans_output);
   make_dir($tissue_expression_elegans_output_directory);
 } 
}
exit 0;

sub make_dir {
   my $directory = shift;
   if (-e $directory){
    print "directory $directory created already\n";
   }
   else { 
    print "creating directory $directory\n";   
    unless(mkdir $directory) {
     die "Unable to create $directory\n";
    }
   }
}
