#!/usr/bin/env perl
use ConciseDescriptions;
use File::Slurp; 
use List::MoreUtils qw(uniq);
use strict;
use warnings;
#
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions
#
#
my $AND = "AND";
my $species_project_name_prefix = $ARGV[0];
chomp($species_project_name_prefix);
my ($species, $project, $species_name, $species_prefix) = split(/$AND/, $species_project_name_prefix);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $species_name =~ s/^\s+//;
 $species_name =~ s/\s+$//;
 $species_prefix =~ s/\s+$//;
 $species_prefix =~ s/^\s+//;

# The path for the output file is $functions_dir; 
# the path for the individual gene descriptions is $individual_path.
#
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my $html = ConciseDescriptions::get_html_dir();
my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
# 
#
if ($species){
 my $concise_descriptions = $html . "concise_descriptions/";
 my $double_space = "  ";
 my $space = " ";
 my $dash_line="";
 my $semantic_concise_descriptions = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/descriptions/";
 my $output_path = $semantic_concise_descriptions;
 my $summary = $output_path . "WBGenes_concise_descriptions.txt";
 my $report  = $semantic_concise_descriptions . "concise_descriptions_report.txt";
#
if (-e $output_path){
   print "$output_path exists\n"; 
} else {
   mkdir $output_path or die "could not create $output_path";
}
my $individual_path = $output_path . "individual_gene_descriptions/";
if (-e $individual_path){
   print "$individual_path exists\n"; 
} else {
   mkdir $individual_path or die "could not create $individual_path";
}
#
# if the output file exists delete it.
if (-e $summary){
   my @args = ("/bin/rm", "-f", $summary);
   system(@args) == 0 or die("could not delete file $summary\n");
}
if (-e $report){
   my @args = ("/bin/rm", "-f", $report);
   system(@args) == 0 or die("could not delete file $report\n");
} 
#
my $homology_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/orthology/output_files/individual_gene_sentences/";
my $go_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_gene_sentences/";
my $tissue_expression_directory = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/tissue_expression/output_files/individual_gene_sentences/";
my $homology_go_from_elegans_directory="";
my @homology_go_from_elegans_files=();
my @tissue_expression_files = ();
if ($species !~/elegans/){
 $homology_go_from_elegans_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/orthology/output_files/individual_gene_sentences_from_GO_elegans/";
 @homology_go_from_elegans_files= glob("$homology_go_from_elegans_directory/WBGene*");
} else {
 @tissue_expression_files = glob("$tissue_expression_directory/WBGene*");
}

#
# if the output files for individual gene sentences exist delete them.
   my @individual_files = glob("$individual_path/WBGene*");
      foreach my $individual (@individual_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
#
my @homology_files  = glob("$homology_directory/WBGene*");
my @go_files  = glob("$go_directory/WBGene*");
#
# Create a list of genes
#
my @list = (@homology_files, @go_files, @tissue_expression_files);
#
 my @unsorted_files = ();
 foreach my $item (@list){
     my @items = split("/",$item);
     my $size_of_array = @items;
     my $file = $items[$size_of_array-1];
     push(@unsorted_files, $file);
 }
 my @sort_unique_list = uniq(sort(@unsorted_files));
#
# Keep a counter for each data 
#
my $homology_size =0;
my $go_size  =0;
my $tissue_size = 0;
#
foreach my $file (@sort_unique_list){
 write_file($summary, {append => 1 }, $dash_line);
 write_file($summary, {append => 1 }, $file);
 write_file($summary, {append => 1 }, "\n");
#
 my $homology_no_go_file  = $homology_directory  . $file;
 my $homology_go_from_elegans_file  = $homology_go_from_elegans_directory  . $file;
 my $go_file  = $go_directory  . $file;
 my $tissue_file = $tissue_expression_directory . $file;
 my $homology ="";
 my $homology_go = "";
 my $go ="";
 my $tissue = "";

 if (-e $tissue_file ){
   $tissue = read_file($tissue_file);
   $tissue_size++;
 }

  if (-e $homology_no_go_file) {
     $homology = read_file($homology_no_go_file);
      $homology =~ s/^\s+//;
      $homology =~ s/\s+$//;
   if (-e $homology_go_from_elegans_file) {
      $homology_go = read_file($homology_go_from_elegans_file);
      $homology_go =~ s/^\s+//;
      $homology_go =~ s/\s+$//;
      $homology .= " " . $homology_go;
     }
     $homology_size++;
     }
 if (-e $go_file) {
     $go = read_file($go_file);
     $go_size++;
     }
     $homology    =~ s/^\s+//;
     $homology    =~ s/\s+$//;
     $go   =~ s/^\s+//;
     $go   =~ s/\s+$//;

     if ($go){
      chomp($go);
      chop($go);
      $go .= "\;";
     }
 my $output;
    $output = $homology . " " . $go . " " . $tissue;
    $output =~ s/  / /g;
    $output =~ s/$double_space/$space/g;

    $output =~ s/^\s+//;
    $output =~ s/\s+$//;

my $last = chop($output);
if ($last=~/\;/){
    $output =~ s/\;$//;
    $output =~ s/^\s+//;
    $output =~ s/\s+$//;
    $output .= "\.\n";
    $output =~ s/^\s+//;
    $output =~ s/\s+$//;
} else{
    $output =~ s/\n$//;
    $output =~ s/^\s+//;
    $output =~ s/\s+$//;
    $output =~ s/\;$//;
    $output =~ s/^\s+//;
    $output =~ s/\s+$//;
    $output .= "\.\n";
    $output =~ s/^\s+//;
    $output =~ s/\s+$//;
}
 
 my $outfile = $individual_path . $file;

 if ($species =~/ratti/){
  if ($output =~/SRAE/){
   my $this_gene = "this gene";
   $output =~s/SRAE\_\S+/$this_gene/g;
  }
 }

 if ($species !~/elegans/){ 
  write_file($outfile, ucfirst($output));
  write_file($summary, {append => 1 }, ucfirst($output));
  write_file($summary, {append => 1 }, "\n\n\n");
 } else {
  write_file($outfile, $output);
  write_file($summary, {append => 1 }, $output);
  write_file($summary, {append => 1 }, "\n\n\n");
 }

}
my @concise_description_files = glob("$individual_path/WBGene*");
my $sum=0;
foreach my $description (@concise_description_files){
 next if ($description=~/WBGene00000000/);
 $sum++;
}
my $total = "Total number of automated descriptions\: " . $sum . "\n";
my $number_homology = "Number of automated descriptions with orthology\: " . $homology_size . "\n";
my $number_go = "Number of automated descriptions with GO information\: " . $go_size . "\n";
my $number_tissue = "Number of automated descriptions with tissue expression information\: " . $tissue_size . "\n";
my $report_output = $total . $number_homology . $number_go . $number_tissue;
 write_file($report, $report_output); 
#
}
exit 0;
