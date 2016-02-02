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
my $go_function_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_function_sentences/";
my $go_process_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_process_sentences/";
my $go_component_directory= $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_component_sentences/";
my $tissue_expression_directory = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/tissue_expression/output_files/individual_gene_sentences/";
my $anatomy_ec_directory = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/anatomy_expression_cluster/output_files/individual_gene_sentences/";
my $molecule_reg_ec_directory = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/molecule_regulation_expression_cluster/output_files/individual_gene_sentences/";
my $gene_reg_ec_directory = $concise_descriptions . "release/$PRODUCTION_RELEASE/$species/gene_regulation_expression_cluster/output_files/individual_gene_sentences/";
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
my @component_files = glob("$go_component_directory/WBGene*");
my @function_files  = glob("$go_function_directory/WBGene*");
my @process_files   = glob("$go_process_directory/WBGene*");
my @anatomy_ec_files = glob("$anatomy_ec_directory/WBGene*");
my @molecule_reg_ec_files = glob("$molecule_reg_ec_directory/WBGene*");
my @gene_reg_ec_files = glob("$gene_reg_ec_directory/WBGene*");
#
# Create a list of genes
#
my @list = (@homology_files, @component_files, @function_files, @process_files, @tissue_expression_files, @anatomy_ec_files, @molecule_reg_ec_files, @gene_reg_ec_files);
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
my $process_size  =0;
my $function_size =0;
my $component_size=0;
my $tissue_size = 0;
my $gene_ec_size = 0;
my $mole_ec_size = 0;
my $anatomy_ec_size = 0;
#
foreach my $file (@sort_unique_list){
 write_file($summary, {append => 1 }, $dash_line);
 write_file($summary, {append => 1 }, $file);
 write_file($summary, {append => 1 }, "\n");
#
 my $homology_no_go_file  = $homology_directory  . $file;
 my $homology_go_from_elegans_file  = $homology_go_from_elegans_directory  . $file;
 my $go_process_file   = $go_process_directory  . $file;
 my $go_component_file = $go_component_directory  . $file;
 my $go_function_file  = $go_function_directory  . $file;
 my $gene_ec_file = $gene_reg_ec_directory . $file;
 my $mole_ec_file = $molecule_reg_ec_directory . $file;
 my $anatomy_ec_file = $anatomy_ec_directory . $file;
 my $tissue_file = $tissue_expression_directory . $file;
 my $homology ="";
 my $homology_go  = "";
 my $go_process   ="";
 my $go_component ="";
 my $go_function  ="";
 my $tissue = "";
 my $gene_ec = "";
 my $mole_ec = "";
 my $anatomy_ec = "";

 if (-e $tissue_file ){
   $tissue = read_file($tissue_file);
   chomp($tissue);
   $tissue_size++;
 }

 if (-e $gene_ec_file ){
   $gene_ec = read_file($gene_ec_file);
   chomp($gene_ec);
   $gene_ec_size++;
 }

 if (-e $mole_ec_file ){
   $mole_ec = read_file($mole_ec_file);
   chomp($mole_ec);
   $mole_ec_size++;
 }

 if (-e $anatomy_ec_file ){
   $anatomy_ec = read_file($anatomy_ec_file);
   chomp($anatomy_ec);
   $anatomy_ec_size++;
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
 if (-e $go_process_file) {
     $go_process = read_file($go_process_file);
     $process_size++;
     }
 if (-e $go_function_file) {
     $go_function = read_file($go_function_file);
     $function_size++;
     }
 if (-e $go_component_file) {
     $go_component = read_file($go_component_file);
     $component_size++;
     }
     $homology    =~ s/^\s+//;
     $homology    =~ s/\s+$//;
     $go_function =~ s/^\s+//;
     $go_function =~ s/\s+$//;
     $go_component=~ s/^\s+//;
     $go_component=~ s/\s+$//;
     $go_process  =~ s/^\s+//;
     $go_process  =~ s/\s+$//;
     $anatomy_ec  =~ s/^\s+//;
     $anatomy_ec  =~ s/\s+$//;
     $gene_ec  =~ s/^\s+//;
     $gene_ec  =~ s/\s+$//;
     $mole_ec  =~ s/^\s+//;
     $mole_ec  =~ s/\s+$//;

     if ($go_process){
      chomp($go_process);
      chop($go_process);
      $go_process .= "\;";
     }
     if ($go_function){
      chomp($go_function);
      chop($go_function);
      $go_function .= "\;";
     }
     if ($go_component){
      chomp($go_component);
      chop($go_component);
      $go_component .= "\;";
     }

 my $output;
    if (($go_process) or ($go_function) or ($go_component) or ($tissue)){
     $output = $homology . " " . $go_process  . " " . $go_function . " " . $tissue . " " . $go_component;
    } elsif ($homology) {
     $output = $homology . " " . $gene_ec . " " . $mole_ec . " " . $anatomy_ec;
    } else {
     $output = ucfirst($gene_ec) . " " . $mole_ec . " " . $anatomy_ec;
    }
    $output =~ s/  / /g;
    $output =~ s/\n//g;
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
my $number_homology  = "Number of automated descriptions with orthology\: " . $homology_size . "\n";
my $number_process   = "Number of automated descriptions with GO process information\: " . $process_size . "\n";
my $number_function  = "Number of automated descriptions with GO molecular function information\: " . $function_size . "\n";
my $number_component = "Number of automated descriptions with GO cellular component information\: " . $component_size . "\n";
my $number_tissue = "Number of automated descriptions with tissue expression information\: " . $tissue_size . "\n";
my $number_gene = "Number of automated descriptions with gene regulation expression cluster information\: " . $gene_ec_size . "\n";
my $number_mole = "Number of automated descriptions with molecule regulation expression cluster information\: " . $mole_ec_size . "\n";
my $number_anatomy = "Number of automated descriptions with anatomy expression cluster information\: " . $anatomy_ec_size . "\n";

my $report_output = $total . $number_homology . $number_process . $number_function . $number_component . $number_tissue;
   $report_output .=  $number_gene . $number_mole . $number_anatomy;
 write_file($report, $report_output); 
#
}
exit 0;
