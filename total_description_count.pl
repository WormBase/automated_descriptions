#!/usr/bin/env perl
use ConciseDescriptions;
use File::Slurp; 
use List::MoreUtils qw(uniq);
use LWP::Simple;
use LWP::UserAgent;
use strict;
use warnings;
#
my $html = ConciseDescriptions::get_html_dir();
my $PRODUCTION_RELEASE =  $ARGV[0];
my $home = $html . "concise_descriptions/";
my $species_file = $home . "species.txt";
my $summary = $home . "release/$PRODUCTION_RELEASE/__number_of_concise_descriptions.txt";
my $manual_descriptions =  $home .  "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/sort.curated_genes.txt";
# Remove previous incarnation of summary
if (-e $summary){
  my @args = ("rm", "-f", $summary);
  system(@args) == 0 or die("could not delete file $summary\n");
}
#
my @lines = read_file($species_file);
#
my @species_array = ();
my %species_name_hash = ();
foreach my $line (@lines){
 chomp($line);
 my ($species, $project, $name, $prefix) = split(/\t/,$line);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $name    =~s/^\s+//;
 $name    =~s/\s+$//;
 $name    =~s/\s/\_/g;
 $prefix  =~s/^\s+//;
 $prefix  =~s/\s+$//;
 $species_name_hash{$species} = $name;
 push(@species_array, $species);
}
#
my $kount = 0;
my $single_breakln = "\n";
my $breakln = "\n\n";

my @manual = read_file($manual_descriptions);
my $count_manual_total=0;
my $count_automated_only_total=0;
foreach my $species (@species_array){
 my $individual_gene_descriptions_directory = "$home/release/$PRODUCTION_RELEASE/$species/descriptions/individual_gene_descriptions/";
 my $orthology_directory ="$home/release/$PRODUCTION_RELEASE/$species/orthology/output_files/individual_gene_sentences/";
 my $GO_function_directory ="$home/release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_function_sentences/";
 my $GO_process_directory ="$home/release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_process_sentences/";
 my $GO_component_directory ="$home/release/$PRODUCTION_RELEASE/$species/gene_ontology/output_files/individual_component_sentences/";
 my $tissue_directory = "$home/release/$PRODUCTION_RELEASE/$species/tissue_expression/output_files/individual_gene_sentences/";
 my $anatomy_ec_directory = "$home/release/$PRODUCTION_RELEASE/$species/anatomy_expression_cluster/output_files/individual_gene_sentences/";
 my $mole_ec_directory = "$home/release/$PRODUCTION_RELEASE/$species/molecule_regulation_expression_cluster/output_files/individual_gene_sentences/";
 my $gene_ec_directory = "$home/release/$PRODUCTION_RELEASE/$species/gene_regulation_expression_cluster/output_files/individual_gene_sentences/";

 if (-e $individual_gene_descriptions_directory ) {
     my @descriptions = read_dir( $individual_gene_descriptions_directory, prefix => 1 );
     my @process   = read_dir($GO_process_directory, prefix => 1 );
     my @function  = read_dir($GO_function_directory, prefix => 1 );
     my @component = read_dir($GO_component_directory, prefix => 1 );
     my @orthology = read_dir($orthology_directory, prefix => 1 );
     my @gene_ec = read_dir($gene_ec_directory, prefix => 1 );
     my @mole_ec = read_dir($mole_ec_directory, prefix => 1 );
     my @anatomy_ec = read_dir($anatomy_ec_directory, prefix => 1 );

     my @tissue=();
     if (-e $tissue_directory){
        @tissue = read_dir($tissue_directory, prefix => 1);
     }
     
     my $manual_count = 0;
     foreach my $d (@descriptions) {
#       if (grep {$_ =~/$d/} @manual){
         chomp($d);
         $d =~ s{\.[^.]+$}{};
         $d =~ s{.*/}{}; 
      foreach my $m (@manual){
         chomp($m);
         $m =~ s{\.[^.]+$}{};
         $m =~ s{.*/}{}; 
       if ($m =~ $d){
        $manual_count++;
        }
       }
     }

     my $orthology_count = @orthology;
     my $process_count   = @process;
     my $function_count  = @function;
     my $component_count = @component;
     my $tissue_count = @tissue;
     my $gene_ec_count = @gene_ec;
     my $mole_ec_count = @mole_ec;
     my $anatomy_ec_count = @anatomy_ec;

     my $count = @descriptions;
     my $name = $species_name_hash{$species};
        $name    =~s/\_/ /g;
     my $name_list = "$name\n";
     my $commify_count = commify($count);
     my $commify_orthology_count = commify($orthology_count);
     my $commify_process_count = commify($process_count);
     my $commify_component_count = commify($component_count);
     my $commify_function_count = commify($function_count);
     my $commify_tissue_count = commify($tissue_count);
     my $commify_gene_ec_count = commify($gene_ec_count);
     my $commify_mole_ec_count = commify($mole_ec_count);
     my $commify_anatomy_ec_count = commify($anatomy_ec_count);

     write_file($summary, {append => 1 }, $name_list);
     my $text = "$commify_count individual gene descriptions\n";
     write_file($summary, {append => 1 }, $text);
     my $commify_manual_count = commify($manual_count);
     my $manual_text = "$commify_manual_count genes have manual descriptions\n";
     $count_manual_total+=$manual_count;
     write_file($summary, {append => 1 }, $manual_text);
     my $subtract_count = $count - $manual_count;
     my $commify_subtract_count = commify($subtract_count);
     $count_automated_only_total+=$subtract_count;
     my $subtract_count_text = "$commify_subtract_count genes have only automated descriptions\n";
     write_file($summary, {append => 1 }, $subtract_count_text);
     my $orthology_text = "$commify_orthology_count orthology sentences\n";
     write_file($summary, {append => 1 }, $orthology_text);
     my $process_text = "$commify_process_count gene ontology process sentences\n";
     write_file($summary, {append => 1 }, $process_text);
     my $function_text = "$commify_function_count gene ontology molecular function sentences\n";
     write_file($summary, {append => 1 }, $function_text);
     my $component_text = "$commify_component_count gene ontology cellular component sentences\n";
     write_file($summary, {append => 1 }, $component_text);
     my $tissue_text = "$commify_tissue_count tissue expression sentences\n";
     write_file($summary, {append => 1 }, $tissue_text);
     my $gene_ec_text = "$commify_gene_ec_count gene expression cluster sentences\n";
     write_file($summary, {append => 1 }, $gene_ec_text);
     my $mole_ec_text = "$commify_mole_ec_count molecule expression cluster sentences\n";
     write_file($summary, {append => 1 }, $mole_ec_text);
     my $anatomy_ec_text = "$commify_anatomy_ec_count anatomy expression cluster sentences\n";
     write_file($summary, {append => 1 }, $anatomy_ec_text);

     write_file($summary, {append => 1 }, $single_breakln);

     $kount += $count;
   }
}
     my $commify_kount = commify($kount);
     write_file($summary, {append => 1 }, $breakln);
     my $text = "Total number of concise descriptions is $commify_kount\n";
     write_file($summary, {append => 1 }, $text);
     my $commify_automated_total = commify($count_automated_only_total);
     my $text_automated = "Total number of concise descriptions that are automated only is $commify_automated_total\n";
     write_file($summary, {append => 1 }, $text_automated);
     my $commify_manual_total = commify($count_manual_total);
     my $text_manual = "Total number of concise descriptions that are manually curated as well is $commify_manual_total\n";
     write_file($summary, {append => 1 }, $text_manual);
#
exit 0;

sub commify {
        local $_  = shift;
        1 while s/^(-?\d+)(\d{3})/$1,$2/;
        return $_;
}
