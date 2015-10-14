#!/usr/bin/env perl
#
# Author: J. Done, California Institute of Technology
# This perl script is written to address the needs of 
# writing concise descriptions for tissue expression.
#
use strict;
use diagnostics;
use DBI;
use File::Slurp;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use Data::Dumper;
use ConciseDescriptions;
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $tissue_dir = $home_elegans . "tissue_expression/";
my $gene_list_dir = $home_elegans . "gene_lists/" ;
my $input_path = $tissue_dir . "input_files/";
my $output_path = $tissue_dir . "output_files/";
my $gene_dir = $output_path . "individual_gene_sentences/";
# Define the ace file and ontology for anatomy
my $input_ace_file = $input_path . "anatomy_obo_terms.ace";
my $output_ace_file = $input_path . "anatomy_obo_terms_ace_names.txt";
my $ontology = "WBbt";
#
my %anatomy_ontology;
my %reverse_anatomy_ontology;
# $anatomy_file is created by the script, get_obo_terms_only.pl, for anatomy ontology 
my $anatomy_file = $input_path . "anatomy_terms.txt";
my @anat_lines = read_file($anatomy_file); 
my @keys_anatomy;
foreach my $anat_line (@anat_lines){
 chomp($anat_line);
 my $key   = $anat_line;
 my $value = $anat_line;
 $key  =~ /\((WBbt.+?)\)/;
 $key = $1; 
 $value =~ s/\(.*//g;
 chomp($value);
 $value =~ s/^\s+//;
 $value =~ s/\s+$//;
 chomp($key);
 $key =~ s/^\s+//;
 $key =~ s/\s+$//;
 push(@keys_anatomy, $key);
 $anatomy_ontology{$key} = $value;
 $reverse_anatomy_ontology{$value} = $key;
}
#
my @out_lines=();
my @ace_lines = read_file($input_ace_file);
foreach my $line (@ace_lines){
 chomp($line);
 my $key;
 my $out_line;
 if ($line =~/\"(WBbt.+?)\"/){
      $key = $1;
      $key =~ s/^\s+//;
      $key =~ s/\s+$//;
  $out_line = $line . "\t" . $anatomy_ontology{$key} . "\n"; 
 } else {
  $out_line = $line . "\n";
 }
  push(@out_lines, $out_line);
}
#
    write_file($output_ace_file, @out_lines);
#
exit 0;


