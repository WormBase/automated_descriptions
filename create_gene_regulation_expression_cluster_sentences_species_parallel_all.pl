#!/usr/bin/env perl
use warnings;
use strict;
use ConciseDescriptions;
use File::Slurp;
#
my @args;
my $status;
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $parallel_file = "./parallel_path.txt";
my $parallel_path = read_file($parallel_file);
 chomp($parallel_path);
 $parallel_path =~s/^\s+//;
 $parallel_path =~s/\s+$//;
my $parallel_exec = $parallel_path . "parallel";
my $three_colons = "\:\:\:";
my $j_flag = "\-j";
my $percent = "85\%";
#
my $script = "./create_sentence_gene_regulation_expression_cluster_species_parallel.pl";
#
# $species_file holds a list of species and project numbers studied by Wormbase.
my $species_file = $home . "species.txt";
my @lines = read_file($species_file);
#
my $elegans="";
my @species_array=();
foreach my $line (@lines){
 my ($species, $project, $name, $prefix) = split(/\t/,$line);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $name    =~s/^\s+//;
 $name    =~s/\s+$//;
 $name    =~s/ /\_/g;
 $prefix  =~s/^\s+//;
 $prefix  =~s/\s+$//;
 my $term = $species . "AND" . $project . "AND" .  $name . "AND" . $prefix;
 if ($species){
  push(@species_array, $term);
 }
}

 @args=($parallel_exec, $j_flag, $percent, $script, $three_colons, @species_array);
 $status=system(@args);

#
exit 0;
