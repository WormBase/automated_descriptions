#!/usr/bin/env perl
use warnings;
use strict;
use ConciseDescriptions;
use File::Slurp;
#
# J. Done and R. Kishore, California Institute of Technology, 2015.
# Assumes that GNU Parallel is installed
#
my @args;
my $status;
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $script = "./download_expression_cluster_summary_parallel.pl";
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
# $species_file holds a list of species and project numbers studied by Wormbase.
my $species_file = $home . "species.txt";
my @lines = read_file($species_file);
my @species_array = ();
foreach my $line (@lines){
 my ($species, $project, $name, $prefix) = split(/\t/,$line);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $name    =~s/^\s+//;
 $name    =~s/\s+$//;
 $prefix  =~s/^\s+//;
 $prefix  =~s/\s+$//;
 print "$species $prefix\n";
 my $term = $species . "AND" . $prefix;
 push(@species_array, $term);
}
 @args=($parallel_exec, $j_flag, $percent, $script, $three_colons, @species_array);
 $status=system(@args);
#
exit 0;