#!/usr/bin/env perl
use strict;
use warnings;
use File::Slurp;
use ConciseDescriptions;
#
# J.Done and R. Kishore, California Institute of Technology, 2015.
#
# Requires GNU parallel installed
#
my $directory_creation_script = "./directory_creation.pl";
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
my $production_release = $ARGV[0];
print "creating directories for $production_release\n";
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $species_file = $home . "species.txt";
my @lines = read_file($species_file);
#
my $species_list="";
my @species_array=();
#
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

 my $term = $species . "AND" . $production_release;
 $species_list .= $term . " ";
 push(@species_array, $term);
 print "$species_list\n";
}

 $species_list  =~s/^\s+//;
 $species_list  =~s/\s+$//;
 my @args=($parallel_exec, $j_flag, $percent, $directory_creation_script, $three_colons, @species_array);
 my $status=system(@args);

exit 0;
