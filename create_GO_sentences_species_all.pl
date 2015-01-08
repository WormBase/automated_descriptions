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
my $parse_script = "./create_GO_sentences_species.pl";
my $parse_elegans_script = "./create_GO_sentences_elegans.pl";
#
# $species_file holds a list of species and project numbers studied by Wormbase.
my $species_file = $home . "species.txt";
my @lines = read_file($species_file);
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
 if ($prefix ne ""){
 print "$species\t$project\n";
 @args=($parse_script, $species, $project, $name, $prefix);
 $status=system(@args);
 } else {
 print "$species\t$project\n";
 @args=($parse_elegans_script, $species, $project, $name, $prefix);
 $status=system(@args);
 }
}
#
exit 0;
