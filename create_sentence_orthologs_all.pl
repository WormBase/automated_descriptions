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
my $parse_script = "./create_sentence_orthologs.pl";
my $parse_elegans_script = "./create_sentence_orthologs_elegans.pl";
#
# $species_file holds a list of species and project numbers studied by Wormbase.
my $species_file = $home . "species.txt";
my @lines = read_file($species_file);
#
foreach my $line (@lines){
 my ($species, $project) = split(/\t/,$line);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;

 if ($species !~/c_elegans/){
  @args=($parse_script, $species, $project);
  $status=system(@args);
 } elsif ($species =~/c_elegans/) {
  @args=($parse_elegans_script, $species, $project);
  $status=system(@args);
 }
}
#
exit 0;
