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
my $elegans_script = "./create_sentence_multiple_orthologs_elegans_parallel.pl";
my $parasite_script = "./create_sentence_multiple_orthologs_parasites.pl";
my $species_script = "./create_sentence_multiple_orthologs_species_parallel.pl";
#
my $AND = "AND";
my $species_project_name_prefix = $ARGV[0];
chomp($species_project_name_prefix);
my ($species, $project, $name, $prefix) = split(/$AND/, $species_project_name_prefix);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $name =~ s/^\s+//;
 $name =~ s/\s+$//;
 $prefix  =~s/^\s+//;
 $prefix  =~s/\s+$//;
#
if ($species =~/elegans/){
  @args=($elegans_script, $species, $project);
  $status=system(@args);
 } elsif (($species =~/malayi/) or ($species =~/volvulus/) or ($species =~/ratti/)) {
  @args=($parasite_script, $species, $project, $name, $prefix);
  $status=system(@args); 
 } else {
  @args=($species_script, $species, $project, $name, $prefix);
  $status=system(@args);
 }

exit 0;

