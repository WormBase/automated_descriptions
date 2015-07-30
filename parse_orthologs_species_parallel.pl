#!/usr/bin/env perl
use warnings;
use strict;
#
my @args;
my $status;
#
my $sapiens_script = "./parse_orthologs_create_sapiens_file.pl";
my $orthologs_script="./parse_elegans_orthologs_species.pl";
my $malayi_script = "./parse_malayi_orthologs_species.pl";
my $volvulus_script = "./parse_volvulus_orthologs_species.pl";
#
my $AND = "AND";
my $species_project_name = $ARGV[0];
chomp($species_project_name);
my ($species, $project, $name) = split(/$AND/, $species_project_name);
 $species =~s/^\s+//;
 $species =~s/\s+$//;
 $project =~s/^\s+//;
 $project =~s/\s+$//;
 $name =~ s/^\s+//;
 $name =~ s/\s+$//;

if ($species =~/elegans/){
  @args=($sapiens_script, $species, $name);
  $status=system(@args);
 } elsif ($species =~/malayi/){
  @args=($volvulus_script, $species);
  $status=system(@args); 
  @args=($orthologs_script, $species);
  $status=system(@args);
 } elsif ($species =~/volvulus/){
  @args=($malayi_script, $species);
  $status=system(@args);
  @args=($orthologs_script, $species);
  $status=system(@args);
 } elsif ($species =~/ratti/){
  @args=($volvulus_script, $species);
  $status=system(@args); 
  @args=($malayi_script, $species);
  $status=system(@args);
  @args=($orthologs_script, $species);
  $status=system(@args);
 } else {
  @args=($orthologs_script, $species);
  $status=system(@args);
 }

exit 0;
