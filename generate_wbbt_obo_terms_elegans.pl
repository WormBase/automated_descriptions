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
my $parse_script = "./01_get_obo_terms_only.pl";
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $tissue_dir = $home_elegans . "tissue_expression/input_files/";
#
 my $input =  $tissue_dir . "WBbt.obo";
 my $output = $tissue_dir . "anatomy_terms.txt";
 @args=($parse_script, $input, $output);
 $status=system(@args);
exit 0;
