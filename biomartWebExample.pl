#!/usr/bin/env perl
# an example script demonstrating the use of BioMart webservice
use strict;
use LWP::UserAgent;
use File::Slurp;
use ConciseDescriptions;

my $html = ConciseDescriptions::get_html_dir();
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
#
my $outpath = $html . "concise_descriptions/semantic_categories/homology/";
my $species = "c_elegans";
my $outpath  = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/$species/orthology/input_files/";

my $output_file = $outpath . "HumanIDs_mart_export.txt";
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}

open (FH,$ARGV[0]) || die ("\nUsage: biomartWebExample.pl Query.xml\n\n");

my $xml;
while (<FH>){
    $xml .= $_;
}
close(FH);


my $path="http://www.biomart.org/biomart/martservice?";
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
my $ua = LWP::UserAgent->new;

my $response;

$ua->request($request, 
	     sub{   
		 my($data, $response) = @_;
		 if ($response->is_success) {
                  if ($data !~/ERROR/){
                        write_file($output_file, {append => 1 }, $data);
                     }
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);

exit 0;
