#!/usr/bin/env perl 
use ConciseDescriptions;
use Carp;
use LWP::Simple;
use OBO::Parser::OBOParser;
use File::Slurp;
use strict;
use warnings;
#
# Authors: J. Done and R. Kishore, California Institute of Technology, 2014. 
# http://wiki.wormbase.org/index.php/Generation_of_automated_descriptions
#
# The path for the output file is $functions_dir; 
# the path for the individual gene descriptions is $individual_path.
#
my $html = ConciseDescriptions::get_html_dir();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $path = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $url = 'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo';
my $text = get $url;

my $output = $path . "goid_altid.txt";
if (-e $output){
   my @args = ("/bin/rm", "-f", $output);
   system(@args) == 0 or die("could not delete file $output\n");
}

$text =~s/expand\_assertion\_to\:\s\"Class\:\s\<http\:\/\/www\.w3\.org\/2002\/07\/owl\#Nothing\>\sEquivalentTo\:\s\?X\sand\s\(RO\_0002162\ssome\s\?Y\)\"\s\[\]//gi;
$text =~s/property\_value\: propformat\-version \"1\.2\" xsd\:string//gi;

my @terms;
my $my_term="";
my @alts = ();
my $id = 0;
my $old_id=0; 
my $name="";
 my @lines = split(/\n/,$text);
   foreach (@lines){
     my $line = $_;
     chomp($line);
       if ($line =~/^name\:/){
           $name = $line;
           $name =~ s/name\://g;
           $name =~ s/^\s+//;
           $name =~ s/\s+$//;
       }
       if (($line =~/^id\:/)or ($line =~/^alt\_id\:/)){
           if ($line=~/^id\:/){
                   $id = $line;
                   $id =~s/id\://g;
                   $id =~ s/^\s+//;
                   $id =~ s/\s+$//;
               if ($old_id ne $id){
                   my $alt_string = join("\-", @alts); 
                   if ($alt_string ne "" ) {
                         my $output_string = "$old_id\t$alt_string\t$name\n";
                         write_file($output, {append => 1 }, $output_string); 
                    }
                   $name ="";
                   $old_id = $id;
                   @alts=();
                   $alt_string="";
                  }
               }
          if ($line =~/^alt\_id\:/){
              my $alt = $line;
              $alt =~s/alt\_id\://g;
              $alt=~ s/^\s+//;
              $alt =~ s/\s+$//;
              push(@alts, $alt);
              }
     }
 }

exit 0;
