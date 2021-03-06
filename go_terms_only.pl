#!/usr/bin/env perl 
use Carp;
use LWP::Simple;
use OBO::Parser::OBOParser;
use File::Slurp;
use ConciseDescriptions;
use strict;
use warnings;
#
my $html = ConciseDescriptions::get_html_dir();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $path = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";

my $url= 'http://geneontology.org/ontology/go.obo';
my $input_file = $path . 'go.obo';
my $output_file = $path . "go_terms.txt";

getstore($url, $input_file);

my $file_size= -s $input_file;

if ($file_size gt 0) {

my $text = read_file( $input_file ) ;
$text =~s/expand\_assertion\_to\:\s\"Class\:\s\<http\:\/\/www\.w3\.org\/2002\/07\/owl\#Nothing\>\sEquivalentTo\:\s\?X\sand\s\(RO\_0002162\ssome\s\?Y\)\"\s\[\]//gi;
$text =~s/property\_value\: propformat\-version \"1\.2\" xsd\:string//gi;
$text =~s/property\_value\: IAO\:0000589 \"cell and encapsulating structures\" xsd\:string//gi;
$text =~s/expand\_assertion\_to\: \"Class\: \?X DisjointWith\: RO\_0002162 some \?Y\" \[\]//gi;
$text =~s/EXACT \[\S+\]/EXACT \[\]/g;
$text =~s/NARROW \[\d+\, PMID\:\d+\]/NARROW \[\]/g;
$text =~s/NARROW \[PMID\:\d+\]/NARROW \[\]/g;
$text =~s/RELATED \[\d+\, PMID\:\d+\]/RELATED \[\]/g;

write_file( $input_file, \$text ) ;
			
my $my_parser = OBO::Parser::OBOParser->new();
my $ontology  = $my_parser->work($input_file);

my @terms;
my $my_term="";

foreach my $term (sort {$a->id() cmp $b->id()} @{$ontology->get_terms()}) {
        if ( defined ($term->id()) && defined($term->name()) ) {
               my $id = $term->id();
               my $name = $term->name();
            if ( $term->is_obsolete() ){
                  print "$name\t$id is obsolete\.\n";
#                 $my_term = $term->name() . "\t \(" . $term->id() . "\)" . " \( Obsolete \) \n";
               }
               else{
                    if (lc ($name) !~/obsolete/){
                      $my_term = $name . "\t \(" . $id . "\)\n";
                    }
                   }
	    push @terms, $my_term;
        }
}

my @output_terms = sort(@terms);

write_file( $output_file, @output_terms );

}

exit 0;
