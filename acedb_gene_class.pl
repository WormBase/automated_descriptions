#!/usr/bin/env perl
use strict;
use warnings;
use ConciseDescriptions;
use Ace;        # Import the AcePerl library
use File::Slurp;

my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $elegans_gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
    
my $output_file = $elegans_gene_list_dir . "acedb_gene_class.txt";
if (-e $output_file) {
 my @args = ("/bin/rm", "-f", $output_file);
 my $status = system(@args);
 }
my $host = read_file("./acedb_host.txt");
chomp($host);
my $port = read_file("./acedb_port.txt");
chomp($port);
my $db = Ace->connect( -host => $host, -port => $port )
        or die "Can't connect to '$host' on port '$port' : ", Ace->error;

#if ($db) {
#    print "Connection great !\n";
# }

my $aql_query = "select item, item->Genes from item in class Gene_Class";
my $count   = $db->aql($aql_query);
my @objects = $db->aql($aql_query);

#print "AceDB obtained $count objects\n"; 

foreach my $obj (@objects){
# In this example, AceDB returned an array of array references
 my @out_obj = @{$obj};
 while (my ($el1, $el2) = splice(@out_obj, 0, 2)) {
# In this example, 
        $el1 =~s/^\s+//;
        $el1 =~s/\s+$//;
        $el2 =~s/^\s+//;
        $el2 =~s/\s+$//;
        my $item = "$el1\t$el2\n";
#   print "$item";
        write_file($output_file, {append => 1}, $item);
 }
}
$db->close();

exit 0;
