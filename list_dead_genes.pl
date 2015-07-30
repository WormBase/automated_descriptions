#!/usr/bin/env perl
use warnings;
use strict;
use File::Slurp;
#use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use diagnostics;
use DBI;
use ConciseDescriptions;

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
#my $gene_list_dir = "./";
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";

my $ip_file = "./db_ip.txt";
my $db_ip = read_file($ip_file);
chomp($db_ip);
$db_ip =~ s/^\s+//;
$db_ip =~ s/\s+$//;

my $output_file =  $gene_list_dir . "sort.dead_genes.txt";
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}

my $dbh = DBI->connect ( "dbi:Pg:dbname=testdb;host=$db_ip", "acedb", "") or die "Cannot connect to database!\n"; 

 print STDERR "Processing gin_dead...\n";
 my %gin_dead = ();
# my $sql_query = "SELECT gin_wbgene.joinkey, gin_dead, gin_wbgene FROM gin_dead, gin_wbgene WHERE (gin_dead.joinkey=gin_wbgene.joinkey) ORDER by gin_wbgene;";
my $sql_query = "SELECT gin_wbgene.joinkey, gin_dead, gin_wbgene FROM gin_dead, gin_wbgene WHERE (gin_dead.joinkey=gin_wbgene.joinkey) ORDER by gin_wbgene;";
 my $result = $dbh->prepare($sql_query);
 $result->execute() or die "Cannot prepare statement: $DBI::errstr\n"; 
 while (my @row = $result->fetchrow) {
    if (($row[1]) and ($row[2])) {
	print "$row[2] is $row[1]\n";
        my $out = $row[2] . "\n";
        write_file($output_file, {append => 1 }, $out);
    }
   
 }

exit 0;
