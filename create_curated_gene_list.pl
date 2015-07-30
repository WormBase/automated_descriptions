#!/usr/bin/env perl
use strict;
use diagnostics;
use DBI;
use List::MoreUtils qw(uniq);
use ConciseDescriptions;
use File::Slurp;

my $ip_file = "./db_ip.txt";
my $db_ip = read_file($ip_file);
chomp($db_ip);
$db_ip =~ s/^\s+//;
$db_ip =~ s/\s+$//;

my $dbh = DBI->connect ( "dbi:Pg:dbname=testdb;host=$db_ip", "acedb", "") or die "Cannot connect to database!\n";

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
my $gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $outfile = $gene_list_dir . "sort.curated_genes.txt";
# if the output file exists delete it.
if (-e $outfile){
   my @args = ("/bin/rm", "-f", $outfile);
   system(@args) == 0 or die("could not delete file $outfile\n");
}
print "Getting data from postgres DB...\n";
my $sql_query = "SELECT con_wbgene.joinkey, con_wbgene, con_desctype,  con_desctext FROM con_wbgene, con_desctype, con_desctext WHERE (con_desctype.joinkey=con_wbgene.joinkey) AND (con_desctext.joinkey=con_wbgene.joinkey) ORDER by con_wbgene;";
my $result = $dbh->prepare($sql_query);
$result->execute() or die "Cannot prepare statement: $DBI::errstr\n";
my $term;
my %genes=();
my @curated_genes;
while (my @row = $result->fetchrow) {
        my $con_wbgene   = $row[1];
        my $con_desctype = $row[2];
        my $con_desctext = $row[3];
#        print "$con_wbgene\t$con_desctype\t$con_desctext$\n";
        if (($con_desctype =~ /Concise/) and ($con_desctext ne "")) {
         push(@curated_genes, $con_wbgene);
       } 
    }

my @unique_curated_genes = uniq(@curated_genes);
for my $gene (sort @unique_curated_genes){
    my $out = "$gene\n";
    write_file($outfile, {append => 1 }, $out);
}

print "Output stored in $outfile\n";
1;
