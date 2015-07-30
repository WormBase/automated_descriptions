#!/usr/bin/env perl
use warnings;
use strict;
use File::Slurp;
use List::MoreUtils qw(uniq);
use diagnostics;
use DBI;
use POSIX qw/strftime/;
use ConciseDescriptions;
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $outpath = $home . "gene_lists/";
 
my $output_file = $outpath . "wb_gene_list.txt";
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
my $ambiguous_file = $outpath . "ambiguous.txt";
# if the output file exists delete it.
if (-e $ambiguous_file){
   my @args = ("/bin/rm", "-f", $ambiguous_file);
   system(@args) == 0 or die("could not delete file $ambiguous_file\n");
}
my $ambiguity_gsa_file = $outpath . "ambiguity_gsa_file.txt";
# if the output file exists delete it.
if (-e $ambiguity_gsa_file){
   my @args = ("/bin/rm", "-f", $ambiguity_gsa_file);
   system(@args) == 0 or die("could not delete file $ambiguity_gsa_file\n");
}
# 
# Define the IP of the database used
# The IP address will not be hard coded into the script 
# for reasons of security and maintainability.  
# Ask Juancarlos or the PostgreSQL DB manager for IP address if necessary.
my $ip_file = "./db_ip.txt";
my $db_ip = read_file($ip_file);
chomp($db_ip);
$db_ip =~ s/^\s+//;
$db_ip =~ s/\s+$//;
#
my $wbgeneid_hash_ref = get_wbgeneid_hash($db_ip, $output_file, $ambiguous_file, $ambiguity_gsa_file);
my %wbgeneid_hash = %$wbgeneid_hash_ref;

#my $file = "./ambiguity";
#my @names = read_file($file);
#foreach my $name (@names){
# chomp($name);
# $name =~ s/^\s+//;
# $name =~ s/\s+$//;
#}
exit 0;

sub get_wbgeneid_hash{
 
 my $db_ip = shift;
 my $output_file=shift;
 my $ambiguous_file=shift;
 my $ambiguity_gsa_file=shift;
 my %wbgeneid_hash=();
 my @ambiguity_array=();
#
# This subroutine returns a hash of gene ids given the gene name
#
 my $dbh = DBI->connect ( "dbi:Pg:dbname=testdb;host=$db_ip", "acedb", "") or die "Cannot connect to database!\n"; 

 print STDERR "Processing gin_synonyms...\n";
 my %gin_synonyms = ();
 my %gin_syntype = ();
 my $result = $dbh->prepare("SELECT * FROM gin_synonyms WHERE gin_synonyms IS NOT NULL");
 $result->execute() or die "Cannot prepare statement: $DBI::errstr\n"; 
 while (my @row = $result->fetchrow) {
    if ($row[1]) {
	$gin_synonyms{$row[0]}{$row[1]} = 1;
#        print "row 1 is $row[1]\n";
    }
    if ($row[2]) {
	$gin_syntype{$row[0]}{$row[2]} = 1;
#        print "row 2 is $row[2]\n";
    }
 }

 print STDERR "Processing gin_locus...\n";
 my %gin_locus = ();
 $result = $dbh->prepare("SELECT * FROM gin_locus WHERE gin_locus IS NOT NULL");
 $result->execute() or die "Cannot prepare statement: $DBI::errstr\n"; 
 while (my @row = $result->fetchrow) {
    if ($row[1]) {
	$gin_locus{$row[0]}{$row[1]} = 1;
    }
 }

 print STDERR "Processing gin_wbgene...\n";
 my %gin_wbgene = ();
 $result = $dbh->prepare("SELECT * FROM gin_wbgene WHERE gin_wbgene IS NOT NULL");
 $result->execute() or die "Cannot prepare statement: $DBI::errstr\n"; 
 while (my @row = $result->fetchrow) {
    if ($row[1]) {
	$gin_wbgene{$row[0]}{$row[1]} = 1;
#        print "$row[0]\t$row[1]\t$row[2]\t$row[3]\n";
    }
 }
 foreach my $jk (keys % gin_synonyms) {
    my @locus = keys % {$gin_locus{$jk}};
    if (scalar (@locus) == 1) {
        my @aux = (keys % {$gin_wbgene{$jk}}, $locus[0], keys % {$gin_synonyms{$jk}});
        my $out = join (",", @aux);
        my $out_string = $aux[0] . "\t" . $aux[1] . "\n";
#	if (is_known_gene($locus[0])) {
	if ($locus[0]) {
        foreach my $name (@aux){
#         print "locus 0 is $locus[0]\n";
         if ($wbgeneid_hash{$name}){
          my @name_array = split(/,/,$wbgeneid_hash{$name});
          my @unique_name_array = uniq(@name_array);
          foreach my $element (@unique_name_array){
           if ($element ne $aux[0]){
              my $ambiguity="$element\t$name\t$aux[0]\n";
              push(@ambiguity_array, $name);
             if ($ambiguous_file){
               write_file($ambiguous_file, {append => 1 }, $ambiguity);
              }
           if ($name ne $aux[0]){
               $wbgeneid_hash{$name} .="\,".$aux[0];
           }
          }
         }
        } else{
           if ($name ne $aux[0]){
            $wbgeneid_hash{$name} = $aux[0];
          }
        }
        }
        if ($output_file){
          write_file($output_file, {append => 1 }, $out_string);
#          write_file($output_file, {append => 1 }, "\n");
     }
  }
 }
}

$dbh->disconnect;

my @unique_ambiguity_array = uniq(sort(@ambiguity_array));
foreach my $element (@unique_ambiguity_array){
        if ($ambiguity_gsa_file){
          write_file($ambiguity_gsa_file, {append => 1 }, $element);
          write_file($ambiguity_gsa_file, {append => 1 }, "\n");
        }
}

return \%wbgeneid_hash;

}

sub is_known_gene {
    my $s = shift;
    if ($s =~ /^(WBGene|GENEPREDICTION|Cr|Cbr|Cbg|Cbn|Cjp|Hpa|Oti|Ppa|Cja|Ovo|Bma|Cre)/i) {
#    if ($s =~ /^(WBGene|GENEPREDICTION)/i) {
#        print "$s\n";
	return 0;
    }
    return 1;
}
