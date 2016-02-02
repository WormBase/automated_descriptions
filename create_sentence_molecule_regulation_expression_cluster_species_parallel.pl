#!/usr/bin/env perl
#
# Author: J. Done, California Institute of Technology
# This perl script is written to address the needs of 
# writing concise descriptions for tissue expression.
#
use strict;
use diagnostics;
use DBI;
use File::Slurp;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use Data::Dumper;
use ConciseDescriptions;
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
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
 if ($species=~/elegans/){
   $prefix = "ce";
 } else {
   $prefix = lc $prefix;
 }
#
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/$species/";
my $mr_ec_dir = $home . "molecule_regulation_expression_cluster/";
my $elegans_gene_list_dir = $home_elegans . "gene_lists/" ;
my $gene_list_dir = $home . "gene_lists/" ;
my $input_path = $mr_ec_dir . "input_files/";
my $output_path = $mr_ec_dir . "output_files/";
my $gene_dir = $output_path . "individual_gene_sentences/";
my $input_file = $input_path . $prefix . "ECsummary_molReg.$RELEASE.txt";

if (-e $input_file) {
my $input_file_size = -s $input_file;
if ($input_file_size > 1){
my @ec_lines = read_file($input_file); 
my %gene_relate=();
my %gene_exp=();
my %molecule_regulator=();

my @gene_list;
   foreach (@ec_lines){
     my $ec_line = $_;
     chomp($ec_line);
     next if ($ec_line =~/GeneID/);
      my ($no, $name, $relate, $regulator, $exp) = split(/\t/,$ec_line);
      if ($no){
          $gene_relate{$no} = $relate;
          $gene_exp{$no} = $exp;
          $molecule_regulator{$no} = $regulator;
      } 
      if ($no){
       push(@gene_list, $no);
      }
    }
#
# if the output files for individual gene sentences exist delete them.
   my @individual_files = glob("$gene_dir/WBGene*");
      foreach my $individual (@individual_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
my $output_log = $output_path . "sentences_for_molecule_regulation_expression_cluster.txt";
# if the output file exists delete it.
if (-e $output_log){
   my @args = ("/bin/rm", "-f", $output_log);
   system(@args) == 0 or die("could not delete file $output_log\n");
}
my $output_test_log = $output_path . "test_sentences_for_molecule_regulation_expression_cluster.txt";
# if the output file exists delete it.
if (-e $output_test_log){
   my @args = ("/bin/rm", "-f", $output_test_log);
   system(@args) == 0 or die("could not delete file $output_test_log\n");
}
my $output = $output_path . "molecule_regulation_expression_cluster.txt";
my $output_test = $output_path . "test_molecule_regulation_expression_cluster.txt";
my $output_genes = $output_path . "gene_ids_molecule_regulation_expression_cluster.txt";
my $output_test_genes = $output_path . "test_gene_ids_molecule_regulation_expression_cluster.txt";
# if the output file exists delete it.
if (-e $output){
   my @args = ("/bin/rm", "-f", $output);
   system(@args) == 0 or die("could not delete file $output\n");
}
if (-e $output_test){
   my @args = ("/bin/rm", "-f", $output_test);
   system(@args) == 0 or die("could not delete file $output_test\n");
}
if (-e $output_genes){
   my @args = ("/bin/rm", "-f", $output_genes);
   system(@args) == 0 or die("could not delete file $output_genes\n");
}
if (-e $output_test_genes){
   my @args = ("/bin/rm", "-f", $output_test_genes);
   system(@args) == 0 or die("could not delete file $output_test_genes\n");
}
# infile is the sorted file containing the list of gene ids with no concise description
my $wbgene_id_name_file = $gene_list_dir . "geneIDs.txt";
my @files = read_file($wbgene_id_name_file); 
my %gene_name=();
   foreach (@files){
     my $file_line = $_;
     chomp($file_line);
      my ($no, $file, $name, $alt, $live) = split(/\,/,$file_line);
          $name =~ s/^\s+//;
          $name =~ s/\s+$//;
          $alt =~ s/^\s+//;
          $alt =~ s/\s+$//;
          $file =~ s/^\s+//;
          $file =~ s/\s+$//;
      if ($name){
          $gene_name{$file} = $name;
      } elsif ($alt) {
          $gene_name{$file} = $alt;
      }
    }
my $dead_gene_list = $elegans_gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);
my @live_gene_array = ();
foreach my $test (@gene_list){
 my $keep =0;
 foreach my $dead (@dead_genes){
   if ($dead =~/$test/){
    $keep = 1;
   }
 }
 if ($keep ==0){
   push(@live_gene_array, $test);
 }
}
my @sorted_live_gene_array = sort(@live_gene_array);
foreach my $gene_id (@sorted_live_gene_array) {
    chomp($gene_id);
          $gene_id =~ s/^\s+//;
          $gene_id =~ s/\s+$//;
    my $name = $gene_name{$gene_id};
    if ($name) {
     sleep (0.001)
    } else {
        print "gene is $gene_id has no name\n";
    }
    if ($name){
    my $gene_file = $gene_dir . $gene_id;
    my $sentence = "";
    my $relation = lc $gene_relate{$gene_id};
    my $experiment = $gene_exp{$gene_id};
    my @exps = split(/\,/, $experiment);
    next if (@exps <= 0);
    my $experiments = "";
    if (@exps==1){
        $experiments = $exps[0];
    } elsif (@exps == 2){
        $experiments =  "$exps[0] and $exps[1]";
    } elsif (@exps > 2) {
        my $c = 0;
        foreach my $e (@exps){
         $c++;
         if ($c < @exps){
          $experiments .= "$e\, ";
         } else {
          $experiments .= "and $e";
         }
        }
    }   
        $experiments =~s/RNA seq/RNA sequencing/g;
        $experiments =~s/RNA\-seq/RNA sequencing/g;
        $experiments =~s/study//gi;
        $experiments = lc $experiments;
        $experiments =~s/rna sequencing/RNA sequencing/g;
        $experiments =~s/ \,/\,/g;
        $experiments =~s/pcr/PCR/g;

    my $regulator = $molecule_regulator{$gene_id};
    my $molecule = "";
    my @molecules = split(/\,/,$regulator);
        if (@molecules==1){
        $molecule = ucfirst $molecules[0];
    } elsif (@molecules == 2){
        my $m0 = ucfirst $molecules[0];
        my $m1 = ucfirst $molecules[1];
        $molecule = "$m0 and $m1";
    } else {
        my $c = 0;
        foreach my $m (@molecules){
         $c++;
         $m = ucfirst $m;
         if ($c < @molecules){
          $molecule .= "$m\, ";
         } else {
          $molecule .= "and $m";
         }
        }
    }
    $molecule =~s/Adsorbable/adsorbable/g;
    $molecule =~s/nanotube/nanotubes/g;
    $molecule =~s/nanotubess/nanotubes/g;
    $molecule =~s/Substances/substances/g;
    $molecule =~s/Single\-walled/single\-walled/g;
    if (@exps > 0 ) {
     $sentence = "$experiments studies indicate that $name is regulated by $molecule"; 
    }
    $sentence =~s/bodywall/body wall/gi;
    $sentence =~s/coelomocyte/coelomocytes/gi;
    $sentence .="\;\n";
    $sentence =~s/'the all'/'all'/g;
    $sentence =~s/the certain/certain/g;
    $sentence =~s/compound/compounds/gi;
    $sentence =~s/  / /g;

#    $sentence = ucfirst($sentence);
     write_file($output_test_log, {append => 1 }, $sentence);
     write_file($output_test_log, {append => 1 }, "\n\n\n");
    if (@molecules < 21) {
     write_file($output_log, {append => 1 }, $sentence);
     write_file($output_log, {append => 1 }, "\n\n\n");
     write_file($gene_file, $sentence);
    }
   }
  }
 }
}
exit 0;
