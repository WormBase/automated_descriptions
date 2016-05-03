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
#$PRODUCTION_RELEASE = "WS253"; ##liyuling
my $RELEASE = ConciseDescriptions::get_release();
#$RELEASE = "WS252"; ##liyuling
my $html = ConciseDescriptions::get_html_dir();
my $home_elegans = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $tissue_dir = $home_elegans . "tissue_expression/";
my $gene_list_dir = $home_elegans . "gene_lists/" ;
my $input_path = $tissue_dir . "input_files/";
my $output_path = $tissue_dir . "output_files/";
my $gene_dir = $output_path . "individual_gene_sentences/";
#
#
# if the output files for individual gene sentences exist delete them.
   my @individual_files = glob("$gene_dir/WBGene*");
      foreach my $individual (@individual_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
my $output_log = $output_path . "sentences_for_tissue_expression.txt";
# if the output file exists delete it.
if (-e $output_log){
   my @args = ("/bin/rm", "-f", $output_log);
   system(@args) == 0 or die("could not delete file $output_log\n");
}
my $output_test_log = $output_path . "test_sentences_for_tissue_expression.txt";
# if the output file exists delete it.
if (-e $output_test_log){
   my @args = ("/bin/rm", "-f", $output_test_log);
   system(@args) == 0 or die("could not delete file $output_test_log\n");
}
my $output = $output_path . "gene_tissue_expression.txt";
my $output_test = $output_path . "test_gene_tissue_expression.txt";
my $output_genes = $output_path . "gene_ids_tissue_expression.txt";
my $output_test_genes = $output_path . "test_gene_ids_tissue_expression.txt";
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
#
# Define the ace file and ontology for anatomy
my $input_ace_file = $input_path . "anatomy_obo_terms.ace";
my $ontology = "WBbt";
# 
my $gene_association_species_file = $input_path . "anatomy_association.$RELEASE.wb";
my $gene_wbbt_hash_ref = get_gene_wbbt_hash($gene_association_species_file);
my %genes = %$gene_wbbt_hash_ref;

# create arrays for neurons
my $vnc = $input_path . "ventral_cord_neuron.txt";
my $interneuron = $input_path . "interneuron.txt";
my $sensory = $input_path . "sensory_neuron.txt";
my $head_motor =  $input_path . "head_neuron.txt";
my $motor = $input_path . "motor_neuron.txt";
my $mechanosensory = $input_path . "mechanosensory_neuron.txt";
my $inner_labial = $input_path . "inner_labial_neuron.txt";
my $outer_labial = $input_path . "outer_labial_neuron.txt";
my $pn = $input_path . "pharyngeal_neuron.txt";
my $pmn = $input_path . "pharyngeal_motor_neuron.txt";
my $amphid = $input_path . "amphid.txt";
my $pharyngeal_interneuron = $input_path . "pharyngeal_interneuron.txt";

my $head_motor_neurons_ref=get_neurons($head_motor);
my @head_motor_neurons=@{$head_motor_neurons_ref};

my $vnc_neurons_ref=get_neurons($vnc);
my @vnc_neurons=@{$vnc_neurons_ref};

my $sensory_neurons_ref=get_neurons($sensory);
my @sensory_neurons=@{$sensory_neurons_ref};

my $interneurons_ref=get_neurons($interneuron);
my @interneurons=@{$interneurons_ref};

my $motor_neurons_ref=get_neurons($motor);
my @motor_neurons=@{$motor_neurons_ref};

my $mechanosensory_neurons_ref=get_neurons($mechanosensory);
my @mechanosensory_neurons=@{$mechanosensory_neurons_ref};

my $pharyngeal_neurons_ref=get_neurons($pn);
my @pharyngeal_neurons=@{$pharyngeal_neurons_ref};

my $pharyngeal_motor_neurons_ref=get_neurons($pmn);
my @pharyngeal_motor_neurons=@{$pharyngeal_motor_neurons_ref};

my $pharyngeal_interneurons_ref=get_neurons($pharyngeal_interneuron);
my @pharyngeal_interneurons=@{$pharyngeal_interneurons_ref};

my $outer_labial_neurons_ref=get_neurons($outer_labial);
my @outer_labial_neurons=@{$outer_labial_neurons_ref};

my $inner_labial_neurons_ref=get_neurons($inner_labial);
my @inner_labial_neurons=@{$inner_labial_neurons_ref};

my $amphid_neurons_ref=get_neurons($amphid);
my @amphid_neurons=@{$amphid_neurons_ref};
# build parent/child hash for granularity
my ($parents_ref, $children_ref) = ConciseDescriptions::get_ontology_parents_children($input_ace_file, $ontology);
my %parents = %$parents_ref;
my %children = %$children_ref;
my $instance_ref = ConciseDescriptions::get_ontology_instance($input_ace_file, $ontology);
my %instance = %$instance_ref;
# helper verb
my $helper_verb = " is expressed in the ";
my $helper_verb_widely = " is expressed widely";
my $helper_verb_several = " is expressed in several tissues including the ";
my $helper_verb_neuron  = " nervous system";
# infile is the sorted file containing the list of gene ids with no concise description
my $wbgene_id_name_file = $gene_list_dir . "wormbase_gene_id_name.list";
my @files = read_file($wbgene_id_name_file); 
my %gene_name=();
my @gene_list;
   foreach (@files){
     my $file_line = $_;
     chomp($file_line);
      my ($file, $name) = split(/\t/,$file_line);
      $gene_name{$file} = $name;
      push(@gene_list, $file);
    }
#my $uncurated_gene_file = $gene_list_dir . "sort.uncurated_genes.txt";
my $dead_gene_list = $gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);
#my @gene_array = read_file($uncurated_gene_file);

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
my $curated_gene_list = $gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);
#my @uncurated_genes_array = ();
#foreach my $test (@live_gene_array){
# my $keep =0;
# foreach my $curated (@curated_genes){
#   if ($curated =~/$test/){
#    $keep = 1;
#   }
# }
# if ($keep ==0){
#   push(@uncurated_genes_array, $test);
# }
#}
#my @sorted_uncurated_genes_array = sort(@uncurated_genes_array); 
# Create WBbt hash so that terms can be referenced by WBbt ID.
my %anatomy_ontology;
my %reverse_anatomy_ontology;
# $anatomy_file is created by the script, get_obo_terms_only.pl, for anatomy ontology 
my $anatomy_file = $input_path . "anatomy_terms.txt";
my @anat_lines = read_file($anatomy_file); 
my @keys_anatomy;
foreach my $anat_line (@anat_lines){
 chomp($anat_line);
 my $key   = $anat_line;
 my $value = $anat_line;
 $key  =~ /\((WBbt.+?)\)/;
 $key = $1; 
 $value =~ s/\(.*//g;
 chomp($value);
 $value =~ s/^\s+//;
 $value =~ s/\s+$//;
 chomp($key);
 $key =~ s/^\s+//;
 $key =~ s/\s+$//;
 push(@keys_anatomy, $key);
 $anatomy_ontology{$key} = $value;
 $reverse_anatomy_ontology{$value} = $key;
}
#
foreach my $gene_id (@sorted_live_gene_array) {
    chomp($gene_id);
          $gene_id =~ s/^\s+//;
          $gene_id =~ s/\s+$//;
    my $gene_file = $gene_dir . $gene_id;
    my $anatomy = $genes{$gene_id};
    my $name = $gene_name{$gene_id};
   # print "$name\t$anatomy\n";  liyuling
    if (($anatomy) and ($name)){

    my @wbbt_elements = split(/\,/,$anatomy);
    my @anatomies=();
    my @child_array=();
    my @parent_array=();

        foreach my $element (@wbbt_elements){
          chomp($element);
          $element =~s/WBbt\:0003679/WBbt\:0005735/g;
#          $element =~s/WBbt\:0003679/WBbt\:0005735/g;
         }
#
# Check if all, certain or none of the neuron list here.
# Substitute names of neurons
#
# 
my $skip;
my $substitute_neuron="";
my $new_wbbt_elements_ref;
my @new_wbbt_elements;
my %quantity_neuron_hash=();
my $quantity_neuron;
# For head motor neurons
$substitute_neuron = "WBbt\:0008610";
$quantity_neuron_hash{$substitute_neuron}="";
# search for existing element
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@head_motor_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For vnc neurons
 $substitute_neuron = "WBbt\:0005300";
 $quantity_neuron_hash{$substitute_neuron}="";
# search for existing element
 $skip = 0;
 $quantity_neuron_hash{$substitute_neuron}="";
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 } 
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@vnc_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
# For sensory neurons
$substitute_neuron = "WBbt\:0005759";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@sensory_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For interneurons
$substitute_neuron = "WBbt\:0005113";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@interneurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For motor_neurons
$substitute_neuron = "WBbt\:0005409";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@motor_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For mechanosensory neurons
$substitute_neuron = "WBbt\:0008431";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@mechanosensory_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For pharyngeal neurons
$substitute_neuron = "WBbt\:0005439";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@pharyngeal_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For pharyngeal motor neurons
$substitute_neuron = "WBbt\:0003677";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@pharyngeal_motor_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
   $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For pharyngeal interneurons
$substitute_neuron = "WBbt\:0003668";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@pharyngeal_interneurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
  $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For amphid neurons
$substitute_neuron = "WBbt\:0005394";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@amphid_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
  $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For outer labial neurons
$substitute_neuron = "WBbt\:0006801";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@outer_labial_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
  $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# For inner labial neurons
$substitute_neuron = "WBbt\:0005117";
$quantity_neuron_hash{$substitute_neuron}="";
 $skip = 0;
foreach my $element (@wbbt_elements){
 if ($element =~/$substitute_neuron/){
  $skip = 1;
  $quantity_neuron_hash{$substitute_neuron}="all";
 }
}
($quantity_neuron, $new_wbbt_elements_ref) =replace_neurons(\@wbbt_elements, \@inner_labial_neurons, $substitute_neuron);
   @new_wbbt_elements = @{$new_wbbt_elements_ref};
   @wbbt_elements = @new_wbbt_elements;
if ($skip ==0) {
 $quantity_neuron_hash{$substitute_neuron}=$quantity_neuron;
}
#
# Check @unique_anatomies; keep higher level terms and discard children.
# For each anatomy term, find the parents and the children
         foreach my $a (@wbbt_elements){
          my @c_array=();
          my @p_array=();
          chomp($a);
          my $ancestors = $parents{$a};
          my $child     = $children{$a};
          if ($child) {$child =~ s/\,$//; @c_array = split/\,/,$child;}
          if ($ancestors) {$ancestors =~ s/\,$//; @p_array = split/\,/,$ancestors;}
# Do not include "The Cell" as either a parent or a child
          foreach my $c (@c_array){
            if ($c !~/WBbt\:0004017/){
            push(@child_array, $c);
           }
          }
          foreach my $p (@p_array){
            if ($p !~/WBbt\:0004017/){
            push(@parent_array, $p);
           }
          }
        }
my @unique_child_array = uniq(@child_array);
my @unique_parent_array = uniq(@parent_array);
#
# Now use the new array with the lower redundant children removed.
#
foreach my $p (@unique_child_array){
my $pointer = 0;
while ($pointer <= $#wbbt_elements) {
#   print "\@wbbt_elements[$pointer] is $wbbt_elements[$pointer]\n";

   if ($wbbt_elements[$pointer] eq $p) {
#      print "parents must be removed.\n";
      splice(@wbbt_elements, $pointer, 1);
   }
   else {
      $pointer++;
   }
 }
}
         my @element_instance=();
         foreach my $element (@wbbt_elements){
          my $instance_of = $instance{$element};
          my @i = split(/\,/,$instance_of);
          foreach my $i (@i){
           push(@element_instance, $i);
          }
         }
         my @anatomy_instances=();
         my @unique_anatomy_instances=();
         my @unique_element_instance=uniq(@element_instance);
         foreach my $i (@unique_element_instance){
          my $a = $anatomy_ontology{$i};
          push(@anatomy_instances, $a); 
         }
          @unique_anatomy_instances=uniq(@anatomy_instances);
         my $Anatomy_Instance = join("\,", @unique_anatomy_instances);   

         foreach my $element (@wbbt_elements){
          my $a = $anatomy_ontology{$element};
          if ($a =~/neuron/){
           if (($a=~/I3/)or($a=~/I4/)or($a=~/I5/)or($a=~/I6/)or($a=~/M1/)or
               ($a=~/M4/)or($a=~/M5/)or($a=~/MI/)){
#            print "exception\t$a\n";
          } else {
            if (($element =~/WBbt\:0008610/) or ($element =~/WBbt\:0005300/) or ($element =~/WBbt\:0005759/) or 
                ($element =~/WBbt\:0005113/) or ($element =~/WBbt\:0005409/) or ($element =~/WBbt\:0008431/) or 
                ($element =~/WBbt\:0005439/) or ($element =~/WBbt\:0003677/) or ($element =~/WBbt\:0003668/) or 
                ($element =~/WBbt\:0005394/) or ($element =~/WBbt\:0006801/) or ($element =~/WBbt\:0005117/) ){
                $a = $quantity_neuron_hash{$element} . " " . $anatomy_ontology{$element};
#                print " a \=\t$element\t$quantity_neuron_hash{$element}\t$anatomy_ontology{$element}\n";
            }
            $a =~s/neuron/neurons/g;
          }
         }
           push(@anatomies, $a); 
        }
          my @unique_anatomies = uniq(@anatomies);
#
#
    my $Anatomy = join("\,", @unique_anatomies);   
    print "524 $gene_id\n";
    print "525 $name\n";
	print "526 $Anatomy\n";
	print "527 $Anatomy_Instance\n";
    my $summary = "$gene_id\t$name\t$Anatomy\n\t\t\t\t$Anatomy_Instance\n\n\n";
    my $gene_output = "$gene_id\n";

    my $size = @unique_anatomies;
    my $sentence = $name . $helper_verb;
    if ($size == 1){
      $sentence .= $unique_anatomies[0];
    if ($sentence =~/\bCell\b/){
        $sentence = $name . $helper_verb_widely;
    }
    } elsif ($size ==2){
      $sentence .= $unique_anatomies[0] . " and the " . $unique_anatomies[1];
      if ($sentence =~/\bCell\b/){
        my $not_cell = "";
        if ($unique_anatomies[0] =~/Cell/){
            $not_cell = $unique_anatomies[1];
        } else {
            $not_cell = $unique_anatomies[0];
        }
        $sentence = $name . $helper_verb_several . $not_cell;
      }
    } else {
       my $count = 0;
       foreach my $a (@unique_anatomies){
        $count++;
        my $index = $count-1;
        if ($count == $size){
         $sentence .= "\, and the " . $unique_anatomies[$index];
        } elsif ($count == 1) {
         $sentence .= $unique_anatomies[$index];
       } else {
         $sentence .= "\, " . $unique_anatomies[$index];
       }
      }
      if ($sentence =~/\bCell\b/){
          $sentence = $name . $helper_verb_several;
          $count = 0;
          my @no_cell = ();
       foreach my $a (@unique_anatomies){
          next if ($a =~/Cell/);
          push(@no_cell, $a);
       }
       my $new_size = @no_cell;
       if ($new_size > 2){
       foreach my $a (@no_cell){
        $count++;
        my $index = $count-1;
        if ($count == $new_size){
         $sentence .= "\, and the " . $no_cell[$index];
        } elsif ($count == 1) {
         $sentence .= $no_cell[$index];
        } else {
         $sentence .= "\, " . $no_cell[$index];
        }
       }
      } else {
      $sentence .= $no_cell[0] . " and the " . $no_cell[1];
      }
     }
    }
    $sentence .="\;\n";
    $sentence =~s/the all/all/g;
    $sentence =~s/the certain/certain/g;
    $sentence =~s/\, male\,/\, in the male\,/g;
    $sentence =~s/\, and the male\;/\, and in the male\;/g;
    $sentence =~s/and the male\;/and in the male\;/g;
    $sentence =~s/including the male\;/including in the male\;/g;
    $sentence =~s/\, female\,/\, in the female\,/g;
    $sentence =~s/\, and the female\;/\, and in the female\;/g;
    $sentence =~s/and the female\;/and in the female\;/g;
    $sentence =~s/including the female\;/including in the female\;/g;
    $sentence =~s/\, hermaphrodite\,/\, in the hermaphrodite\,/g;
    $sentence =~s/\, and the hermaphrodite\;/\, and in the hermaphrodite\;/g;
    $sentence =~s/and the hermaphrodite\;/and in the hermaphrodite\;/g;
    $sentence =~s/including the hermaphrodite\;/including in the hermaphrodite\;/g;

    my @sentence_terms = split(/\,/, $sentence);
	print $summary;
        print $sentence;
        print $gene_output; 

     write_file($output_test, {append => 1 }, $summary);
     write_file($output_test_genes, {append => 1 }, $gene_output);
     write_file($output_test_log, {append => 1 }, $gene_output);
     write_file($output_test_log, {append => 1 }, $sentence);
     write_file($output_test_log, {append => 1 }, "\n\n\n");

    if (@sentence_terms < 21) {
     write_file($output, {append => 1 }, $summary);
     write_file($output_genes, {append => 1 }, $gene_output);
     write_file($output_log, {append => 1 }, $gene_output);
     write_file($output_log, {append => 1 }, $sentence);
     write_file($output_log, {append => 1 }, "\n\n\n");
     write_file($gene_file, $sentence);
    }
 }
}
exit 0;
sub get_neurons{
my $input_file = shift;
my @neurons_name_id = read_file($input_file); 
my %neuron_id_name=();
my %neuron_name_id=();
my @neuron_list;
   foreach (@neurons_name_id){
     my $file_line = $_;
     chomp($file_line);
      my ($name, $id) = split(/\t/,$file_line);
      chomp($id);
      chomp($name);
      $name=~s/^\s+//; 
      $name=~s/\s+$//;
      $id=~s/^\s+//; 
      $id=~s/\s+$//;  
      $neuron_id_name{$id} = $name;
      $neuron_id_name{$name} = $id;
      push(@neuron_list, $id);
    }
return \@neuron_list;
}
sub replace_neurons{
my $wbbt_elements_array_ref = shift;
my $neuron_array_ref = shift;
my $substitute_neuron = shift;

my @neuron_array = @{$neuron_array_ref};
my $neuron_size = @neuron_array;
my @wbbt_elements= @{$wbbt_elements_array_ref};

my $neuron_quantifier="";
my $neuron_count = 0; 
        foreach my $element (@wbbt_elements){
          chomp($element);
          foreach my $neuron (@neuron_array){
            chomp($neuron);
            if ($element =~/$neuron/){
              $neuron_count++;
              $element = $substitute_neuron;
            }
          }
         }

if ($neuron_count > 0){
 if ($neuron_count < $neuron_size){
     $neuron_quantifier = "certain";
  } elsif ($neuron_count >= $neuron_size){
     $neuron_quantifier = "all";
  }
# print "neuron count\t$neuron_count\t$neuron_quantifier\t$substitute_neuron\n";
}
my @new_wbbt_elements = sort(uniq(@wbbt_elements));
return ($neuron_quantifier, \@new_wbbt_elements);

}
sub get_gene_wbbt_hash{
 my $file = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/Uncertain/);
   next if ($line =~/NOT/);
   next if ($line =~/Life\_stage/);
   next if ($line =~/WB\_REF\:WBPaper00040986/);
   next if ($line !~/taxon\:6239/);
   next if ($line !~/IDA/);
   next if ($line =~ /Enriched/); ##liyuling 20160418	
   #next if ($line !~/IEP/);  #liyuling
   chomp($line);
#   print "line is $line\n";
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1];
   my $gene_name = $fields[2];
   my $qualifier = $fields[3];
   my $wbbt = $fields[4];
   my $paper = $fields[5];
   $qualifier =~ s/^\s+//;
   $qualifier =~ s/\s+$//;
   $wbbt =~ s/^\s+//;
   $wbbt =~ s/\s+$//;
   $gene_id =~ s/^\s+//;
   $gene_id =~ s/\s+$//;
 
   if ($wbbt !~/WBbt/){
     print "ERROR\: $gene_id\t$gene_name\t$wbbt\n";
   }

    if ($hash{$gene_id}){
        $hash{$gene_id} .= "\," . $wbbt;
     } else {
        $hash{$gene_id} = $wbbt;
     }
#   print "hash is $gene_id\t$hash{$gene_id}\n";
  } 
 return \%hash;
}
