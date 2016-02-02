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
my $anatomy_ec_dir = $home . "anatomy_expression_cluster/";
my $elegans_gene_list_dir = $home_elegans . "gene_lists/" ;
my $gene_list_dir = $home . "gene_lists/" ;
my $input_path = $anatomy_ec_dir . "input_files/";
my $output_path = $anatomy_ec_dir . "output_files/";
my $gene_dir = $output_path . "individual_gene_sentences/";
my $input_file = $input_path . $prefix . "ECsummary_anatomy.$RELEASE.txt";

my $tissue_dir = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/tissue_expression/input_files/";
# create arrays for neurons
my @neurons = ();
my %neuron_names = ();

my $vnc = $tissue_dir . "ventral_cord_neuron.txt";
my $vnc_neurons_ref=get_neurons($vnc);
my @vnc_neurons=@{$vnc_neurons_ref};
foreach my $n (@vnc_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $interneuron = $tissue_dir . "interneuron.txt";
my $interneurons_ref=get_neurons($interneuron);
my @interneurons=@{$interneurons_ref};
foreach my $n (@interneurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $sensory = $tissue_dir . "sensory_neuron.txt";
my $sensory_neurons_ref=get_neurons($sensory);
my @sensory_neurons=@{$sensory_neurons_ref};
foreach my $n (@sensory_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $head_motor =  $tissue_dir . "head_neuron.txt";
my $head_motor_neurons_ref=get_neurons($head_motor);
my @head_motor_neurons=@{$head_motor_neurons_ref};
foreach my $n (@head_motor_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $motor = $tissue_dir . "motor_neuron.txt";
my $motor_neurons_ref=get_neurons($motor);
my @motor_neurons=@{$motor_neurons_ref};
foreach my $n (@motor_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $mechanosensory = $tissue_dir . "mechanosensory_neuron.txt";
my $mechanosensory_neurons_ref=get_neurons($mechanosensory);
my @mechanosensory_neurons=@{$mechanosensory_neurons_ref};
foreach my $n (@mechanosensory_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $inner_labial = $tissue_dir . "inner_labial_neuron.txt";
my $inner_labial_neurons_ref=get_neurons($inner_labial);
my @inner_labial_neurons=@{$inner_labial_neurons_ref};
foreach my $n (@inner_labial_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
   push(@neurons, $n);
   $neuron_names{$n} = $n;
 }
}

my $outer_labial = $tissue_dir . "outer_labial_neuron.txt";
my $outer_labial_neurons_ref=get_neurons($outer_labial);
my @outer_labial_neurons=@{$outer_labial_neurons_ref};
foreach my $n (@outer_labial_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
  push(@neurons, $n);
  $neuron_names{$n} = $n;
 }
}

my $pn = $tissue_dir . "pharyngeal_neuron.txt";
my $pn_ref=get_neurons($pn);
my @pn_neurons=@{$pn_ref};
foreach my $n (@pn_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
   push(@neurons, $n);
   $neuron_names{$n} = $n;
 }
}

my $pmn = $tissue_dir . "pharyngeal_motor_neuron.txt";
my $pharyngeal_motor_neurons_ref=get_neurons($pmn);
my @pharyngeal_motor_neurons=@{$pharyngeal_motor_neurons_ref};
foreach my $n (@pharyngeal_motor_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
   push(@neurons, $n);
   $neuron_names{$n} = $n;
 }
}

my $amphid = $tissue_dir . "amphid.txt";
my $amphid_neurons_ref=get_neurons($amphid);
my @amphid_neurons=@{$amphid_neurons_ref};
foreach my $n (@amphid_neurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
   push(@neurons, $n);
   $neuron_names{$n} = $n;
 }
}

my $pharyngeal_interneuron = $tissue_dir . "pharyngeal_interneuron.txt";
my $pharyngeal_interneurons_ref=get_neurons($pharyngeal_interneuron);
my @pharyngeal_interneurons=@{$pharyngeal_interneurons_ref};
foreach my $n (@pharyngeal_interneurons){
 $n =~s/^\s+//;
 $n =~s/\s+$//;
 if ($n){
   push(@neurons, $n);
   $neuron_names{$n} = $n;
 }
}


if (-e $input_file) {

my $input_file_size = -s $input_file;

if ($input_file_size > 1){

my @ec_lines = read_file($input_file); 
my %gene_relate=();
my %gene_exp=();
my %gene_anatomy=();

my @gene_list;
   foreach (@ec_lines){
     my $ec_line = $_;
     chomp($ec_line);
     next if ($ec_line =~/GeneID/);
      my ($no, $name, $relate, $anatomy, $exp) = split(/\t/,$ec_line);
      if ($no){
          $gene_relate{$no} = $relate;
          $gene_exp{$no} = $exp;
          $gene_anatomy{$no} = $anatomy;
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
my $output_log = $output_path . "sentences_for_anatomy_expression_cluster.txt";
# if the output file exists delete it.
if (-e $output_log){
   my @args = ("/bin/rm", "-f", $output_log);
   system(@args) == 0 or die("could not delete file $output_log\n");
}
my $output_test_log = $output_path . "test_sentences_for_anatomy_expression_cluster.txt";
# if the output file exists delete it.
if (-e $output_test_log){
   my @args = ("/bin/rm", "-f", $output_test_log);
   system(@args) == 0 or die("could not delete file $output_test_log\n");
}
my $output = $output_path . "gene_anatomy_expression_cluster.txt";
my $output_test = $output_path . "test_gene_anatomy_expression_cluster.txt";
my $output_genes = $output_path . "gene_ids_anatomy_expression_cluster.txt";
my $output_test_genes = $output_path . "test_gene_ids_anatomy_expression_cluster.txt";
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
# helper verb
my $helper_verb = " is expressed in the ";
my $helper_verb_widely = " is expressed widely";
my $helper_verb_several = " is expressed in several tissues including the ";
my $helper_verb_neuron  = " nervous system";
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
my $male_specific = "male\-specific";
my $hermaphrodite_specific = "hermaphrodite\-specific";
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
    my $herm = 0;
    my $male = 0;
    my $sentence = "";
    my $relation = lc $gene_relate{$gene_id};
    my $experiment = $gene_exp{$gene_id};
    my @exps = split(/\,/, $experiment);
    next if (@exps <= 0);
    my $experiments = "";
    my $studies_indicate = "";
    if (@exps==1){
        $experiments = $exps[0];
        $studies_indicate = "studies indicate";
    } elsif (@exps == 2){
        $experiments =  "$exps[0] and $exps[1]";
        $studies_indicate = "studies indicate";
    } elsif (@exps > 2) {
        $studies_indicate = "studies indicate";
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
    $experiments =~s/study//gi;
    $experiments =~s/RNA seq/RNA sequencing/g;
    $experiments =~s/RNA\-seq/RNA sequencing/g;
    $experiments = lc $experiments;
    $experiments =~s/rna sequencing/RNA sequencing/g;
    $experiments =~s/ \,/\,/g;

    my $anatomy = $gene_anatomy{$gene_id};
    if (($anatomy =~/$male_specific/) or ($anatomy =~/$hermaphrodite_specific/)){
     if (($anatomy =~/$male_specific/) and ($anatomy =~/$hermaphrodite_specific/)){
      $anatomy =~s/$male_specific//g;
      $anatomy =~s/$hermaphrodite_specific//g;
      $anatomy .= "\,hermaphrodite";
      $anatomy .= "\,$male_specific tissues";
     } elsif ($anatomy =~/$male_specific/) {
      $anatomy =~s/$male_specific//g;
      $anatomy .= "\,$male_specific tissues";
     } elsif ($anatomy =~/$hermaphrodite_specific/) {
      $anatomy =~s/$hermaphrodite_specific//g;
      $anatomy .= "\,$hermaphrodite_specific tissues";
     }
      $anatomy =~s/\,\,/\,/g;
      $anatomy =~s/\,$//g;
      $anatomy =~s/^\,//g;
    }
    my $tissues = "";
    my $neurons = "";
    my @no_neurons = ();
    my @nerve_cells = ();
    my @anatomies = split(/\,/,$anatomy);
    my $including=0;
    foreach my $a (@anatomies){
     if (($neuron_names{$a}) or ($a=~/neuron/i) or ($a=~/ganglion/i)){
       if ($a ne 'neuron'){
        push(@nerve_cells, $a);
       } else {
          $including = 1;
       }
     } else {
        push(@no_neurons, $a);
     }
    }
    my $neuron_count = 0;
    if (@nerve_cells > 0){
     if ((@nerve_cells==1) and (($nerve_cells[0] eq 'neuron') or ($nerve_cells[0] eq 'neurons'))){
         $neurons = "neurons";
         $neuron_count=1;
      } elsif (@nerve_cells==1){
           $neurons = $nerve_cells[0];
           $neurons =~s/neuron//gi;
           $neuron_count=2;
      } elsif (@nerve_cells == 2){
           $neuron_count=2;
         if (($nerve_cells[0] ne 'neuron') and ($nerve_cells[1] ne 'neuron')){
           $neurons =  "$nerve_cells[0] and $nerve_cells[1]";
           $neurons =~s/neuron//gi;
         } elsif ($nerve_cells[0] ne 'neuron') {
          $neurons =  $nerve_cells[0];
          $neurons =~s/neuron//gi;
         } elsif ($nerve_cells[1] ne 'neuron') {
          $neurons =  $nerve_cells[1];
          $neurons =~s/neuron//gi;
         }
      } elsif (@nerve_cells > 2) {
        my $c = 0;
        foreach my $n (@nerve_cells){
         $c++;
         if ($c < (@nerve_cells)){
          $neurons .= "$n\, ";
         } elsif ($c==@nerve_cells) {
          $neurons .= "and $n";
         }
        }
         $neuron_count = $c;
         $neurons =~s/neuron//gi;
       }
      }   
         $neurons =~s/ \,/\,/gi;         
        if (@no_neurons==1){
        $tissues = $no_neurons[0];
    } elsif (@no_neurons == 2){
        $tissues =  "$no_neurons[0] and $no_neurons[1]";
    } else {
        my $c = 0;
        foreach my $a (@no_neurons){
         $c++;
         if ($c < @no_neurons){
          $tissues .= "$a\, ";
         } elsif ($c==@no_neurons) {
          $tissues .= "and $a";
         }
        }
    }
   
     my @gs = ();
     my @ns = (); 
     $neurons =~s/'\, and '/'\,'/gi;
     $neurons =~s/' and '/'\,'/gi;
     $neurons =~s/and/\,/gi;
     my @nervous = split(/\,/, $neurons);
      foreach my $n (@nervous){
          $n =~ s/^\s+//;
          $n =~ s/\s+$//;
       if ($n =~/ganglion/){
        push(@gs, $n);
       } else {
        push(@ns, $n);
       }
      }
    my $nerve = "";
    my @nerves = ();
    if ((@ns > 0) and (@gs > 0)){
       @nerves = (@ns, @gs);
    } elsif (@ns > 0){
       @nerves = (@ns);
    } elsif (@gs > 0){
       @nerves = (@gs);
    }
    if (@nerves==1){
        $nerve = $nerves[0];
        $nerve =~ s/^\s+//;
        $nerve =~ s/\s+$//;
    } elsif (@nerves == 2){
        $nerve =  "$nerves[0] and $nerves[1]";
        $nerve =~ s/^\s+//;
        $nerve =~ s/\s+$//;
    } elsif (@nerves > 2) {
        my $c = 0;
        foreach my $n (@nerves){
        $n =~ s/^\s+//;
        $n =~ s/\s+$//;
         $c++;
         if ($c < @nerves){
          $nerve .= "$n\, ";
         } elsif ($c==@nerves) {
          $nerve .= "and $n";
         }
        }
    }
    if ($nerve){
        $nerve =~s/^\s+//;
        $nerve =~s/\s+$//;
        $nerve =~s/\, \,/\,/g;
#    print "nerve is $nerve\n";
    }
    if ((@anatomies==1) and ($anatomy =~/neuron/)){
     $sentence = "$experiments $studies_indicate that $name is enriched in neurons";
    } elsif ((@no_neurons > 0) and (@nerve_cells==0)) {
     if ($anatomy=~/neuron/){
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues and in neurons"; 
     } else {
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues"; 
     }
    } elsif ((@no_neurons==0) and (@nerve_cells > 0)){
         if (($nerve=~/ganglion/) and ($including==1)){
           $sentence = "$experiments $studies_indicate that $name is enriched in the neurons including the $nerve";
         } elsif ($including==1) {
           $sentence = "$experiments $studies_indicate that $name is enriched in neurons including the $nerve neurons";
         } elsif ($nerve =~/ganglion/) {
           $sentence = "$experiments $studies_indicate that $name is enriched in the $nerve";
         } else {
           $sentence = "$experiments $studies_indicate that $name is enriched in the $nerve neurons";
         }
    }  elsif ((@no_neurons > 0) and (@nerve_cells > 0)){
         if (($nerve=~/ganglion/) and ($including==1)) {
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues and in neurons including $nerve"; 
         } elsif ($including==1) {
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues and in neurons including the $nerve neurons"; 
         } elsif ($nerve=~/ganglion/) {
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues and in the $nerve"; 
         } else {
     $sentence = "$experiments $studies_indicate that $name is enriched in the $tissues and in the $nerve neurons"; 
         }
    }
    $sentence =~s/bodywall/body wall/gi;
    $sentence =~s/coelomocyte/coelomocytes/gi;
    $sentence .="\;\n";
    $sentence =~s/the all/all/g;
    $sentence =~s/the certain/certain/g;
    $sentence =~s/  / /g;

#    $sentence = ucfirst($sentence);

    my @sentence_terms = split(/\,/, $anatomy);

#     write_file($output_test, {append => 1 }, $summary);
#     write_file($output_test_genes, {append => 1 }, $gene_output);
#     write_file($output_test_log, {append => 1 }, $gene_output);
     write_file($output_test_log, {append => 1 }, $sentence);
     write_file($output_test_log, {append => 1 }, "\n\n\n");

    if (@sentence_terms < 21) {
#     write_file($output, {append => 1 }, $summary);
#     write_file($output_genes, {append => 1 }, $gene_output);
#     write_file($output_log, {append => 1 }, $gene_output);
     write_file($output_log, {append => 1 }, $sentence);
     write_file($output_log, {append => 1 }, "\n\n\n");
     write_file($gene_file, $sentence);
    }
   }
  }
 }
}
exit 0;
sub get_neurons{
my $input_file = shift;
my @neurons_name_id = read_file($input_file); 
my %neuron_id_name=();
my %neuron_name_id=();
my @neuron_list=();
   foreach (@neurons_name_id){
     my $file_line = $_;
     chomp($file_line);
      my ($name, $id) = split(/\t/,$file_line);
      $name=~s/^\s+//; 
      $name=~s/\s+$//;
      $id=~s/^\s+//; 
      $id=~s/\s+$//;
      $neuron_id_name{$id} = $name;
      $neuron_id_name{$name} = $id;
      push(@neuron_list, $name);
    }
foreach my $n (@neuron_list){
      my $out = "$n\n";
      my $neuron_log = "./neuron.log";
      write_file($neuron_log, {append => 1 }, $out);
}
return \@neuron_list;
}
