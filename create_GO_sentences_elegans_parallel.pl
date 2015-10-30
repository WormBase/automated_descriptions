#!/usr/bin/env perl
use ConciseDescriptions;
use LWP::Simple;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use Switch;
use warnings;
use strict;
my $AND = "AND";
my $species_term = $ARGV[0];
my ($species, $project, $species_name, $species_prefix) = split(/$AND/, $species_term);

chomp($species);
chomp($project);
chomp($species_name);

 $species_name =~s/\_/ /g;

if ($species_prefix){
 chomp($species_prefix);
 $species_prefix =~ s/^\s+//;
 $species_prefix =~ s/\s+$//;
}

$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;
$species_name =~ s/^\s+//;
$species_name =~ s/\s+$//;

if ($species =~/elegans/){
#
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $gene_list_dir = $home . "gene_lists/";
# Define GO term file
my $ontology = "GO";
my $go_file = $gene_list_dir . "go_terms.txt";
my $go_altid =  $gene_list_dir . "goid_altid.txt";
my $input_ace_file = $gene_list_dir . "go_terms.ace";
my $db_gene_list  = $gene_list_dir . "wb_gene_list.txt";
# Define component adjectives
my @component_adjectives = ("GO\:0005622");
#
# Create a file containing the GO terms and their respective parents
# build parent/child hash for granularity
my ($parents_ref, $children_ref) =  ConciseDescriptions::get_ontology_parents_children($input_ace_file, $ontology);
my %parents = %$parents_ref;
my %children = %$children_ref;
# Create GO hash so that terms can be referenced by GO ID.
my $gene_ontology_ref = ConciseDescriptions::get_ontology_hash($ontology, $go_file, $go_altid);
my %gene_ontology = %$gene_ontology_ref;
my $go_elegans_dir = $home . "gene_ontology/input_files/";
my $gene_association_elegans_file = $go_elegans_dir . "gene_association.wb";
my $elegans_orthology = $home . "orthology/";
my $go_path = $home . "gene_ontology/output_files/";
my $output_process_file   = $go_path . "sentences_for_process.txt";
my $output_function_file  = $go_path . "sentences_for_function.txt";
my $output_component_file = $go_path . "sentences_for_component.txt";
#
my $wbgene_elegans_name_hash_ref = get_wbgene_name_elegans_hash($gene_association_elegans_file, $db_gene_list);
my %wbgene_elegans_name_hash = %$wbgene_elegans_name_hash_ref;
#
my $doublecomma = "\,\,";
my $the_the = "the the";
my $the = "the";
my $is_is = "is is";
my $is = "is";
my $is_a_is_a = "is a is a";
my $is_a = "is a";
my $blank_comma = "\, \,";
my $comma = "\, ";
my $nadp =" NAD or NADP as acceptor\,";
my $component_watch_string = "integral component of";
my $component_watch_string_the = "integral component of the";
my $structural="structural constituent";
my $is_structural="is a structural constituent";
#
my $synaptic_1_0 = "synaptic transmission\, GABAergic";
my $synaptic_1_1 = "GABAergic synaptic transmission";

my $synaptic_2_0 = "synaptic transmission\, cholinergic";
my $synaptic_2_1 = "cholinergic synaptic transmission";

my $synaptic_3_0 = "synaptic transmission\, dopaminergic";
my $synaptic_3_1 = "dopaminergic synaptic transmission";

my $synaptic_4_0 = "synaptic transmission\, glutamatergic";
my $synaptic_4_1 = "glutamatergic synaptic transmission";

my $synaptic_5_0 = "synaptic transmission\, glycinergic";
my $synaptic_5_1 = "glycinergic synaptic transmission";

my $molting =   "the molting cycle";
my $molting_1 = "molting cycle\, chitin\-based cuticle";
my $molting_2 = "molting cycle\, collagen and cuticulin\-based cuticle";
my $molting_3 = "molting cycle\, protein\-based cuticle";
my $molting_0 = "molting cycle";

my $the_cell_0 = "localized to the cell\;";
my $the_cell_1 = "expressed widely\;";

my $embryo   = "embryo development";
my $embryo_1 = "embryo development ending in birth or egg hatching";

my $growth_0 = "multicellular organism growth";
my $growth_1 = "growth";
#
if (-e $go_path){
   my $go_path_flag = 1;
} else {
   mkdir $go_path or die "could not create $go_path";
}
my $individual_process_path = $go_path . "individual_process_sentences/";
if (-e $individual_process_path){
   my $individual_process_path_flag = 1;
} else {
   mkdir $individual_process_path or die "could not create $individual_process_path";
}
my $individual_component_path = $go_path . "individual_component_sentences/";
if (-e $individual_component_path){
   my $individual_component_path_flag = 1;
} else {
   mkdir $individual_component_path or die "could not create $individual_component_path";
}
my $individual_function_path = $go_path . "individual_function_sentences/";
if (-e $individual_function_path){
   my $individual_function_path_flag = 1;
} else {
   mkdir $individual_function_path or die "could not create $individual_function_path";
}
# if the output file exists delete it.
if (-e $output_function_file){
   my @args = ("/bin/rm", "-f", $output_function_file);
   system(@args) == 0 or die("could not delete file $output_function_file\n");
}
# if the output file exists delete it.
if (-e $output_component_file){
   my @args = ("/bin/rm", "-f", $output_component_file);
   system(@args) == 0 or die("could not delete file $output_component_file\n");
}
# if the output file exists delete it.
if (-e $output_process_file){
   my @args = ("/bin/rm", "-f", $output_process_file);
   system(@args) == 0 or die("could not delete file $output_process_file\n");
}
# if the output files for individual gene sentences exist delete them.
   my @individual_process_files = glob("$individual_process_path/WBGene*");
      foreach my $individual (@individual_process_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
# if the output files for individual gene sentences exist delete them.
   my @individual_component_files = glob("$individual_component_path/WBGene*");
      foreach my $individual (@individual_component_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
# if the output files for individual gene sentences exist delete them.
   my @individual_function_files = glob("$individual_function_path/WBGene*");
      foreach my $individual (@individual_function_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }
my $dead_gene_list = $gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);
my $gene_array_ref = get_gene_array($gene_association_elegans_file);
my @gene_array = @$gene_array_ref;

my @live_gene_array = ();
foreach my $test (@gene_array){
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
my @uncurated_genes_array = ();
foreach my $test (@sorted_live_gene_array){
 my $keep =0;
# foreach my $curated (@curated_genes){
#   if ($curated =~/$test/){
#    $keep = 1;
#   }
# }
 if ($keep ==0){
   push(@uncurated_genes_array, $test);
 }
}
my @sorted_uncurated_genes_array = sort(@uncurated_genes_array); 

my $gene_name_hash_ref = get_gene_name_hash($gene_association_elegans_file, $db_gene_list);
my %gene_name_hash = %$gene_name_hash_ref;
my $go = "P";
my $elegans_process_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_process_hash = %$elegans_process_hash_ref;
   $go = "C";
my $elegans_component_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_component_hash = %$elegans_component_hash_ref;
   $go = "F";
my $elegans_function_hash_ref = get_gene_go_hash($gene_association_elegans_file, $go, 1);
my %elegans_function_hash = %$elegans_function_hash_ref;

foreach my $gene_id (@sorted_uncurated_genes_array){
 chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;
 my @gene_processes = ();
 my @gene_functions = ();
 my @gene_components = ();

 my @GO_processes  = ();
 my @GO_functions  = ();
 my @GO_components = ();
 
 my $function_string = "";
 my $component_string = "";
 my $process_string = "";

 my $gene_name = $gene_name_hash{$gene_id};
 
 my $gene_process = $elegans_process_hash{$gene_id};
 if ($gene_process){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_process;
# check for granularity
 my $processes_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
 my @processes = uniq(sort(@{$processes_ref}));
 if (@processes){
# convert the ids into terms and put them in the array
    my $process_count = 0;
    foreach my $process (@processes){
     next if ($process =~/\[IEA\]/);
     next if ($process =~/\[ISS\]/);
     next if ($process =~/\[IBA\]/);
     next if ($process =~/\[IBD\]/);
     next if ($process =~/$gene_id/);
     next if ($process =~/$gene_name/);


     $process_count++;
     if ($process_count > 0){
        $process =~ s/^\s+//;
        $process =~ s/\s+$//;
     my $go_term = "";
             my ($evidence) = $process =~ /\[(\w+)\]/;
             if ($evidence) {
                 $process =~ s/\[$evidence\]//g;
             }
             $process =~ s/^\s+//;
             $process =~ s/\s+$//;
           
     $go_term = ConciseDescriptions::get_ontology_term($process, \%gene_ontology);
     if ($go_term){
      if ( lc ($go_term) !~/obsolete/) {
       if (grep {$_ =~/\Q$go_term\E/} @GO_processes) {
        print "$go_term ignored\; already in array\n";
       } else {
         $go_term .= "\[$evidence\]";
         push(@GO_processes, $go_term);
        }
       }
      }
     }
    }
     if (@GO_processes){
         @GO_processes = uniq(sort(@GO_processes));
#
# yet untested, $imp="IMP" would identify processes with evidence code separately; else it includes them without 
# adding information that it was "based on mutant phenotype"
#
     my $imp = "";
     my $process_goterm = get_verb_process_goterm(\@GO_processes, $imp);
     if ($process_goterm) {   
        $process_string = $process_goterm;
        $process_string =~ s/^\s+//;
        $process_string =~ s/\s+$//;
        $process_string =~ s/ \,/\,/g;
        $process_string =~ s/ +/ /g;
#
# Remove evidence codes
#
#        my ($evidence_to_remove) = $process_string =~ /\[(.*)\]/;
#        my ($evidence_to_remove) = $process_string =~ /\[(\w+)\]/;
#        if ($evidence_to_remove) {
#             $process_string =~ s/\[[^\]]*\]//g;
#           }
         my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
# Remove evidence codes
         foreach my $ec (@evidences){
	    $process_string =~ s/\[$ec\]//g;
         }

         $process_string =~ s/$molting_0/$molting/g;
         $process_string =~ s/$molting_1/$molting/g;
         $process_string =~ s/$molting_2/$molting/g;
         $process_string =~ s/$molting_3/$molting/g;

         $process_string =~ s/$embryo_1/$embryo/g;

         $process_string =~ s/$the_the/$the/g;
         $process_string =~ s/$blank_comma/$comma/g;

        $process_string =~ s/$synaptic_1_0/$synaptic_1_1/g;
        $process_string =~ s/$synaptic_2_0/$synaptic_2_1/g;
        $process_string =~ s/$synaptic_3_0/$synaptic_3_1/g;
        $process_string =~ s/$synaptic_4_0/$synaptic_4_1/g;
        $process_string =~ s/$synaptic_5_0/$synaptic_5_1/g;
        $process_string =~ s/ +/ /g;
     }# process_goterm
    } # GO_processes
 } # @processes
} # gene_process
 my $gene_function = $elegans_function_hash{$gene_id};
 if ($gene_function){
 my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_function;
# check for granularity
 my $functions_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "", "", \%parents);
 my @functions = uniq(sort(@{$functions_ref}));
 my $exp_functions_ref = get_exp_go_array(\@functions);
 my @exp_functions=uniq(sort(@{$exp_functions_ref}));
 my $iss_functions_ref = get_stat_go_array(\@functions);
 my @iss_functions=uniq(sort(@{$iss_functions_ref}));
 my $iea_functions_ref = get_iea_go_array(\@functions);
 my @iea_functions=uniq(sort(@{$iea_functions_ref}));
 my $iea_function_string="";
 my $iss_function_string="";
 my $exp_function_string="";
 
 if (@iea_functions){
     @iea_functions = sort(uniq(@iea_functions));
      my $ieaf_string="";
      foreach my $ieaf (@iea_functions){
       $ieaf_string .= "$ieaf\t";
     }
     $iea_function_string = get_function_string($iea_functions_ref, "IEA", \%gene_ontology);
 }
 if (@iss_functions){
      @iss_functions = sort(uniq(@iss_functions));
      my $issf_string="";
      foreach my $issf (@iss_functions){
       $issf_string .= "$issf\t";
     }
     $iss_function_string = get_function_string($iss_functions_ref, "ISS", \%gene_ontology);
 }
 if (@exp_functions){
     @exp_functions = sort(uniq(@exp_functions));
     $exp_function_string = get_function_string($exp_functions_ref, "EXP", \%gene_ontology);
 }
 foreach my $function (@iea_functions){
    if (grep {$_ =~/\Q$function\E/} @GO_functions) {
        print "$function ignored\; already in array\n";
      } else {
        push (@GO_functions, $function);
     }
 }
 foreach my $function (@iss_functions){
    if (grep {$_ =~/\Q$function\E/} @GO_functions) {
        print "$function ignored\; already in array\n";
      } else {
        push (@GO_functions, $function);
     }
 }
 foreach my $function (@exp_functions){
    if (grep {$_ =~/\Q$function\E/} @GO_functions) {
        print "$function ignored\; already in array\n";
      } else {
        push (@GO_functions, $function);
     }
 }
 if ((@exp_functions) and (@iss_functions) and (@iea_functions)){
     $function_string = $exp_function_string . " and " . $iss_function_string . " and " . $iea_function_string;
 } elsif ((@exp_functions) and (@iss_functions)) {
     $function_string = $exp_function_string . " and " . $iss_function_string;
 } elsif ((@exp_functions) and (@iea_functions)) {
     $function_string = $exp_function_string . " and " . $iea_function_string;
 } elsif ((@iss_functions) and (@iea_functions)) {
     $function_string = $iss_function_string . " and " . $iea_function_string;
 } elsif (@exp_functions) {
     $function_string = $exp_function_string;
 } elsif (@iss_functions) {
     $function_string = $iss_function_string;
 } elsif (@iea_functions){
     $function_string = $iea_function_string;
 }
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
      $function_string =~ s/ +/ /g;
      $function_string =~ s/$the_the/$the/g;
      $function_string =~ s/$blank_comma/$comma/g;
      $function_string =~ s/$nadp//g;
      $function_string =~ s/$doublecomma/$comma/g;     
      $function_string =~ s/ +/ /g;

} # if gene function string
 my $gene_component = "";
 my @gene_component_results = ();
 my $gene_component_result = $elegans_component_hash{$gene_id};
 if ($gene_component){
    @gene_component_results = split(/\,/, $gene_component_result);
 }
 my @unique_component_results = uniq(@gene_component_results);
 if (@unique_component_results > 1){
  my @c=();
  foreach my $c (@unique_component_results){
   if ($c !~/GO\:0005623/){ 
     push(@c, $c);
   }
  }
   $gene_component = join(",",@c);
 } else {
   $gene_component = $gene_component_result;
 }
 if ($gene_component){
   my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
  my $gene_row = $gene_id . "\," . $gene_name . "\," . $gene_component;
 if ($gene_component){
    @gene_components = uniq(split(/\,/, $gene_component));
 }
   
   my $components_ref = ConciseDescriptions::get_granular_array($ontology, $gene_row, "IEA", "ISS", \%parents);
   my @components = uniq(sort(@{$components_ref}));
#
# Check if any components are adjectives 
#
   my @new_component_adjectives=();
   if (@components){
     my $pointer = 0;
      while ($pointer <= $#components) {
       my $component = $components[$pointer];
       foreach my $adjective (@component_adjectives){
        if ($component =~/$adjective/){
           push(@new_component_adjectives, $component);
           splice(@components, $pointer, 1);
        } else {
           $pointer++;
        } 
       } 
     }
   }
 my $adj_phrase = " is ";
 my $adj_sentence;
 my $adj_count = @new_component_adjectives; 

 if ($adj_count > 0){
   my $count = 0;
   $adj_sentence = $adj_phrase;

   foreach my $adj (@new_component_adjectives){
#             my ($evidence) = $adj =~ /\[(.*)\]/;

#             my ($evidence) = $adj =~ /\[(\w+)\]/;
#             if ($evidence) {
#                 $adj =~ s/\[$evidence\]//g;
#             }
# Remove evidence codes
         foreach my $ec (@evidences){
	    $adj =~ s/\[$ec\]//g;
         }
             $adj =~ s/^\s+//;
             $adj =~ s/\s+$//;
     my $go_term = "";
     $go_term = ConciseDescriptions::get_ontology_term($adj, \%gene_ontology);
     next if ($go_term eq "");
     next if (lc($go_term) =~/obsolete/);
     my $count++;
     if (($count == 1) and ($adj_count <= 2)){
      $adj_sentence .= $go_term;
     } elsif (($count < ($adj_count-1)) and ($adj_count >= 3)){
       $adj_sentence .= $go_term . "\,";
     } elsif (($count == $adj_count) and ($adj_count >= 2)){
       $adj_sentence .= " and " . $go_term;
     } # if count
  } # foreach adj
 } # if adj_count
   my $component_count = 0;
   if (@components) {
    foreach my $component (@components){
    $component_count++;
    if ($component_count > 0){
     next if ($component =~/$gene_id/);
     next if ($component =~/$gene_name/);
     next if ($component =~/\[IEA\]/);
     next if ($component =~/\[ISS\]/);
     next if ($component =~/\[IBD\]/);
     next if ($component =~/\[IBA\]/);
#             my ($evidence) = $component =~ /\[(.*)\]/;
#             my ($evidence) = $component =~ /\[(\w+)\]/;
#             if ($evidence) {
#                 $component =~ s/\[$evidence\]//g;
#             }
         my $evidence;
         foreach my $ec (@evidences){
            if ($component =~/\[$ec\]/){
                $evidence = $ec;
	        $component =~ s/\[$ec\]//g;
            }
         }
             $component =~ s/^\s+//;
             $component =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($component, \%gene_ontology);
     if ($go_term){
      if ( lc ($go_term) !~/obsolete/){
        $go_term .="\[$evidence\]";
      if (grep {$_ =~/\Q$go_term\E/} @GO_components) {
        print "$go_term ignored\; already in array\n";
      } else {
        push(@GO_components, $go_term);
      }
     }
    }# 
   } # if component count
  } # foreach component
       if (@GO_components){
        @GO_components = sort(uniq(@GO_components));
        $component_string = get_verb_component_goterm(\@GO_components);
         if ($adj_sentence){
          $component_string .= " and " . $adj_sentence;
         }
        $component_string =~ s/$structural/$is_structural/g;
        $component_string =~ s/$component_watch_string/$component_watch_string_the/g;
        $component_string =~ s/$the_the/$the/g;
        $component_string =~ s/$blank_comma/$comma/g;

my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
# Remove evidence codes
         foreach my $ec (@evidences){
	    $component_string =~ s/\[$ec\]//g;
         }
        $component_string =~ s/ +/ /g;
        $component_string =~ s/^\s+//;
        $component_string =~ s/\s+$//;
        $component_string =~ s/ \,/\,/g;
        $component_string =~ s/ +/ /g;
      }
   } # if components array
 }
#  my $sentence = "";
#
#  if ((@GO_components) and (@GO_processes) and (@GO_functions)){
#    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $function_string . "\; " . $gene_name . " " . $component_string;
#  } elsif ((@GO_processes) and (@GO_functions)){
#    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $function_string;
#  } elsif ((@GO_components) and (@GO_functions)){
#    $sentence = $gene_name . " " . $function_string . "\; " . $gene_name . " " . $component_string;
#  } elsif ((@GO_components) and (@GO_processes)){
#    $sentence = $gene_name . " " . $process_string . "\; " . $gene_name . " " . $component_string;
#  } elsif (@GO_components){
#    $sentence = $gene_name . " " . $component_string;
#  } elsif (@GO_processes){
#    $sentence = $gene_name . " " . $process_string;
#  } elsif (@GO_functions){
#    $sentence = $gene_name . " " . $function_string;
#  } else {
#    $sentence = "";
#  }
  if (($process_string) or ($function_string) or ($component_string)){
          my $doublespace = "  ";
          my $space = " ";
# 
          $component_string =~s/$growth_0/$growth_1/g;
          $component_string =~s/$structural/$is_structural/g;
          $function_string =~s/$growth_0/$growth_1/g;
          $function_string =~s/$structural/$is_structural/g;
          $process_string =~s/$growth_0/$growth_1/g;
          $process_string =~s/$structural/$is_structural/g;
# remove redundant is
          $function_string  =~s/$is_is/$is/g;
          $function_string  =~s/$is_a_is_a/$is_a/g;
          $process_string   =~s/$is_is/$is/g;
          $process_string   =~s/$is_a_is_a/$is_a/g;
          $component_string =~s/$is_is/$is/g;
          $component_string =~s/$is_a_is_a/$is_a/g;
# remove any double spaces
          $component_string =~s/$doublespace/$space/g;
          $function_string  =~s/$doublespace/$space/g;
          $process_string   =~s/$doublespace/$space/g;
# add semi-colon at end
      if (length($process_string) > 0){
          $process_string .= "\;";
        my $out = $individual_process_path . $gene_id;
        my $sentence = $gene_name . " " . $process_string;
           $sentence =~s/$doublespace/$space/g;
        write_file($out, $sentence);
        write_file($output_process_file, {append => 1 }, $gene_id);
        write_file($output_process_file, {append => 1 }, "\n");
        write_file($output_process_file, {append => 1 }, $sentence);
        write_file($output_process_file, {append => 1 }, "\n\n\n");
         }
      if (length($component_string) > 0){
          $component_string .= "\;";
          if ($component_string=~/$the_cell_0/){
              my $location = index ( $component_string, $the_cell_0 );
              if ($location != -1){
                $component_string =~s/$the_cell_0/$the_cell_1/g;
              }
          }
        my $out = $individual_component_path . $gene_id;
        my $sentence = $gene_name . " " . $component_string;
           $sentence =~s/$doublespace/$space/g;
        write_file($out, $sentence);
        write_file($output_component_file, {append => 1 }, $gene_id);
        write_file($output_component_file, {append => 1 }, "\n");
        write_file($output_component_file, {append => 1 }, $sentence);
        write_file($output_component_file, {append => 1 }, "\n\n\n");
         }
      if (length($function_string) > 0){
          $function_string .= "\;";
        my $sentence = $gene_name . " " . $function_string;
        my $out = $individual_function_path . $gene_id;
           $sentence =~s/$doublespace/$space/g;
        write_file($out, $sentence);
        write_file($output_function_file, {append => 1 }, $gene_id);
        write_file($output_function_file, {append => 1 }, "\n");
        write_file($output_function_file, {append => 1 }, $sentence);
        write_file($output_function_file, {append => 1 }, "\n\n\n");
         }
  }
 }
} # if elegans
exit(0);

sub get_wbgene_name_elegans_hash{
 my $file = shift;
 my $db_gene_list = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[2];
  my $gene_id = $fields[1];
  if (($gene_name) and ($gene_id)){
       $hash{$gene_name} = $gene_id;
   }
 }
 my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);
     if (($WBGene) and ($Gene)){ 
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        chomp($Gene);
        $hash{$Gene} = $WBGene;
       }
     }
 }
 return \%hash;
}
sub get_gene_name_hash{
 my $file = shift;
 my $db_gene_list = shift;
 my %hash=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   chomp($line);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1];
   my $gene_name = $fields[2];
   if (($gene_name) and ($gene_id)){
      $gene_name =~ s/^\s+//;
      $gene_name =~ s/\s+$//;
      $gene_id   =~ s/^\s+//;
      $gene_id   =~ s/\s+$//;
      $hash{$gene_id} = $gene_name;
  }
 }
 my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);         
     if (($WBGene) and ($Gene)){
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $hash{$WBGene} = $Gene;
      }
     }
 }
 return \%hash;
}
sub get_gene_go_hash{
 my $file = shift;
 my $go = shift;
 my $ec = shift;
 my %hash=();
 
 $go =~ s/^\s+//;
 $go =~ s/\s+$//;
 my @exp_ec = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP");
 my @other_ec = ();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   next if ($line =~/UniProt/);
   next if ($line =~/IBA/);
   next if ($line =~/IBD/);
   next if ($line =~/GO\:0003674/); # veto molecular_function
   next if ($line =~/GO\:0005554/); # veto molecular_function alt_id
   next if ($line =~/GO\:0005575/); # veto cellular_component
   next if ($line =~/GO\:0008372/); # veto cellular_component alt_id
   next if ($line =~/GO\:0008150/); # veto biological_process
   next if ($line =~/GO\:0000004/); # veto biological_process alt_id
   next if ($line =~/GO\:0007582/); # veto biological_process alt_id
   next if ($line =~/GO\:0005488/); # veto binding
   next if ($line =~/GO\:0005515/); # veto protein binding
   chomp($line);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1]; 
   my $go_term = $fields[4]; 
   my $term;
   $go_term =~ s/^\s+//;
   $go_term =~ s/\s+$//;
   $gene_id =~ s/^\s+//;
   $gene_id =~ s/\s+$//;
   my $evidence_code = $fields[6]; 
   my $with = $fields[7];
   my $GO = $fields[8]; 
   if (($with =~/WBPhenotype/g) and ($evidence_code=~/IEA/) and ($GO =~/P/) ) {
         $evidence_code = "IMP";
   }
   if ($go_term){
   if ($go =~/$GO/){
    if ($hash{$gene_id}){
     if ($hash{$gene_id}!~/$go_term/g){
       if (($evidence_code) and ($ec)) {
         $term = $go_term . "\[" . $evidence_code . "\]";
       } else {
         $term = $go_term;
       }
        if ($hash{$gene_id} !~/$term/g){
            $hash{$gene_id} .= "\, " . $term;
        }
     } elsif ((grep {$_ =~ $evidence_code} @exp_ec) and ($ec)) {
         my $temp_string = $hash{$gene_id};
         my ($evidence) = $temp_string =~ /$go_term\[(\w+)\]/;
         if (grep {$_ =~ $evidence} @exp_ec) {
          print "$gene_id\t$go_term\t$evidence\t$hash{$gene_id}\n";
         } elsif ($evidence !~ $evidence_code) {
          my $old =  $go_term . "\[" . $evidence . "\]";
          my $replace = $go_term . "\[" . $evidence_code . "\]";
            $temp_string =~ s/$old/$replace/g;
            $hash{$gene_id} = $temp_string;
         }
     }
    } else {
       if (($evidence_code) and ($ec)) {
         $term = $go_term . "\[" . $evidence_code . "\]";
       } else{
         $term = $go_term;
       }
        $hash{$gene_id} = $term;
    }
   }
  }
 } 
 return \%hash;
}
sub get_verb_process_goterm{
    my $goterms_ref = shift;
    my $imp = shift;
    my $imp_string = " based on mutant phenotype";
    my @goterms = uniq(sort(@{$goterms_ref}));
    my $verb_goterm="";
# 
# The following was added (untested) to incorporate an IMP, "based on mutant phenotype".
# cluster the GO terms
#
   if ($imp !~/IMP/){
       $verb_goterm = get_verb_goterm_process(\@goterms);
 } else {     
    my @IMP=();
    my @other_ec =();
    if (@goterms){
       foreach my $goterm (@goterms){
         if ($goterm =~/\[IMP\]/){
          if (grep {$_ =~/\Q$goterm\E/} @IMP) {
              print "$goterm ignored\; already in array\n";
          } else {
              push(@IMP, $goterm);
          }
         } else {
            if (grep {$_ =~/\Q$goterm\E/} @other_ec) {
                print "$goterm ignored\; already in array\n";
          } else {
                push(@other_ec, $goterm);
          }
         }
      } 
    }
 my $phrase_imp = "";
 if (@IMP){
   $phrase_imp = get_verb_goterm_process(\@IMP);
   $phrase_imp .= $imp_string;
 }
 my $phrase_other_ec="";
 if (@other_ec){
   $phrase_other_ec = get_verb_goterm_process(\@other_ec);
 }
 if (($phrase_imp) and ($phrase_other_ec)){
    $verb_goterm = $phrase_other_ec . " and " . $phrase_imp;
 } elsif ($phrase_imp){
    $verb_goterm = $phrase_imp;
 } elsif ($phrase_other_ec){
    $verb_goterm = $phrase_other_ec;
  }
 } 
 return $verb_goterm;
}
sub get_verb_goterm_process{
    my $goterms_ref = shift;
    my $helper_string_0 = " is involved in ";
    my $helper_string_1 = " functions in ";
    my $watch_string = "involved in";
    my @goterms = uniq(sort(@{$goterms_ref}));
    my $helper_string="";
    my $verb_goterm="";
#
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
        foreach my $goterm (@goterms){
        if (@goterms == 1){
         $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         $verb_count++;
         $verb_goterm =  $helper_string . $goterm;      
       } elsif (@goterms == 2){ 

           $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
           if ($count > 0) {
            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

           } else {
             $verb_goterm = $helper_string . $goterm;
           }

       } elsif (@goterms > 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count == 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count < (@goterms-1)) {

            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count >= (@goterms-1)) {
            if ($verb_goterm =~ $helper_string){
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }
         } 
        }
         $count++;
        } # goterms loop
       } # goterms gt 0

 return $verb_goterm; 
}
sub get_verb_function_goterm{
       my $goterms_ref = shift;
       my $ec = shift;

       my $helper_string_before_iea = " is predicted to have ";
       my $helper_string_before_iea_1 = " is predicted to be ";
       my $helper_string_after_iea = "\, based on protein domain information";

       my $helper_string_before_iss = " is predicted to have ";
       my $helper_string_after_iss = "\, based on sequence information";

       my $helper_string_before_exp = " exhibits ";
#      as per R. Kishore's changes
#       my $helper_string_after_exp = "\, based on experimental evidence";
       my $helper_string_after_exp = "";

       my $helper_string_before=""; 
       my $helper_string_after=""; 

       if ($ec =~ /IEA/){
        $helper_string_before=$helper_string_before_iea; 
        $helper_string_after=$helper_string_after_iea;
       } 
       if ($ec =~ /ISS/){
        $helper_string_before=$helper_string_before_iss; 
        $helper_string_after=$helper_string_after_iss;
       } 
       if ($ec =~ /EXP/){
        $helper_string_before=$helper_string_before_exp; 
        $helper_string_after=$helper_string_after_exp;
       } 

       my $verb_goterm="";
       my $watch_string_0 = "activity";
       my $watch_string_1 = "structural constituent";

       my $helper_string_0=$helper_string_before;
       my $helper_string_1=$helper_string_after;

       my @goterms = uniq(sort(@{$goterms_ref}));
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
         my $goterm="";
        foreach $goterm (@goterms){
        chomp($goterm);
        if (@goterms == 1){

         $helper_string_0=$helper_string_before;
         $helper_string_1=$helper_string_after;

           if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = " is a ";
                $helper_string_1 = $helper_string_after;
           }

         $verb_count++;
         $verb_goterm =  $helper_string_0 . $goterm . $helper_string_1;      
       } elsif (@goterms == 2){ 

         $helper_string_0=$helper_string_before;
         $helper_string_1=$helper_string_after;

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

           if ($count > 0) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }

            if ($verb_goterm =~ $helper_string_0){
             $verb_goterm .= " and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;        
             $verb_goterm .= $helper_string_1; 
            }

           } else {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
            }

           if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = $helper_string_after;
           }

               $verb_goterm = $helper_string_0 . $goterm ;
           }

       } elsif (@goterms > 2) {

           if ( $goterm !~/$watch_string_0/){
                $helper_string_0 = "";
                $helper_string_1 = "";
            }


         if ($count == 0) {
            if ( $goterm =~/$watch_string_1/){
                $helper_string_0 = "";
                $helper_string_1 = " and " . $helper_string_before;
           }

           if ( $goterm =~/$watch_string_1/){
                            $verb_goterm = $helper_string_0 . $goterm . $helper_string_1;
           } else {
                             $verb_goterm = $helper_string_0 . $goterm . "\, ";
           }
         } elsif ($count < (@goterms-1)){

             if (($verb_goterm =~ $helper_string_1) and ($goterm !~ $watch_string_1)) {
                $verb_goterm .= $goterm . "\, ";
            } elsif ($verb_goterm =~ $helper_string_0) {
             $verb_goterm .= "\, ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= "\, and "; 
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }
         } elsif ($count >= (@goterms-1)) {
                  $helper_string_0=$helper_string_before;
                  $helper_string_1=$helper_string_after;

            if (($verb_goterm =~ $helper_string_0) and ($goterm !~ $watch_string_1)) {
                 $verb_goterm .= "\, and ";
                 $verb_goterm .= $goterm;
                 $verb_goterm .= $helper_string_1; 
            } elsif ($verb_goterm =~ $helper_string_0){
             $verb_goterm .= "\, and ";
             $verb_goterm .= $goterm;
             $verb_goterm .= $helper_string_1; 
            } elsif (($verb_goterm !~ $helper_string_1) or ($goterm =~ $watch_string_1)) {
             $verb_goterm .= "\, and ";
             $verb_goterm .= $helper_string_0;    
             $verb_goterm .= $goterm;             
             $verb_goterm .= $helper_string_1;  
            }

         } 
        }
         $count++;
        } # goterms loop
       } # goterms gt 0 
 return $verb_goterm;
}
sub get_verb_component_goterm{
       my $goterms_ref = shift;
       my @goterms = uniq(sort(@{$goterms_ref}));
       my $verb_goterm="";
       my $helper_string = "";
       my $watch_string = "integral component of";
       my $watch_string_the = "integral component of the";
       my $helper_string_0 = " is localized to the ";
       my $helper_string_1 = " is an ";
       if (@goterms){
         my $count = 0;
         my $verb_count = 0;
         my $goterm="";
        foreach $goterm (@goterms){
        if (@goterms == 1){
         $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         $verb_count++;
         $verb_goterm =  $helper_string . $goterm;      
       } elsif (@goterms == 2){ 

           $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
           if ($count > 0) {
            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= " and the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

           } else {
             $verb_goterm = $helper_string . $goterm;
           }

       } elsif (@goterms > 2) {
            $helper_string = $helper_string_0;
           if ( $goterm =~/$watch_string/){
                $helper_string = $helper_string_1;
           }
         if ($count == 0) {
             $verb_goterm = $helper_string . $goterm;
         } elsif ($count <  @goterms-1) {

            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= "\, the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and "; 
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }

         } elsif ($count >= (@goterms-1)) {
            if ($verb_goterm =~/$helper_string/){
             $verb_goterm .= " and the ";
             $verb_goterm .= $goterm;
            } else {
             $verb_goterm .= " and ";
             $verb_goterm .= $helper_string;    
             $verb_goterm .= $goterm;               
            }
         } 
        }
         $count++;
        } # goterms loop
       } # goterms gt 0 
 return $verb_goterm;
}
#
sub get_exp_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(sort(@{$go_array_ref}));
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[EXP\]/) or 
  ($go =~/\[IDA\]/) or ($go =~/\[IPI\]/) or 
  ($go =~/\[IMP\]/) or ($go =~/\[IMP\]/) or 
  ($go =~/\[IGI\]/) or ($go =~/\[IEP\]/) ) {
    if (grep {$_ =~/\Q$go\E/} @array) {
        print "$go ignored\; already in array\n";
      } else {
        push(@array, $go);
      }  
  }
 }
  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
sub get_stat_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(sort(@{$go_array_ref}));
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if ( ($go =~/\[ISS\]/) or 
  ($go =~/\[ISO\]/) or ($go =~/\[ISA\]/) or 
  ($go =~/\[ISM\]/) or ($go =~/\[IGC\]/) or 
  ($go =~/\[IBA\]/) or ($go =~/\[IBD\]/) or 
  ($go =~/\[IKR\]/) or ($go =~/\[IRD\]/) or
  ($go =~/\[RCA\]/) ) {
    if (grep {$_ =~/\Q$go\E/} @array) {
        print "$go ignored\; already in array\n";
      } else {
        push(@array, $go);
      }
  }
 }
  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
sub get_iea_go_array{
 my $go_array_ref = shift;
 my @array=();
 my @go_array = uniq(sort(@{$go_array_ref}));
  foreach my $go (@go_array){
        $go =~ s/^\s+//;
        $go =~ s/\s+$//;
   if (($go =~/\[IEA\]/) or 
  ($go =~/\[TAS\]/) or ($go =~/\[NAS\]/) or 
  ($go =~/\[IC\]/) or ($go =~/\[ND\]/)) {
   if (grep {$_ =~/\Q$go\E/} @array) {
       print "$go ignored\; already in array\n";
      } else {
       push(@array, $go);
      } 
  }
 }
  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
sub get_function_string {
 my $functions_array_ref = shift;
 my $ec = shift;
 my $gene_ontology_ref = shift;
 my $function_string = "";
 my @functions = uniq(sort(@{$functions_array_ref}));
 my @GO_functions = ();
 
 if (@functions){
    my $function_count=0;
    foreach my $function (@functions){
          next if ($function =~/GO\:0005488/); # ignore binding
          next if ($function =~/GO\:0005515/); # ignore protein binding
          $function_count++;
          if ($function_count > 0){
          my ($evidence) = $function =~ /\[(.*)\]/;
             if ($evidence) {
                 $function =~ s/\[$evidence\]//g;
             }
             $function =~ s/^\s+//;
             $function =~ s/\s+$//;
     my $go_term = ConciseDescriptions::get_ontology_term($function, $gene_ontology_ref);
     if ($go_term){
      if (lc ($go_term) !~/obsolete/){  
     $go_term .= "[$evidence]";
      if ( $go_term !~ /activity/) { 
          $go_term =~s/binding/binding activity/g;
       }
      }  
     }
      if (grep {$_ =~/\Q$go_term\E/} @GO_functions) {
                print "$go_term ignored\; already in array\n";
      } else {
              push(@GO_functions, $go_term);
      }
    }
    } # foreach functions
    my @uniq_GO_functions = uniq(sort(@GO_functions));
    if (@uniq_GO_functions){
       my $gof_string = "";
       foreach my $gof (@uniq_GO_functions){
        $gof_string .= "$gof\t";
       }
       $function_string = get_verb_function_goterm(\@uniq_GO_functions, $ec);
      }
my @evidences=("EXP","IMP","IGI","IPI","IDA","IEP","IEA","ISS","ISA","ISO","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","ND","IC");
# Remove evidence codes
         foreach my $ec (@evidences){
	    $function_string =~ s/\[$ec\]//g;
         }
      $function_string =~ s/^\s+//;
      $function_string =~ s/\s+$//;
      $function_string =~ s/ \,/\,/g;
# Remove evidence codes
 my $size = @functions;

 } # if gene_functions array

 return $function_string;
}
sub get_gene_array{
 my $file = shift;
 my @array=();
 my @lines = read_file($file);
  foreach my $line (@lines){
   next if ($line =~/\!/);
   my @fields = split(/\t/, $line);
   my $gene_id = $fields[1];
      $gene_id =~ s/^\s+//;
      $gene_id =~ s/\s+$//;
   if ($gene_id){
    push(@array, $gene_id);  
  }
 }
  my @sorted_array = uniq(sort(@array));
  
 return \@sorted_array;
}
