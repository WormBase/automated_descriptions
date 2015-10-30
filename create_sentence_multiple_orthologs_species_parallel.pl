#!/usr/bin/env perl
use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use List::Util qw/first/;
use File::Slurp;
use Text::CSV;
use warnings;
use strict;

my $species = $ARGV[0];
my $project = $ARGV[1];
my $name = $ARGV[2];
my $prefix= $ARGV[3];

chomp($species);
chomp($project);
chomp($name);
chomp($prefix);
       
$project =~ s/^\s+//;
$project =~ s/\s+$//;
$species =~ s/^\s+//;
$species =~ s/\s+$//;
$name    =~ s/^\s+//;
$name    =~ s/\s+$//;
$prefix  =~ s/^\s+//;
$prefix  =~ s/\s+$//;

my $RELEASE = ConciseDescriptions::get_release();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/";
# The list of genes with no concise descriptions is in $infile
my $elegans_gene_list_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $db_gene_list  = $elegans_gene_list_dir . "wb_gene_list.txt";

my $curated_gene_list = $elegans_gene_list_dir . "sort.curated_genes.txt";
my @curated_genes = read_file($curated_gene_list);

my $orthology = $home . "release/$PRODUCTION_RELEASE/$species/orthology/";
my $orthology_elegans_dir = $home . "release/$PRODUCTION_RELEASE/c_elegans/orthology/";

my $popularity_gene_id_name_file = $elegans_gene_list_dir . "textpresso_gene_popularity.text";
my $popularity_gene_id_name_hash_ref = get_gene_popularity_hash($popularity_gene_id_name_file);

my $wormbase_gene_id_name_file = $elegans_gene_list_dir . "wormbase_gene_id_name.list";
my $wormbase_gene_id_name_hash_ref = get_wbgene_name_hash($wormbase_gene_id_name_file);

my $gene_class_file = $elegans_gene_list_dir . "acedb_gene_class.txt";
my $gene_class_hash_ref = ConciseDescriptions::get_gene_class_hash($gene_class_file);

my $sapiens_orthology = $orthology . "input_files/orthologs.$name.H_Sapiens.txt";
my $orthology_elegans = $orthology_elegans_dir . "input_files/c_elegans.orthologs.txt";
my $orthology_elegans_species = $orthology . "input_files/orthologs.$species.elegans.txt";
my $orthology_species = $orthology . "input_files/$species.orthologs.txt";
my $output_file = $orthology . "output_files/sentences_for_orthology.txt";
#
if (-e $orthology){
   print "$orthology exists\n"; 
} else {
   mkdir $orthology or die "could not create $orthology";
}
my $output_path = $orthology;
my $individual_path = $output_path . "output_files/individual_gene_sentences/";
if (-e $individual_path){
   print "$individual_path exists\n"; 
} else {
   mkdir $individual_path or die "could not create $individual_path";
}
# if the output file exists delete it.
if (-e $output_file){
   my @args = ("/bin/rm", "-f", $output_file);
   system(@args) == 0 or die("could not delete file $output_file\n");
}
# if the output files for individual gene sentences exist delete them.
   my @individual_files = glob("$individual_path/WBGene*");
      foreach my $individual (@individual_files){
      my @args = ("/bin/rm", "$individual");
      system(@args) == 0 or die("could not delete $individual\n");
   }

my ($gene_elegans_id_hash_ref, $gene_elegans_hash_ref) = get_id_name_hash($orthology_elegans_species, $db_gene_list);
my %gene_elegans_hash= %{$gene_elegans_hash_ref};
my %gene_elegans_id_hash= %{$gene_elegans_id_hash_ref};

my $gene_name_hash_ref = get_gene_name_hash($orthology_species, $db_gene_list);
my %gene_name_hash = %$gene_name_hash_ref;

my $gene_array_ref = get_gene_array($orthology_species);
my @gene_array = @{$gene_array_ref};

my $and = "and";
my $AND = "AND";
my $comma = "\,";
my $is_a = " is an ortholog of human ";
my $is_a_ce = " is an ortholog of C\. elegans ";
my $human = " and an ortholog of human ";

my $doublespace = "  ";
my $space = " ";
my $spacesemicolon = " \;";
my $semicolon = "\;";

# Define the various model organisms used
my @mods=("S\. cerevisiae","yeast","S\. pombe","mouse","Drosophila","Chlamydomonas","zebrafish","Rat","rat","chicken","Strongylocentrotus purpuratus","Arabidopsis","Dictyostelium","bacterial","E\. coli","H\. influenzae","Dictyostelium");
 my $elegans = "C\. elegans";

my $dead_gene_list = $elegans_gene_list_dir . "sort.dead_genes.txt";
my @dead_genes = read_file($dead_gene_list);

my @uncurated_genes_array = ();
foreach my $test (@gene_array){
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
my @uncurated_live_genes_array = ();
foreach my $test (@sorted_uncurated_genes_array){
 my $keep =0;
 foreach my $dead (@dead_genes){
   if ($dead =~/$test/){
    $keep = 1;
   }
 }
 if ($keep ==0){
   push(@uncurated_live_genes_array, $test);
 }
}
my @sorted_uncurated_live_genes_array = uniq(sort(@uncurated_live_genes_array)); 
 
foreach my $gene_id (@sorted_uncurated_live_genes_array){
  my $sentence  = "";
  chomp($gene_id);
       $gene_id =~ s/^\s+//;
       $gene_id =~ s/\s+$//;
  my $out = $individual_path . $gene_id;
  my $gene_name = $gene_name_hash{$gene_id};
#
# Fill array, @elegans_array, with list of WBGeneID orthologs.
#
  my $elegans_list = "";
  my @elegans_array = ();
  if ($gene_elegans_id_hash{$gene_id}){
      $elegans_list = $gene_elegans_id_hash{$gene_id};
   if ($elegans_list =~/$AND/){
      @elegans_array = split(/$AND/, $elegans_list);
    } else {
     $elegans_array[0] = $elegans_list;
   }
  }

  my $elegans_phrase = "";
  my $count_elegans = 0;
  my $elegans_size = @elegans_array;
#  print "$gene_name\t$elegans_size\n";
  if ($elegans_size < 6){

  foreach my $elegans (@elegans_array){
       $elegans =~ s/^\s+//;
       $elegans =~ s/\s+$//;
    my $gene_elegans = $gene_elegans_hash{$elegans};
    if ($gene_elegans){
    if ($count_elegans == 0){
     $elegans_phrase = $gene_elegans;  
    } elsif (($count_elegans > 0) and ($count_elegans < ($elegans_size -1))){
     $elegans_phrase .= "\, " . $gene_elegans;  
    } elsif (($count_elegans > 0) and ($count_elegans == ($elegans_size -1))){
     $elegans_phrase .= " and " . $gene_elegans;  
    }
   }
    $count_elegans++;
  } # end foreach my $elegans

  } else {

my $list = join(",",@elegans_array);
#print "$gene_name\t$list\n";
$sentence = $gene_name;
$sentence = reduce_ortholog_sentence($sentence, $list,  $popularity_gene_id_name_hash_ref, $wormbase_gene_id_name_hash_ref,  $gene_elegans_hash_ref, $gene_elegans_id_hash_ref, $gene_name_hash_ref, $gene_class_hash_ref);
$sentence .= "\;\n";
  }

if (($elegans_phrase) and ($elegans_size < 6)) {
      $sentence  = $gene_name . $is_a_ce . $elegans_phrase;
      $sentence .= "\;\n";
}

if ($sentence) {
      write_file($output_file, {append => 1 }, $gene_id);      
      write_file($output_file, {append => 1 }, "\n");
      write_file($output_file, {append => 1 }, $sentence);
      write_file($output_file, {append => 1 }, "\n\n\n");
      write_file($out, $sentence);
} 
 
} # end foreach my $gene_id (@gene_array)
exit 0;

sub reduce_ortholog_sentence{

my $sentence=shift;
my $ortholog_list = shift;
my $popularity_gene_id_name_hash_ref = shift;
my $wormbase_gene_id_name_hash_ref = shift;
my $gene_elegans_hash_ref = shift;
my $gene_elegans_id_hash_ref = shift;
my $gene_name_hash_ref = shift;
my $gene_class_hash_ref = shift;

my %popularity_gene_id_name_hash = %{$popularity_gene_id_name_hash_ref};
my %wbgene_name_hash = %{$wormbase_gene_id_name_hash_ref};
my %gene_elegans_hash= %{$gene_elegans_hash_ref};
my %gene_elegans_id_hash= %{$gene_elegans_id_hash_ref};
my %gene_name_hash = %$gene_name_hash_ref;
my %gene_class_hash = %{$gene_class_hash_ref};
my @orthologs = split(/\,/,$ortholog_list);

my $is_an_ortholog = " is an ortholog of C\. elegans ";
my $is_an_ortholog2 = " is an ortholog of the members of the C\. elegans ";
my $and_members = " and members of the C\. elegans ";

my %family_hash=();
my %family_elements_hash=();
my %family_popularity_hash=();
my %family_popularity_number_hash=();
my %family_popularity_order_hash=();
my @families=();

if (@orthologs >= 5){
  foreach my $ortholog (@orthologs){
   my $family;
   if ($gene_class_hash{$ortholog}){
       $family = $gene_class_hash{$ortholog};
   } else {
       $family = "cosmids";
   }
   push(@families, $family);
   if ($family_hash{$family}){
     $family_elements_hash{$family} += 1; 
     $family_hash{$family} .= "\,$ortholog";
 } else {
     $family_elements_hash{$family} = 1; 
     $family_hash{$family} = $ortholog;
   }   
 }
}

@families = uniq(@families);

#if ($family_hash{"cosmids"}){
# my $c = $family_hash{"cosmids"};
#}
my %sorted_gene_family=();
foreach my $f (@families){
  my @genes = split(/\,/, $family_hash{$f});
  my %popularity=();
  foreach my $id (@genes){
   my $pop = 0;
   if ($popularity_gene_id_name_hash{$id}){
      $popularity{$id} = $popularity_gene_id_name_hash{$id};
   } else {
      $popularity{$id} = 0;
   }
  }
  my @sorted_genes;
  my $count = $family_elements_hash{$f};
  if ($count > 1) {
   my $sorted_genes_ref = by_popularity(\%popularity);
   @sorted_genes = @{$sorted_genes_ref};
  } elsif ($count ==1) {
   @sorted_genes = $family_hash{$f};
  }
  $sorted_gene_family{$f}=join(',', @sorted_genes);
}

my @singlets=();
my @multiples=();
my @cosmids=();

foreach my $f (@families){
 my $count;
 next if ($f=~/cosmids/);
   $count = $family_elements_hash{$f};
   if ($count == 1){
    push (@singlets, $f);
   } elsif ($count > 1) {
    push (@multiples, $f);
  }
 }

@singlets = uniq(@singlets);
@multiples = uniq(@multiples);

my $s_size=0;
foreach my $s (@singlets){
  next if ($s =~/cosmids/);
  $s_size++;
}

if ($s_size > 0){
 $sentence .= $is_an_ortholog;
}

my $kount = 0;
my $size = @singlets;

foreach my $s (@singlets){
  next if ($s =~/cosmids/);
  my $gene_id = $family_hash{$s};
  my $gene = $wbgene_name_hash{$gene_id};
  if ($kount == 0){
     $sentence .= $gene;  
    } elsif (($kount > 0) and ($kount < ($s_size-1))){
     $sentence .= "\, " . $gene;  
    } elsif (($kount > 0) and ($kount >= ($s_size-1)) and ($s_size > 1) ){
     $sentence .= " and " . $gene;  
    }
 $kount++;
}

$kount = 0;
$size = @multiples;

$kount = 0;
my $m_count = 0;
my $limit = 0;
my $m_size=0;
foreach my $m (@multiples){
 next if ($m =~/cosmids/);
 $m_size++;
}

if (($m_size > 0) and ($s_size > 0)){
 $sentence .= " and members of the C. elegans ";
} elsif ($m_size > 0){
 $sentence .= $is_an_ortholog2;
}
foreach my $m (@multiples){
  next if ($m =~/cosmids/);
  if ($kount == 0){
     $sentence .= $m;  
    } elsif (($kount > 0) and ($kount < ($m_size-1))){
     $sentence .= "\, " . $m;  
    } elsif (($kount > 0) and ($kount >= ($m_size-1)) and ($m_size > 1)){
     $sentence .= " and " . $m;  
    }
 $kount++;
}

my $gene_class_statement;

 if ($m_size ==1) {
  $gene_class_statement = " gene class including ";
 } elsif ($m_size > 1){
  $gene_class_statement = " gene classes including ";  
 }

if ($m_size > 0){
 $sentence .= $gene_class_statement;
}
if ($m_size == 1){
  $limit = 3;
 } elsif ($m_size > 1){
  $limit = 1;
 }
foreach my $m (@multiples){
  next if ($m =~ /cosmids/);
  $m_count++;
  my $gene_array_string = $sorted_gene_family{$m};
  my @sorted_genes = split(/\,/, $gene_array_string);
  my $size = @sorted_genes;
  my $newlimit = 0;
  if ($size < $limit){
    $newlimit = $size;
  }
  my $g_count = 0;
  $kount = 0;
  foreach my $g (@sorted_genes){
     chomp($g);
     $g_count++;
     next if ($g_count > $limit);
   my $gene = $wbgene_name_hash{$g};
   if ($newlimit == 0){
   if ($kount == 0){
     $sentence .= $gene;  
    } elsif (($kount > 0) and ($kount < ($limit-1))){
     $sentence .= "\, " . $gene;  
    } elsif (($kount > 0) and ($kount == ($limit-1)) and ($limit > 1) ){
     $sentence .= " and " . $gene;  
    }
   } else {
    if ($kount == 0){
     $sentence .= $gene;  
    } elsif (($kount > 0) and ($kount < ($newlimit-1))){
     $sentence .= "\, " . $gene;  
    } elsif (($kount > 0) and ($kount == ($newlimit-1)) and ($newlimit > 1) ){
     $sentence .= " and " . $gene;  
    }
   }
 $kount++;

 }

    if (($m_size > 1) and ($m_count < $m_size)) {
      $sentence .= " and ";
    }
}
my $COSMIDS = "cosmids";
my $cosmid_family = $family_hash{$COSMIDS};
if ($cosmid_family){
my $cosmid_array_string = $sorted_gene_family{$COSMIDS};
my @cosmids = split(/\,/, $cosmid_array_string);

if ((@cosmids) and ($m_size > 0)){
 $sentence .= "\, and ";
} elsif ((@cosmids) and ($s_size > 0)){
 $sentence =~s/ and /\, /g;
 $sentence .= "\, ";
} else {
 $sentence .= $is_an_ortholog;
}

$kount=0;
my $size = @cosmids;
foreach my $c (@cosmids){
      $c =~s/^\s+//;
      $c =~s/\s+$//;
   my $gene = $wbgene_name_hash{$c};
  if ($kount == 0){
     $sentence .= $gene;  
    } elsif (($kount > 0) and ($kount < ($size-1))){
     $sentence .= "\, " . $gene;  
    } elsif (($kount > 0) and ($kount == ($size-1)) and ($size > 1)){
     $sentence .= " and " . $gene;  
    }
 $kount++;
}
}

 return $sentence;
}

sub get_gene_popularity_hash{
my $file = shift;
my %hash=();

my @lines = read_file($file);
foreach my $line (@lines){
next if ($line =~/id/);
chomp($line);
my ($wbgene_id, $name, $number) = split(/\t/, $line);
        $wbgene_id =~s/^\s+//;
        $wbgene_id =~s/\s+$//;
        $number =~s/^\s+//;
        $number =~s/\s+$//;
 if (($number) and ($wbgene_id)){
   $hash{$wbgene_id} = $number;
 }
}

return \%hash;
}

sub get_wbgene_name_hash{
my $file = shift;
my %hash=();
my @lines = read_file($file);
foreach my $line (@lines){
chomp($line);
$line =~s/\"//g;
my ($wbgene_id, $name) = split(/\t/, $line);
 if (($name) and ($wbgene_id)){
      $name    =~s/^\s+//;
      $name    =~s/\s+$//;
      $wbgene_id =~s/^\s+//;
      $wbgene_id =~s/\s+$//;
   $hash{$wbgene_id} = $name;
 }
}
return \%hash;
}

#sub get_gene_class_hash{
#my $file = shift;
#my %hash=();
#
#my @lines = read_file($file);
#foreach my $line (@lines){
#chomp($line);
#$line =~s/\"//g;
#my ($family, $wbgene_id) = split(/\t/, $line);
# if (($family) and ($wbgene_id)){
#      $family    =~s/^\s+//;
#      $family    =~s/\s+$//;
#      $wbgene_id =~s/^\s+//;
#      $wbgene_id =~s/\s+$//;
#   $hash{$wbgene_id} = $family;
# }
#}
#
#return \%hash;
#}

sub by_popularity{
my $hash_ref=shift;
my %hash = %{$hash_ref};
my @array=();

   foreach my $key (reverse sort keys %hash) {
#    print "$key \: " . $popularity{$key} . "\n";
     push(@array, $key);
   }

 return \@array;
}

sub get_id_name_hash{
 my $orthology=shift;
 my $db_gene_list  = shift;

 my @lines = read_file($orthology);
 my $old_gene="";
 my %gene_name_hash=();
 my %gene_id_hash=();
 my @multiple_source_array=();
 my @single_source_array=();
 my $oldgene="";
 my $newgene; 
 foreach my $line (@lines){ 
    chomp($line);   
    if ($line =~/WBGene/){
     my ($wb_gene_id, $gene_name, $species_name, $id, $name, $source) = split(/\t/, $line); 
        $wb_gene_id =~ s/^\s+//;
        $wb_gene_id =~ s/\s+$//;
        $newgene = $wb_gene_id;
     if ($newgene ne $oldgene){
         my $single_size = @single_source_array;
         my $multiple_size = @multiple_source_array;

         if ($multiple_size > 0){
             foreach my $a_line (@multiple_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
                       $gene_name_hash{$a_id}  = $a_name;
 
               if ($gene_id_hash{$a_wb_gene_id}){
                   $gene_id_hash{$a_wb_gene_id} .= "AND" .  $a_id;
               } else{
                   $gene_id_hash{$a_wb_gene_id} = $a_id;
               }
            }
          } elsif ($single_size > 0) {
             foreach my $a_line (@single_source_array){
               my ($a_wb_gene_id, $a_gene_species_name, $a_id, $a_name, $a_source) = split(/\t/, $a_line);
               my ($a_gene_name, $a_species_name) = split(/ /, $a_gene_species_name);
#                   $gene_name_hash{$a_id}  = $a_name;
#                   if not ($gene_name_hash{$a_id}){
                       $gene_name_hash{$a_id}  = $a_name;
#                    } 

               if ($gene_id_hash{$a_wb_gene_id}){
                   $gene_id_hash{$a_wb_gene_id} .= "AND" .  $a_id;
               } else{
                   $gene_id_hash{$a_wb_gene_id} = $a_id;
               }
          }
         }
         @multiple_source_array=();
         @single_source_array=();
         $oldgene=$newgene;
     }
         
         if ($line =~/\;/){
          my $multiple = @multiple_source_array;
           push(@multiple_source_array, $line);
       } elsif ($line) {
          my $single = @single_source_array;
           push(@single_source_array, $line);
       }
  }
 }

 my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
     }
 }

 return \%gene_id_hash, \%gene_name_hash;
}
sub get_gene_array{
 my $orthology = shift;
 my @gene_array=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  chomp($line);
  if ($line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        push(@gene_array, $WBGene);
     }
 }
      my @sorted_gene_array = uniq(sort(@gene_array));
      return \@sorted_gene_array;
}
sub get_gene_name_hash{
 my $orthology = shift;
 my $db_gene_list = shift;
 my %gene_name_hash=();
 my @lines = read_file($orthology);
 foreach my $line (@lines){
  next if ($line =~ /\#/);
  next if ($line =~ /\=/);
  chomp($line);
  if ($line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
     }
 }
 my @db_lines = read_file($db_gene_list);
 foreach my $db_line (@db_lines){
  chomp($db_line);
  if ($db_line =~/WBGene/){
     my ($WBGene, $Gene) = split(/\t/, $db_line);         
        $WBGene =~ s/^\s+//;
        $WBGene =~ s/\s+$//;
        $Gene =~ s/^\s+//;
        $Gene =~ s/\s+$//;
        $gene_name_hash{$WBGene} = $Gene;
     }
 }
      return \%gene_name_hash;
}
