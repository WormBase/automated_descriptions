#!/usr/bin/env perl
use strict;
use warnings;
use ConciseDescriptions;
use List::MoreUtils qw(uniq);
use File::Slurp;
use LWP::Simple;

my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
my $RELEASE = ConciseDescriptions::get_release();
my $html = ConciseDescriptions::get_html_dir();
my $home = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/";
my $gene_list_dir = $home . "gene_lists/";
my $db_gene_list  = $gene_list_dir . "wb_gene_list.txt";

my $synList ="http://textpresso-dev.caltech.edu/celegans/tdb/celegans/synonymList/c_elegansGeneSynList";
my $synonym_hash_ref=get_synonym_hash($synList);
my %synonym_hash = %{$synonym_hash_ref};

my $deadList = "http://brahma.textpresso.org/concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/gene_lists/sort.dead_genes.txt";
my $dead_hash_ref = get_dead_hash($deadList);
my %dead_hash = %{$dead_hash_ref};
my $output_dir = $html . "concise_descriptions/textpresso/";
my $output = $output_dir . "textpresso_papers_results_genes.all.txt";
if (-e $output){
  my @args = ("rm", "-f", $output);
  system(@args) == 0 or die("could not delete file $output\n");
}
my $header="  paper id\t\tgenes\n";
write_file($output, {append => 1 }, $header);
 
 my %hash=();
 my %occurence_hash=();
 my %unique_gene_list_hash=();
 my @genes=();
 my @all_papers=();
 my @sorted_genes=(); 

 if (-e $db_gene_list){

 my @lines = read_file($db_gene_list);
  foreach my $line (@lines){
  chomp($line);
  my @fields = split(/\t/, $line);
  my $gene_name = $fields[1];
  my $gene_id = $fields[0];
  push(@genes, $gene_id);
  $hash{$gene_id} = $gene_name;
 }
  @sorted_genes = sort(@genes);
foreach my $gene_id (@sorted_genes){
 my @papers=();
 next if ($dead_hash{$gene_id});
 my $gene_list = $hash{$gene_id};
 if ($synonym_hash{$gene_id}){
  $gene_list .= "\,$synonym_hash{$gene_id}";
 }
 my @genes = split(/\,/, $gene_list);
 my @unique_genes = uniq(@genes);
 my $unique_gene_list = join(',', @unique_genes);
 foreach my $gene (@unique_genes){
  if (($gene=~/\S+\-\d+/) or ($gene=~/\S+\.\S+/)){
  my $first_letter = substr($gene, 0, 1);
  my $second_letter = substr($gene, 1, 1);
# first case: lc first . lc second . lc rest of gene
$first_letter = lc($first_letter);
$second_letter = lc($second_letter);
$gene = lc($gene);
my ($count_hash_ref, $wbpapers_ref) = get_papers($first_letter, $second_letter, $gene);
my @wbpapers = @{$wbpapers_ref};
my %count_hash = %{$count_hash_ref};
foreach my $wb (@wbpapers){
 my $newid=$wb . "AND" . $gene_id;
 my $count = $count_hash{$wb};
 if ($occurence_hash{$newid}){
  $occurence_hash{$newid}+=$count;
 } else {
  $occurence_hash{$newid}=$count;
 }
 push(@papers, $wb);
 }
# second case: uc first . uc second . uc rest of gene
$first_letter = uc($first_letter);
$second_letter = uc($second_letter);
$gene = uc(lc $gene);
  ($count_hash_ref, $wbpapers_ref) = get_papers($first_letter, $second_letter, $gene);
  @wbpapers = @{$wbpapers_ref};
  %count_hash = %{$count_hash_ref};
foreach my $wb (@wbpapers){
 my $newid=$wb . "AND" . $gene_id;
 my $count = $count_hash{$wb};
 if ($occurence_hash{$newid}){
  $occurence_hash{$newid}+=$count;
 } else {
  $occurence_hash{$newid}=$count;
 }
 push(@papers, $wb);
 }
# first case: uc first . lc second . lc rest of gene
$first_letter = uc($first_letter);
$second_letter = lc($second_letter);
$gene = ucfirst(lc $gene);
  ($count_hash_ref, $wbpapers_ref) = get_papers($first_letter, $second_letter, $gene);
  @wbpapers = @{$wbpapers_ref};
  %count_hash = %{$count_hash_ref};
  foreach my $wb (@wbpapers){
 my $count = $count_hash{$wb};
 my $newid=$wb . "AND" . $gene_id;
 if ($occurence_hash{$newid}){
  $occurence_hash{$newid}+=$count;
 } else {
  $occurence_hash{$newid}=$count;
 }
   push(@papers, $wb);
  }
#print "There are $count publications\.\n";
$gene = lc($gene);
} # if gene matches xxx-# or xxx.x
} # for each gene (each synonym of gene id)
my @unique_papers = uniq(@papers);
foreach my $p (@unique_papers){
 push(@all_papers, $p);
}
$unique_gene_list =~s/\#//g;
$unique_gene_list =~s/\s+//g;
$unique_gene_list_hash{$gene_id}=$unique_gene_list;
#my $output_string ="$gene_id\t$unique_gene_list\t$count\n";
#write_file($output, {append => 1 }, $output_string);

} # foreach gene id
} # if db gene list exists

my @unique_all_papers = sort(uniq(@all_papers));
foreach my $p (@unique_all_papers){
 my $string ="$p\t";
 write_file($output, {append => 1 }, $string);
 foreach my $gene_id (@sorted_genes){
 next if ($dead_hash{$gene_id});
 my $newid = $p . "AND" . $gene_id;
 if ($occurence_hash{$newid}){
 my $count = $occurence_hash{$newid};
 if ($count > 0){
  my $name = $hash{$gene_id};
  my $output_string ="$name\($count\) ";
  write_file($output, {append => 1 }, $output_string);
  }
 }
 }
 write_file($output, {append => 1 }, "\n");
}

exit 0;

sub valid_url {
my $url = shift;
my $status = 1;
if (head($url)) {
  $status = 0;
} else {
   $status = 1;
}
return $status;
}

sub get_papers {
 my $first_letter = shift;
 my $second_letter = shift;
 my $gene = shift;
 my @papers=();
 my @sections = qw(results);
 my %count_hash=();
 foreach my $section (@sections){
  my $url_base = "http://textpresso-dev.caltech.edu/celegans/tdb/celegans/ind/$section/keyword/";
  my $webpage = $url_base . $first_letter . "/" . $second_letter . "/" . $gene;
  my $status = valid_url($webpage);
   if ($status ==0){
    my $content = ConciseDescriptions::getwebpage($webpage);
     if ($content !~/404 Not Found/){
       my @lines = split(/WBPaper/, $content);
       foreach my $line (@lines){
        chomp($line);
        $line =~s/^\s+//;
        $line =~s/\s+$//;
        next if ($line =~/^\s*$/);
        my ($paper, $place) = split(/\#/, $line);
        $paper =~s/^\s+//;
        $paper =~s/\s+$//;
        $place =~s/^\s+//;
        $place =~s/\s+$//;
        my @places = split(/ /,$place);
        my $count = @places;
#        $line =~s/\#.*//g;
#        $line =~s/\.*//g;
#        $line =~s/\.//g;
        if ($paper =~/^(.*)\./){
            $paper = $1;
        }
        if ($paper !~/^WBPaper/){
         $paper = "WBPaper" . $paper;
        }
         $count_hash{$paper} = $count;
         push(@papers, $paper);
       }
      }
     }
    }

 return (\%count_hash, \@papers);
}
sub get_synonym_hash{
my $webpage = shift;
my %hash = ();
my $status = valid_url($webpage);
   if ($status ==0){
    my $content = ConciseDescriptions::getwebpage($webpage);
     if ($content !~/404 Not Found/){
        my @lines = split(/\n/, $content);  
        foreach my $line (@lines){
         chomp($line);
         # remove leading whitespace
         $line =~ s/^\s+//;
         # remove trailing whitespace
         $line =~ s/\s+$//;
         my $id = substr($line, 0, index($line, ","));
         # remove leading whitespace
         $id =~ s/^\s+//;
         # remove trailing whitespace
         $id =~ s/\s+$//;
         my ($rest) = $line; 
             $rest =~s/$id\,//;
         if ($rest=~/his/){
           $rest =~s/his /his\-/g;
          }
         # remove leading whitespace
         $rest =~ s/^\s+//;
         # remove trailing whitespace
         $rest =~ s/\s+$//;
         $hash{$id} = $rest;
      }  
     }
    }
return \%hash;
}
sub get_dead_hash{
my $webpage = shift;
my %hash = ();
my $status = valid_url($webpage);
   if ($status ==0){
    my $content = ConciseDescriptions::getwebpage($webpage);
     if ($content !~/404 Not Found/){
        my @lines = split(/\n/, $content);  
        foreach my $line (@lines){
         chomp($line);
         my $id = $line;
         $hash{$id} = $id;
      }  
     }
    }
return \%hash;
}
