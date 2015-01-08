#!/usr/bin/env perl
#
# convert obo file to .ace format according to wiki instructions
# http://wiki.wormbase.org/index.php/Updating_go.ace_file
# 
use warnings;
use strict;
use ConciseDescriptions;
use LWP::Simple;
use File::Slurp;

# get from URL
# html path is defined by reading the an input file (e.g., docroot in RH, /var/www in Ubuntu, etc.)
my $html = ConciseDescriptions::get_html_dir();
my $PRODUCTION_RELEASE = ConciseDescriptions::get_production_release();
#
my $path = $html . "concise_descriptions/release/$PRODUCTION_RELEASE/c_elegans/gene_lists/";
my $url = 'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo';
my $local_file = $path . "gene_ontology_ext.obo";

my $obo_file = get $url;
my $outfile = $path . 'go_terms.ace';
open (OUT, ">$outfile") or die "Cannot create $outfile : $!";

my (@terms) = split/\[Term\]/, $obo_file;
my $header = shift @terms;
my $version = 'unknown';
if ($header =~ m/remark: cvs version: \$(.*?)\$/) { $version = $1; }
# print "$version\n";

my %hash;

foreach my $entry (@terms) {
  next unless ($entry =~ m/id: (GO:\d+)/);
  my $id = $1;
#   print "ID $id\n";
  $hash{$id}{main} .= qq(GO_term: "$id"\n); 
  if ($entry =~ m/is_obsolete: true/) { $hash{$id}{main} .= qq(Status\t"Obsolete"\n); }
    else { $hash{$id}{main} .= qq(Status\t"Valid"\n); }
  if ($entry =~ m/def: "(.*?)"/) { $hash{$id}{main} .= qq(Definition\t"$1"\n); }
  if ($entry =~ m/name: (.*)/) { $hash{$id}{main} .= qq(Term\t"$1"\n); }
  if ($entry =~ m/synonym: "(.*?)" ([A-Z]+)/) { 
    my (@pairs) = $entry =~ m/synonym: "(.*?)" ([A-Z]+)/g;
    while (scalar @pairs > 0) {
      my $syn = shift @pairs; my $type = shift @pairs;
      if ( ($type eq 'BROAD') || ($type eq 'EXACT') || ($type eq 'NARROW') || ($type eq 'RELATED') ) { 
          $type = ucfirst(lc($type));
          $hash{$id}{main} .= qq($type\t"$syn"\n); }
        else { print qq(ERR $type not a valid type in synonym line in $entry\n); } } }
  if ($entry =~ m/namespace: (.*)/) { my $namespace = $1; $namespace = ucfirst($namespace); $hash{$id}{main} .= qq($namespace\n); }
  $hash{$id}{main} .= qq(Version\t"Gene Ontology $version"\n);
  
  if ($entry =~ m/is_a: (GO:\d+)/) { 
    my (@list) = $entry =~ m/is_a: (GO:\d+)/g;
    foreach my $isa (@list) { $hash{$id}{isa}{$isa}++; $hash{$isa}{isad}{$id}++; $hash{$id}{anc}{$isa}++; $hash{$isa}{desc}{$id}++; } }
  if ($entry =~ m/relationship: part_of (GO:\d+)/) {
    my (@list) = $entry =~ m/relationship: part_of (GO:\d+)/g;
    foreach my $pof (@list) { $hash{$id}{pof}{$pof}++; $hash{$pof}{pofd}{$id}++; $hash{$id}{anc}{$pof}++; $hash{$pof}{desc}{$id}++; } }
} # foreach my $entry (@terms)

my %recurse;		# filter through this hash to avoid multiple elements getting the same ancestors/descendants from different terms that share the same anc/desc
foreach my $id (sort keys %hash) {
  if ($hash{$id}{main}) { 
    print OUT qq($hash{$id}{main});
    if ($hash{$id}{isad}) {
      foreach my $isad (sort keys %{ $hash{$id}{isad} }) { print OUT qq(Instance\t"$isad"\n); } }
    if ($hash{$id}{isa}) {
      foreach my $isa (sort keys %{ $hash{$id}{isa} }) { print OUT qq(Instance_of\t"$isa"\n); } }
    if ($hash{$id}{pofd}) {
      foreach my $pofd (sort keys %{ $hash{$id}{pofd} }) { print OUT qq(Component\t"$pofd"\n); } }
    if ($hash{$id}{pof}) {
      foreach my $pof (sort keys %{ $hash{$id}{pof} }) { print OUT qq(Component_of\t"$pof"\n); } }
    if ($hash{$id}{anc}) {
      %recurse = ();			# clear out recurse filter hash
      &recurse($id, 'anc');		# recurse through all depth in tree and add each term to the recurse filter hash
      foreach my $newid (sort keys %recurse) { print OUT qq(Ancestor\t"$newid"\n); }
    }
    if ($hash{$id}{desc}) {
      %recurse = ();
      &recurse($id, 'desc');
      foreach my $newid (sort keys %recurse) { print OUT qq(Descendent\t"$newid"\n); }
    }
    print OUT qq(\n);
  } # if ($hash{$id}{main}) 
} # foreach my $id (sort keys %hash)

close (OUT) or die "Cannot close $outfile : $!";

sub recurse {				# go through each parent/child depeding on direction and recurse through this function to get full depth in that direction.  add all terms to %recurse hash to filter multiple values that are the same 
  my ($id, $direction) = @_;
  if ($hash{$id}{$direction}) {
    foreach my $newid (sort keys %{ $hash{$id}{$direction} }) {
      $recurse{$newid}++;
      &recurse($newid, $direction); } }
} # sub recurse
