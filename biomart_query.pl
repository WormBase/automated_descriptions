#!/usr/bin/env perl
use strict;
use warnings;

my @args = ("./biomartWebExample.pl", "./BioMartQuery.xml");
my $status = system(@args);
exit 0;
