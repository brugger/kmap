#!/usr/bin/perl 
# 
# Memory hungry kmer counter implementation to validate the C program. Wayyyy to goooo
# 
# 
# Kim Brugger (20 Jan 2012), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $fasta_file = shift;

my $fasta_hash = readfasta( $fasta_file );

my $seq = $$fasta_hash{'gid:3052'};

#print "$seq\n";

my $kmer_size = 36;

my %k_count;

for(my $i=0;$i < length($seq)-$kmer_size+1;$i++) {

  $k_count{ substr($seq, $i, $kmer_size)}++; 
  
}

foreach my $k (sort keys %k_count ) {
  print "$k\t$k_count{ $k }\n";
}

#
# Read the fasta files and puts entries into a nice array
#
sub readfasta {
  my ($file) = @_;  

  my %res;

  my $sequence;
  my $header;

  open (my $f, "$file" ) || die "Could not open $file: $1\n";
  while (<$f>) {
    chomp;
    if (/^\>/) {
      if ($header) { # we have a name and a seq
	$res{ $header } =  uc($sequence);
      }
      $header = $_;
      $header =~ s/^\>//;
      $sequence = "";
    }
    else {$sequence .= $_;}
  }

  $res{ $header } = $sequence ;

  close( $f );

  return \%res;
}
