#!/usr/bin/env perl
use strict;

if(@ARGV != 11)
{
        print "Example: perl split_bc.pl barcode.list read_1.fq.gz read_2.fq.gz read_len split_read skip bc_len_select bc_len_redundant bc_pos skip_as_bc bc_trim_gdna \n";
        exit(0);
}
my $read_len = $ARGV[3]; # Get read length
my $skip = $ARGV[5];
my $bc_len_select = $ARGV[6];
my $bc_len_redundant = $ARGV[7];
my $bc_pos = $ARGV[8];
my $skip_as_bc = $ARGV[9];
my $bc_trim_gdna = $ARGV[10];
# print "$bc_pos";
# my ($bc_len_redundant, $bc_len_select, $skip) = (20, 15, 20); # 120=20+80+20
# n1= 15bp BC (random, no A) + 13bp repeat CGCTTTGTTCGTG​ + remaining; n2= 18bp BC; get barcode and determine lengths, V350085271 PE100+140
# SE120=skip20+80gDNA+15bc+5
# ($bc_len_redundant, $bc_len_select, $skip) = (20, 15, 20)= (15bc+redundant, 15bc as bx, 20trim)
# randome BC new design, bs as seq not number, skip 5' 20bp, add 3' 20bp as bc on r2, then alignment for new_bc

my %barcode_hash; # hash map of barcodes

# open IN,"$ARGV[0]" or die "can't open barcode.list"; # open barcode list
my $n = 0; # initialize counter


open IN1,"gzip -dc $ARGV[1] |" or die "cannot open read1"; # open read1
open IN2,"gzip -dc $ARGV[2] |" or die "cannot open read2"; # open read2
open OUT1, "| gzip > $ARGV[4].1.fq.gz" or die "Can't write file";
open OUT2, "| gzip > $ARGV[4].2.fq.gz" or die "Can't write file";
$n = 0; # initialize couner
my $reads_num; # number of reads
my $progress; # n million reads
my %index_hash; # key is c_barcode, val is index
my %index_hash_reverse; # key is index, val is c_barcode
my $split_barcode_num; # num of represented c_barcodes
my $T;
my $id;
my $reads_num;
my @line;
my @Read_num;
$Read_num[0] = 0;
my $split_reads_num;
my $bc;
my $bx;
my $read;
my $additional_bc;

while(<IN2>){ # read lines from Read2
  chomp; # get line
  @line = split; # split line (this seems unnecessary)
  $n ++; # counter += 1
  if($n % 4 == 1){ # if counter modulo 4 = 1
    $reads_num ++; # add 1 to the number of reads processed
    my @A  = split(/\//,$line[0]); # split line by forward slash, save as A
         $id = $A[0]; # id is the readname
         if($reads_num % 1000000 == 1) # if 1mil reads
         {
              print "reads processed $progress (M) reads ...\n"; # print message
              $progress ++; # add 1 mil reads for next time
         }
  }

  if($n % 4 == 2){ # if we have the second line (sequence)

    if($bc_pos==5){
      $bc = substr($line[0], 0, $bc_len_redundant);
      $bx = substr($line[0], 0, $bc_len_select);
      $read = substr($line[0], $bc_len_redundant+$bc_trim_gdna, $read_len);
      }elsif($bc_pos==3){
      ## start from index skip, get length of read_len
      $read = substr($line[0], $skip, $read_len);
      $bc = substr($line[0], $skip+$read_len, $bc_len_redundant);
      $bx = substr($line[0], $skip+$read_len, $bc_len_select);
      $additional_bc = substr($line[0], 0, $skip);
      }

    if(!(exists $index_hash{$bx})){ # if the combined id isn't in the combined hash map
      $split_barcode_num ++; # add another represented barcode
      $index_hash{$bx} = $split_barcode_num; # add the combined id to the hash map
      $index_hash_reverse{$split_barcode_num} = $bx; # add the reverse (key is n combined barcodes, value is c_bc)
      $Read_num[$index_hash{$bx}] = 0; # this is an array of num reads for each represented barcode
                                          # because the barcodes represented start at 0 they can act as the index
    }
    $split_reads_num ++; # add a split read
    $Read_num[$index_hash{$bx}] ++;

      $T = <IN1>; chomp($T); # remove new line
      print OUT1 "$id\#$bc\tBX:Z:$bx\n";
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n"; # JUST PRINT
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n"; # JUST PRINT
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n"; # JUST PRINT

    if($skip_as_bc){
        print OUT2 "$id\#$bc\#$additional_bc\tBX:Z:$bx\n";
      } else {
        print OUT2 "$id\#$bc\tBX:Z:$bx\n";
      }
    print OUT2 "$read\n";
    $T = <IN2>; $n++;chomp($T);
    print OUT2 "$T\n";
    $T = <IN2>; $n++;chomp($T);
    my $qual = substr($T,0,$read_len);
    print OUT2 "$qual\n";
  }

}
close IN1;
close IN2;
close OUT1;
close OUT2;

my $barcode_types = 4**$bc_len_select;

open OUT3, ">split_stat_read1.log" or die "Can't write file";
print OUT3 "Barcode_types = 4** $bc_len_select = $barcode_types\n";
my $r;
$r = 100 *  $split_barcode_num/$barcode_types;
print OUT3 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 *  $split_barcode_num/$reads_num;
print OUT3 "Reads_pair_num  = $reads_num \n";
print OUT3 "Reads_pair_num(after split) = $split_barcode_num ($r %)\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
  print OUT3 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}
close OUT3;

### check diversity
open OUT4, ">bc.diversity.count.log" or die "Can't write file";

my $bc_diversity_split = $split_barcode_num/$reads_num;
print OUT4 "$bc_diversity_split\n";
close OUT4;

print "all done!\n";
