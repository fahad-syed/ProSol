#!/usr/bin/perl

=head1 NAME

  run2wayBlast.pl

=head1 DESCRIPTION

  For running a 2 way blast from a list of file pairs

=head1 ARGUMENTS

Usage:

=head1 Authors

 Mikko Arvas
=cut
use strict;
use warnings;
use Cwd;
use File::Basename;
use Getopt::Long;

my $input='';
GetOptions (
	    'i=s' => \$input,
	   );
my $usage = "Usage:\n-i input list file with structure:\nfile1\tfile1Type\tfile2\tfile2Type\toutFile \n";
unless ($input) {
  die print "$usage";
}

#check if list file exists
unless(-e $input){
  print STDERR "File $input doesn't exist!\n$usage\n ... ";
  exit -1;
}
open INFILE, "<$input";
#read in list file
my $rowc=0;
my %intable = ();
while (<INFILE>) {
  ++$rowc;
  my $row = $_;
  my @rowa=split("\t",$_);
  #print "@rowa\n";
  #parse names and file types
  $intable{'file1'} = $rowa[0];
  $intable{$rowa[0].'Type'} = $rowa[1];
  $intable{'file2'} = $rowa[2];
  $intable{$rowa[2].'Type'} = $rowa[3];
  $intable{'outFile'} =$rowa[4];
  chomp($intable{'outFile'});
  my @files= ($intable{'file1'}, $intable{'file2'});
  my @dbfiles =('hr','in','sd','si','sq');
  foreach my $f (@files) {
    #check that all the files exist 
    unless (-e $intable{'file1'}) {
      print STDERR "File $intable{'file1'} doesn't exist!\n$usage\n ... ";
      exit -1;
    } else {
      #check if data base has to be done and make it if necessary
      my $datatype = 'p';
      if ($intable{$f.'Type'} eq 'P' ) {
      } elsif ($intable{$f.'Type'} eq 'N' ) {
	$datatype = 'n';
      }
 #     print "$f\n";
      my @dbfilesR = ();
      my $nodb = 1;
      for my $i (@dbfiles) {
	(my $pi = $i)=~s/^/$datatype/;
	unless (-e $f.".".$pi) {
	  $nodb=0;
	  print "$f.$pi blast db does not exists\n";
	  last;
	}
      }
      my $pparameter = "T";
      if ($datatype eq 'n') {
	$pparameter = "F";
      }
      if (! $nodb) {
	my $cmd = "formatdb -i $f -p $pparameter  -o T ";
	print "fCMD: $cmd\n";
	system($cmd) == 0 or die "system \'$cmd\' failed: $?\n"; 
      }
    }
  }
  #once dbs are done run the blast
  #check if 1vs2 exist or run
  if ($intable{$intable{'file1'}.'Type'} eq 'P' && $intable{$intable{'file2'}.'Type'} eq 'P' ) {
    if (! -e $intable{'outFile'}.'.1V2' ) {
      my $cmd = " blastall -p blastp -d $intable{'file2'} -i $intable{'file1'} -o  $intable{'outFile'}.1V2 -m 8 -b 1 -v 1  -e 1-e5";
      print "bCMD: $cmd\n";
      system($cmd) == 0 or die "system \'$cmd\' failed: $?\n";
    }
    #check if 2vs1 exist or run
    if (! -e $intable{'outFile'}.'.2V1' ) {
      my $cmd = " blastall -p blastp -d $intable{'file1'} -i $intable{'file2'} -o  $intable{'outFile'}.2V1 -m 8 -b 1 -v 1 -e 1-e5";
      print "bCMD: $cmd\n";
      system($cmd) == 0 or die "system \'$cmd\' failed: $?\n";
    }
  } else {
    print "$intable{'file1'}: $intable{$intable{'file1'}.'Type'} and $intable{'file2'}: $intable{$intable{'file2'}.'Type'} -> unknown types\n"
  }
}

exit;


