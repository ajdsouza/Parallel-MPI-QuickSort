#/usr/bin/perl

use strict;
use warnings;

use File::Spec;
use List::Util qw(first);

my %ld;
my @results;
my @errors;

die " provide directory path as arg1 perl plotgraph <dir,dir,dir>\n" unless @ARGV > 0;

my $dirs = $ARGV[0];
my @dirs = split/,/,$dirs;

for my $dir ( @dirs ) {

opendir(DIR, $dir) or die $!;

#p-1-n-10000.txt

my @files 
        = grep { 
            /^p\-\d+\-n\-\d+.txt/             # p-16-n-300.txt
	    && -f "$dir/$_"   # and is a file
	} readdir(DIR);

closedir(DIR);

# Loop through the array printing out the filenames

     foreach my $file (@files) {
       
        $file =~ s/\s+// if $file;

        my ($p,$n) = ( $file =~ /^p\-(\d+)\-n\-(\d+).txt/ );

        open (FH, File::Spec->catfile($dir,$file)) || die "ERROR Unable to open file: $!\n";

	my $last;
	my $first = <FH>;

	while (<FH>) { if ( $_ =~ /\d+\.\d+/ ){ $last = $_; last;}}
	close FH;

	# for one line files
        $last=$first unless $last;

        $last =~ s/^\s+|\s+$// if $last;
        chomp $last if $last;
        $last =~ s/^\s+|\s+$// if $last;

        my $msg = '';
 
	my $tm=-1;
        unless ($last and $last =~ /\d+\.\d+/) { # sorting took: 0.117762 s
         $msg = $last if $last;
         $last = -1;
         push @errors,"$file=$p,$n=$last,\"$msg\"";
        }
        else {
         ( $tm ) = ( $last =~ /\s*(\d+\.\d+) s\s*/ );
	 print "$file $last=$tm\n" if $tm > 5;
         push @results, "$file=$p,$n=$tm,\"$msg\"";
        }
	print "$file - $last - $tm\n" unless $tm > 0;
        push @{$ld{p}{$p}{$n}},$tm;
        push @{$ld{n}{$n}{$p}},$tm;

    }
}

my $plot = "gplot.txt";
open(FH, '>>', $plot) or die "Failed to open file $plot for writing \n";

for my $res (@results){
   print FH "$res\n";
}
for my $err (@errors){
   print FH "$err\n";
}
close(FH);


# create a speed up hash type too
my %templd;
for my $type ( sort keys %ld ){

   for my $typeval ( sort keys %{$ld{$type}} ) {

    for my $ps ( sort keys %{$ld{$type}{$typeval}} ) {

         my $speedup=-1;
         my $p1speed=-1;
         my $bestspeed=-1;

         if ( ( $type =~ /^p$/ ) and ( $ld{$type}{1}{$ps} ) ){ $p1speed = first { $_ > 0 }  sort {$a <=> $b} @{$ld{$type}{1}{$ps}}; }
         if ( ( $type =~ /^n$/ ) and ( $ld{$type}{$typeval}{1} ) ){ $p1speed = first { $_ > 0 }  sort {$a <=> $b} @{$ld{$type}{$typeval}{1}}; }

         $bestspeed = first { $_ > 0 } sort {$a <=> $b} @{$ld{$type}{$typeval}{$ps}};
         $bestspeed=-1 unless $bestspeed;

         push @{$templd{"speedup_$type"}{$typeval}{$ps}},($p1speed/$bestspeed) if $p1speed and $p1speed > 0 and $bestspeed > 0;

     }
   }
}


for my $type ( sort keys %templd ){
 $ld{$type}=$templd{$type};
}

my %tags = ('p'=>'n','n'=>'p',,'speedup_p'=>'n','speedup_n'=>'p');

# type = p
for my $type ( sort keys %ld ){

   my @pso;
# typeval p=4 etc
   for my $typeval ( sort keys %{$ld{$type}} ) {
     push @pso,keys %{$ld{$type}{$typeval}};
   }
# all the n values
   my %psoh = map { $_ => 1 } @pso ;

   my $plottypeval = "gplot$type.txt";
   open(FH, '>>', $plottypeval) or die "Failed to open file $plottypeval for writing \n";

   
# p,,
  print FH "$type - $tags{$type},";
#  print FH ",";

# n=100
# # sort the column header numeric ascending
  for my $p (  sort { $a <=> $b} keys %psoh  ){
   print FH "$p,";
  }
  print FH "\n";

       
 for my $typeval ( sort {$a <=> $b} keys %{$ld{$type}}){

    print FH "$type=$typeval,";
    
    for my $p ( sort { $a <=> $b} keys %psoh ){

     if ( $ld{$type}{$typeval}{$p} and @{$ld{$type}{$typeval}{$p}}){
      # get the best time > 0
      my $tm=-1;
      if ( $type =~ /speedup/ ) {
       $tm = first { $_ > 0 } sort {$b <=> $a} @{$ld{$type}{$typeval}{$p}};
      } else {
       $tm = first { $_ > 0 } sort {$a <=> $b} @{$ld{$type}{$typeval}{$p}};
      }
      $tm='' unless $tm;
      print FH "$tm,";
     }
     else {
       print FH ",";
     }

   }

   print FH "\n";

 }

  close(FH);

}



exit 0;
