use strict;
use warnings;

=comm
AUTHOR: R.PRADEEP
INSTITUTE: ICRISAT
USAGE: perl ParsePedMap4PairwiseSNPMatrix.pl SNP.plink.ped SNP.plink.map group.txt
OUTPUT: Generates the inputfile.Matrix.txt and Genotype.fasta file have each sample print as fasta sequence
Based on the group ids, the output matrix will have relatvie order of ids

group.txt
IS1004  1
IS10302 1
IS1212  1
IS1219  2
IS1233  1
IS12697 1
IS12937 2
IS12945 2
IS13294 2
IS13444 2
IS13549 1
IS13809 3
IS15401 3
IS15466 3
IS15744 3
IS15845 1
IS15931 4
IS18484 4
IS18758 4

=cut

    my $ped = $ARGV[0];
    open(my $ped_fh, '<:encoding(UTF-8)', $ped)
      or die "Could not open file '$ped' $!";
      
      
my @map;
my %group=();
my %getalleledata=();
my %acc1=();
my %acc2=();
my %comb=();

my %getdata=();
my %panids=();
my %deleteSNPs=();
my %countGtyp4SNP=();
my $accCount=1;

    open(my $log_fh, '>:encoding(UTF-8)', "$ped.Matrix.txt")
      or die "Could not open file '$ped.log' $!";


#getmap();
getgroup();

    open(my $fasta_fh, '>:encoding(UTF-8)', "Genotype.fasta")
      or die "Could not open file 'Genotype.fasta' $!";
print $log_fh "Parsing ped input file with samples..\n";
    while (my $row = <$ped_fh>) {
		#next if ( $row =~ m/^#/xms );
      chomp $row;
		my @elements = split( /\s/, $row );
        my $sample=$elements[1];
        
        my $marknum=0;
        my $seq="";
        if(exists $group{$sample}){
            for(my $i=6; $i<=$#elements; $i++){
                my $allele1=$elements[$i];
                #my $snpindex=($i-6)-$marknum;
                #my $snpid=$map[$snpindex];
                $i++;
                my $allele2=$elements[$i];
                $marknum++;
                #if($allele1 ne "0"){
                #$getalleledata{$snpid}{$sample}=$allele1;
                #push @{$getalleledata{$sample}},$allele1;
                my $groupid=$group{$sample};
                $acc1{$groupid}{$sample}++;
                $acc2{$groupid}{$sample}++;            
                #}
                $seq=$seq.$allele1;
                #print "$snpid\t$sample\t$allele1\n";
                #print "[$i]-";
    
            }
            $getalleledata{$sample}=$seq;
            #print "\n$seq $sample\n";
            print $log_fh "[$accCount] ";
            $seq=~s/0/N/g;
            print $fasta_fh ">$sample\n$seq\n";
            $accCount++;
        }
    }
close($ped_fh);

print $log_fh "\nFinished parsing input file\n";
$accCount=1;

#for my $snpid ( keys %getalleledata ) {
    foreach my $c1group (sort keys %acc1){
	foreach my $cul1 (sort keys %{$acc1{$c1group}}) {
        print $log_fh "[$accCount] ";
        foreach my $c2group(sort keys %acc2){
    		foreach my $cul2 (sort keys %{$acc2{$c2group}}) {
			if(exists $getalleledata{$cul1} && exists $getalleledata{$cul2} ){
				my $combination=$cul1."~~!!~~".$cul2;
                my @arr1=split(//,$getalleledata{$cul1});
                my @arr2=split //,$getalleledata{$cul2};
                #print scalar(@arr2).length($getalleledata{$cul2})."  ===\t===".scalar(@arr1).length($getalleledata{$cul1})."\n";
                #print "Now $combination is in progress\n";
                for(my $indx=0; $indx <= $#arr1 || $indx <= $#arr2; $indx++){ #my $indx (0 .. $#arr2){
                    #print "$arr1[$indx] $arr2[$indx]\n";
                    if((exists $arr1[$indx] && exists $arr2[$indx]) &&($arr1[$indx]=~m/[A|T|G|C]/ ) && ($arr2[$indx]=~m/[A|T|G|C]/) && $arr1[$indx] ne $arr2[$indx]){
                        $comb{$combination}++;
                    }
                }
				
				#print "snpid:$snpid $cul1\t$cul2\t$getalleledata{$snpid}{$cul1}\t$getalleledata{$snpid}{$cul2}\n";

			}
		}
        }
        $accCount++;
	} 
    }
#}
print $log_fh "\n";
print $log_fh "Finished assessing SNP count for all combinations\n";


# Print final output matrix
foreach my $c1group (sort keys %acc1){
	foreach my $cul1 (sort keys %{$acc1{$c1group}}) { 
		print $log_fh "\t$cul1";
	}
}
print $log_fh "\n";

my $printflag="off";
foreach my $c1group (sort keys %acc1){
	foreach my $cul1 (sort keys %{$acc1{$c1group}}) { 
		print $log_fh "$cul1\t";
        foreach my $c2group(sort keys %acc2){
    		foreach my $cul2 (sort keys %{$acc2{$c2group}}) {

				my $combination=$cul1."~~!!~~".$cul2;
                my $combination2=$cul2."~~!!~~".$cul1;
			#if(exists $comb{$combination} && $printflag eq "on"){ # Tringle matrix
            if(exists $comb{$combination} ){  # Square matrix  
				print $log_fh "$comb{$combination}\t";
				#print "snpid:$snpid $cul1\t$cul2\t$getalleledata{$snpid}{$cul1}\t$getalleledata{$snpid}{$cul2}\n";
			}elsif(exists $comb{$combination2}){ # For square matrix
                print $log_fh "$comb{$combination2}\t";
            }
            else{
				print $log_fh "0\t";
			}

			if($cul1 eq $cul2){$printflag="on";}
            } # cult2 loop
        }
		print $log_fh "\n";
		$printflag="off";
	} #cult1 loop
}


sub getgroup{
        my $g_file = $ARGV[2];
    open(my $gp_fh, '<:encoding(UTF-8)', $g_file)
      or die "Could not open file '$g_file' $!";
      while (my $row=<$gp_fh>){
        my @elements = split( /\t/, $row );
        $group{$elements[0]}=$elements[1];
      }
}

sub getmap{

        my $map = $ARGV[1];
    open(my $map_fh, '<:encoding(UTF-8)', $map)
      or die "Could not open file '$map' $!";
      
        while (my $row = <$map_fh>) {
         #next if ( $row =~ m/^#/xms );
        chomp $row;
		my @elements = split( /\t/, $row );
        my $id=$elements[0].":".$elements[3];
        push(@map,$id);
        }
}
