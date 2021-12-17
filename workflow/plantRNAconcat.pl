#!/usr/bin/env perl

##########################################################################################
# Name:		plantRNAconcat.pl
# Author:	Valerie Cognat valerie.cognat@ibmp-cnrs.unistra.fr
# Copyrigth :	IBMP - CNRS
# Description:	Concatenate tRNAscan output, secondary structure tRNAscan output and
#               upstream and downstream sequence
#
# Created : 19 november 2012 - version 1
# Last modification : :
# - 12 february 2020 - version 2 (trnascanSE 2.0.5 + seqKit)
# - 25 March 2020 - version 3 (extract specific columns in tRNAscanSE output, can differ if using -H option)
# - 12 April 2021 - version 3.1 (Add p-1 empty column)
# - 06 September 2021 - version 3.2 (correction of polyT pattern to highlight more than 4 T)
##########################################################################################

use strict;
use Getopt::Long;


########################
## Loading parametres ##
my ($trna, $struct, $splitstruct, $up, $down, $out);
my $usage = "[USAGE] : perl plantRNAconcat.pl -trna <tRNAscan-SE file> -struct <tRNAscan-SE structure file> -split <split tRNA struct file> -up <upstream sequence file> -down <downstream sequence file> -out <csv outfile>\n";
my $optparse = GetOptions ('trna=s' => \$trna,
				'struct:s' => \$struct,
                'split:s' => \$splitstruct,
				'up:s' => \$up,
                'down:s' => \$down,
                'out:s' => \$out);
unless ($trna && $struct && $splitstruct && $up && $down && $out) {
    print $usage."\n";
    exit;
}

open (TRNASCAN, "$trna") || die "\t** Error : can't open $trna -> Exit **\n";
open (TRNASTRUCT, "$struct") || die "\t** Error : can't open $struct -> Exit **\n";
open (SPLITSTRUCT, "$splitstruct") || die "\t** Error : can't open $splitstruct -> Exit **\n";
open (UP, "$up") || die "\t** Error : can't open $up -> Exit **\n";
open (DOWN, "$down") || die "\t** Error : can't open $down -> Exit **\n";
open (OUT, ">$out") || die "\t** Error : can't create $out -> Exit **\n";

open (ERR, ">$out.err") || die "\t** Error : can't create $out.err -> Exit **\n";

my %info; # hash to stock informations about each tRNA


#--------------------------------------
# tRNAscan output treatment
#
# Header : 3 lines
my $tRNAline = <TRNASCAN>;
my $tRNAline = <TRNASCAN>;
my $tRNAline = <TRNASCAN>;
my @colnb = split(/\t/,$tRNAline); # calculate the number of column with the header
print $#colnb;
#print OUT $tRNAline;
while ($tRNAline = <TRNASCAN>) {
    chomp ($tRNAline);
    $tRNAline =~ s/ //g; # supress space used for tRNAscan-SE paging
    #print $tRNAline."\n";
    next if ($tRNAline eq "");
    my @tab = split(/\t/,$tRNAline);
    if (defined ($info{$tab[0].".trna".$tab[1]})) {
        print ERR "Warning tRNA: ".$info{$tab[0].".trna".$tab[1]}." already exists in hash, check $trna for $tRNAline.\n";
    }
    else {
        $tRNAline =~ s/\s/;/g; ##  csv format
        $info{$tab[0].".trna".$tab[1]} = {};
        $info{$tab[0].".trna".$tab[1]}{'tRNA'} = join(";",@tab[0..8]); # take first 9 columns
        $info{$tab[0].".trna".$tab[1]}{'tRNA'} .= ";".@tab[$#colnb];  # add Note in last column
        #print "* $tRNAline\n";
        #print "* ".$info{$tab[0].".trna".$tab[1]}{'tRNA'}."\n";
        #print ERR $tab[0].".trna".$tab[1]." line : ".$tRNAline."\n";
    }
}
close TRNASCAN;


#--------------------------------------
# tRNA struct treatment
#
my ($istart, $istop) = (0,0); #
#print OUT "\tsecondary_structure";
while (my $structline = <TRNASTRUCT>) {
    chomp ($structline);
    next if ($structline eq "");
    #print $structline."\n";
    my @tab = split(/\s/,$structline);
    while ($structline !~ /^Seq:/) {
        if ($structline =~ /^Possible intron:/){
            $structline =~ /Possible intron: (\d+)-(\d+)\s/;
            $istart = $1;
            $istop = $2;
        }
        $structline = <TRNASTRUCT>;
    }
    my @str = split(/\s/,$structline);
    if (! defined ($info{$tab[0]})) {
        print ERR "** Warning SEQ : ".$info{$tab[0]}." doesn't exist in hash, check $struct for $structline. **\n";
    }
    else {
        if (($istart != 0 ) && ($istop != 0)){ # epissage de l'intron
            $info{$tab[0]}{'intron'} = substr($str[1], $istart-1, ($istop-$istart+1));
            $info{$tab[0]}{'seq'}  = substr($str[1], 0, $istart-1);
            $info{$tab[0]}{'seq'}  .= substr($str[1], $istop, length($str[1])-$istop+1 );
            #print "INTRON: ".$info{$tab[0]}{'intron'} ."\nEPISS: ".$info{$tab[0]}{'seq'} ."\n";
        }
        else {
            $info{$tab[0]}{'seq'} = $str[1];
        }
    }
    while ($structline !~ /^Str:/) {
        $structline = <TRNASTRUCT>;
    }
    my @str = split(/\s/,$structline);
    if (! defined ($info{$tab[0]})) {
        print ERR "** Warning STRUCT : ".$info{$tab[0]}." doesn't exist in hash, check $struct for $structline. **\n";
    }
    else {
        if (($istart != 0 ) && ($istop != 0)){ # epissage de l'intron dans la struct
            $info{$tab[0]}{'struct'}  = substr($str[1], 0, $istart-1);
            $info{$tab[0]}{'struct'}  .= substr($str[1], $istop, length($str[1])-$istop+1 );
        }
        else {
            $info{$tab[0]}{'struct'} = $str[1];
        }
    }
    ($istart, $istop) = (0,0); # reinitialisation
}
close TRNASTRUCT;


#--------------------------------------
# Upstream sequence treatment
#
#print ERR "\tupstream\n";
while (my $upline = <UP>) {
    chomp ($upline);
    next if ($upline eq "");
    if ($upline =~ ">") { # fasta header
        my @fastaname = split(/ /,$upline);
        my @trnaid = split(/-/,$fastaname[1]);
        $trnaid[0] =~ s/tRNA/trna/g;  ## WARNING, in tRNAscan bed file it is not in lower case as in struct file
        #print $trnaid[0]."\n";
        if (! defined ($trnaid[0])) {
            print ERR "Warning UP : ".$trnaid[0]." doesn't exist in hash, check $up for $upline.\n";
        }
        else {
            $upline = lc(<UP>); # next line to have up seq
            chomp($upline);
            $upline =~ s/caa/CAA/g;
            $info{$trnaid[0]}{'up'} = $upline;
            #print $trnaid[0]." up : *".$upline."*\n";
        }
    }

}
close UP;


#--------------------------------------
#  Downstream sequence treatment
#
#print OUT "\tdownstream\n";
while (my $downline = <DOWN>) {
    chomp ($downline);
    next if ($downline eq "");
    if ($downline =~ ">") { # fasta header
        my @fastaname = split(/ /,$downline);
        my @trnaid = split(/-/,$fastaname[1]);
        $trnaid[0] =~ s/tRNA/trna/g; ## WARNING, in tRNAscan bed file it is not in lower case as in struct file
        #print $trnaid[0]."\n";
        if (! defined ($trnaid[0])) {
            print ERR "Warning DOWN : ".$trnaid[0]." doesn't exist in hash, check $up for $downline.\n";
        }
        else {
            $downline = lc(<DOWN>); # next line to have up seq
            chomp($downline);
						if ($downline =~ m/(t{4,})/) {  # if polyT, upper case the pattern
							my $polyT= uc($1);
							$downline =~ s/$1/$polyT/;
						}
						print "$downline \n";
            $info{$trnaid[0]}{'down'} = $downline;
            #print $trnaid[0]." down : *".$downline."*\n";
        }
    }

}
close DOWN;

#--------------------------------------
# tRNA split structure treatment
#
while (my $splitline = <SPLITSTRUCT>) {
    chomp($splitline);
    next if ($splitline eq "");
    next if ($splitline =~ "^id");
    # id,Acc-Stem1,P8-9,D-stem1,D-loop,D-stem2,p26,Ac-stem1,Ac-loop,Ac-stem2,V-region,T-stem1,T-loop,T-stem2,Acc-stem2,p73
    my ($id, $accStem1, $p8, $Dstem1, $Dloop, $Dstem2, $p26, $AcStem1, $AcLoop, $AcStem2, $Vregion, $Tstem1, $Tloop, $Tstem2, $AccStem2, $p73) = split(/,/,$splitline);
    #print "$id * $accStem1 * $p8 * $Dstem1 * $Dloop * $Dstem2 * $p26 * $AcStem1 * $AcLoop * $AcStem2 * $Vregion * $Tstem1 * $Tloop * $Tstem2 * $AccStem2 * $p73**\n";
    $info{$id}{'accStem1'} = $accStem1;
    $info{$id}{'p8-9'} = $p8;
    $info{$id}{'Dstem1'} = $Dstem1;
    $info{$id}{'Dloop'} = $Dloop;
    $info{$id}{'Dstem2'} = $Dstem2;
    $info{$id}{'p26'} = $p26;
    $info{$id}{'AcStem1'} = $AcStem1;
    $info{$id}{'AcLoop'} = $AcLoop;
    $info{$id}{'AcStem2'} = $AcStem2;
    $info{$id}{'Vregion'} = $Vregion;
    $info{$id}{'Tstem1'} = $Tstem1;
    $info{$id}{'Tloop'} = $Tloop;
    $info{$id}{'Tstem2'} = $Tstem2;
    $info{$id}{'AccStem2'} = $AccStem2;
    $info{$id}{'p73'} = $p73;
}
close SPLITSTRUCT;

#--------------------------------------
# output
#
print OUT "SequenceName;tRNA#;tRNA_Bounds_Begin;tRNA_Bounds_End;tRNA_Type;Anti_Codon;Intron_Bounds_Begin;Intron_Bounds_End;Cove_Score;Note;upstream;downstream;Sequence;Intron;secondary_structure;p-1;Acc-Stem1;p8-9;D-stem1;D-loop;D-stem2;p26;Ac-stem1;Ac-loop;Ac-stem2;V-region;T-stem1;T-loop;T-stem2;Acc-stem2;p73\n";
for my $trnaseq ( keys %info ) {
    print OUT $info{$trnaseq}{'tRNA'}.";".$info{$trnaseq}{'up'}.";".$info{$trnaseq}{'down'}.";".$info{$trnaseq}{'seq'};
    #print $trnaseq."--".$info{$trnaseq}{'tRNA'}.";".$info{$trnaseq}{'up'}.";".$info{$trnaseq}{'down'}.";".$info{$trnaseq}{'seq'}."\n";
    if (defined $info{$trnaseq}{'intron'}){
        print OUT ";".$info{$trnaseq}{'intron'}; # else empty col
    }
    else {
        print OUT ";";
    }
    print OUT ";".$info{$trnaseq}{'struct'};
		# Add empty column for p-1 before accStem1
    print OUT ";;".$info{$trnaseq}{'accStem1'}.";".$info{$trnaseq}{'p8-9'}.";".$info{$trnaseq}{'Dstem1'}.";".$info{$trnaseq}{'Dloop'}.";".$info{$trnaseq}{'Dstem2'};
    print OUT ";".$info{$trnaseq}{'p26'}.";".$info{$trnaseq}{'AcStem1'}.";".$info{$trnaseq}{'AcLoop'}.";".$info{$trnaseq}{'AcStem2'}.";".$info{$trnaseq}{'Vregion'};
    print OUT ";".$info{$trnaseq}{'Tstem1'}.";".$info{$trnaseq}{'Tloop'}.";".$info{$trnaseq}{'Tstem2'}.";".$info{$trnaseq}{'AccStem2'}.";".$info{$trnaseq}{'p73'};
    print OUT"\n";
}

close OUT;

close ERR;

exit (0);
