#!/usr/bin/perl

##########################################################################################
# Name:		tRNA2gff.pl
# Author:	Valerie Cognat valerie.cognat@ibmp-cnrs.unistra.fr
# Copyrigth :	IBMP - CNRS
# Description:	convert tRNA csv file to tRNA gff file
#
# Created : 8 november 2012
# Version 1.0 - 18 november 2012
# Version 2.0 - 2019, July, 3rd
#		- add mpo genome
#		- correct ID if -ref acc and acc begins by chr or Chr
#		- correct tRNA coordinates
# Version 2.1 - 2019, July the 8th
#		- change sp to string for gff3 file
# Version 2.2 - 2019 , July the 8th
#		- change gff format to have gene / tRNA and intron => add -ft option (all / gene)
#		- add sly genome
# Version 2.3 - 2019, July the 10th
#		- add an option for plantRNA xls (csv) to gff (header change)
#		- add tRNA-Sec
# Version 2.4 - 2019, July the 23rd
#		- modify IUPAC code to add Sec and Pyl (pyrrolysine)
#       - case insensitive for Meti/MetI
#       - usage correction (2020/03/20)
# Version 2.5 - 2020, April the 2nd
#       - Add sequence-region information
# Version 2.6 - 2020 July the 17th
#       - correction of sequence name if not chromosome
# Version 2.7 - 2021 Februar the 11th
#       - correction of sequence name if not chromosome or number
#				- filter sequence-region to have only accession present in the gff3
# Version 2.8 - 2021 October the 29th
#				- add eMet & iMet detection
##########################################################################################

use strict;
use Getopt::Long;
use File::Temp qw/ :POSIX /;

my $version = "version 2.8 - 2021/10/29";

########################
## Loading parametres ##
my ($gff, $csv, , $alias, $specie);
my ($build, $ref, $ft) = ("", "acc", "all"); # default value
my $usage = "[USAGE] : perl tRNA2gff.pl -csv <csv_file> -gff <gff_file> -sp <specie> -build <genome build version> -ref <acc or chr> -ft <feature_type> -alias <genome_alias_file>
\t\t-sp = name of the species
\t\t-build = genome assembly build version - format : source buildName (ex: NCBI B36)
\t\t-ref = acc for accession number or chr for chromosome name # default = acc
\t\t-ft = all or gene (for all features or gene only) # default = all
\t\t-alias = alias file contains accession chr sequence_length - the file is mandatory in trnaflow
[VERSION] : $version\n";
my $optparse = GetOptions ('gff=s' => \$gff,
				'csv=s' => \$csv,
				'sp=s' => \$specie,
                'build:s' => \$build,
                'ref:s' => \$ref,
				'ft:s' => \$ft,
                'alias=s' => \$alias);
unless ($specie && $gff && $csv && $ref && $ft && $alias) {
    die $usage."\n";
}
if ($ref ne "acc" && $ref ne "chr") {
    die "\t** Error : -ref option have to be set to 'acc' or 'chr' **\n";
}
if ($specie eq "" ) {
    die "\t** Error : -sp option have to be set to genome description**\n";
}
if ($ft ne "all" && $ft ne "gene") {
    die "\t** Error : -ft option have to be set to 'all' or 'gene' **\n";
}
unless (-r $csv) {
print "Error, can't read $csv file\n";
}
unless (-r $alias) {
    print "Error, can't read $alias file\n";
}

open(GFF, ">$gff") || die "Error, can't create $gff file\n";
open(CSVhead, "$csv") || die "Error, can't open $csv file\n";

# Traitement entête CSV pour avoir indice des colonnes
my $line = <CSVhead>;
chomp($line);
my @colname = split (/;/, $line);
# indices des colonnes
my ($acci, $starti, $stopi, $typei, $aci, $intronbi, $intronei, $introntypei, $chri); # = ("NA","NA","NA","NA","NA","0","0","NA","NA");
for (my $i=0; $i <= $#colname; $i++) {
	#print $i."*".lc($colname[$i])."*\n";
    if (lc($colname[$i]) eq 'chr') {
        $chri = $i;
    }
    elsif ($colname[$i] eq "SequenceName") {
        $acci = $i;
    }
    elsif (($colname[$i] eq "tRNA_Bounds_Begin") || ($colname[$i] eq "Start")){
        $starti = $i;
    }
    elsif (($colname[$i] eq "tRNA_Bounds_End") || ($colname[$i] eq "Stop")){
        $stopi = $i;
    }
    elsif (($colname[$i] eq "tRNA_Type") || ($colname[$i] eq "AA")){
        $typei = $i;
    }
    elsif (($colname[$i] eq "Anti_Codon") || ($colname[$i] eq "AC")){
        $aci = $i;
    }
    elsif ($colname[$i] eq "Intron_Bounds_Begin") {
        $intronbi = $i;
    }
    elsif ($colname[$i] eq "Intron_Bounds_End") {
        $intronei = $i;
    }
    elsif ($colname[$i] eq "Intron_type") {
        $introntypei = $i;
    }
	# else {
	# 	print $i."**".$colname[$i]."**\n";
	# }
}
#print "$acci * $starti * $stopi * $typei * $aci * $intronbi * $intronei * $introntypei * $chri\n";
if ($chri eq "NA") { # sometimes, depends of the preprocessing step
	$chri = $acci;
}
if ($acci eq "NA") { #plantRNA output
	$acci = $chri;
}
close(CSVhead);

# TREATMENT OF ALIAS FILE
my %refsize; # hash to store size of each sequence (chr, scaffold)
if ($ref eq "acc") {
    my $comm = "tail -n +2 $alias | sort -V -k1,1 > $alias.sort";
	system ($comm) == 0 || die "** Error on $comm **\n";
    open(ALIAS, "$alias.sort") || die "Error, can't open $alias file\n";
    while (my $aline = <ALIAS>) {
        chomp($aline);
        my ($accName, $chrName, $size) = split(/ /, $aline);
        $refsize{$accName} = $size;
    }
    close(ALIAS);
}
else { # chrName
    my $comm = "tail -n +2 $alias | sort -V -k2,2 > $alias.sort";
	system ($comm) == 0 || die "** Error on $comm **\n";
    open(ALIAS, "$alias.sort") || die "Error, can't open $alias file\n";
    while (my $aline = <ALIAS>) {
        chomp($aline);
        my ($accName, $chrName, $size) = split(/ /, $aline);
        $refsize{$chrName} = $size;
    }
    close(ALIAS);
}

# GFF header
# ## is for specific metadata / # is for comments
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = (1900+$year)."-".(1+$mon)."-".$mday;
print GFF "##gff-version 3\n";
print GFF "##species $specie\n";
print GFF "##genome-build $build\n"; #source buildName
print GFF "#Chromosomal coordinates of tDNA extracted from PlantRNA DB http://plantrna.ibmp.cnrs.fr/\n";
print GFF "#date: $date\n";

my %typeNB;
$typeNB{"A"} = 0;
$typeNB{"C"} = 0;
$typeNB{"D"} = 0;
$typeNB{"E"} = 0;
$typeNB{"F"} = 0;
$typeNB{"G"} = 0;
$typeNB{"H"} = 0;
$typeNB{"I"} = 0;
$typeNB{"J"} = 0;
$typeNB{"K"} = 0;
$typeNB{"L"} = 0;
$typeNB{"M"} = 0;
$typeNB{"Mi"} = 0;
$typeNB{"Me"} = 0;
$typeNB{"N"} = 0;
$typeNB{"P"} = 0;
$typeNB{"Q"} = 0;
$typeNB{"R"} = 0;
$typeNB{"S"} = 0;
$typeNB{"T"} = 0;
$typeNB{"U"} = 0;
$typeNB{"V"} = 0;
$typeNB{"W"} = 0;
$typeNB{"Y"} = 0;
$typeNB{"U"} = 0;
$typeNB{"O"} = 0;
$typeNB{"Sup"} = 0;

## CSV file conversion and extract unique sequence-region (chr ou acc)
# sort csv file - (keep header) - sort - write csv.sort
#my $comm = "head -n 1 $csv > $csv.sort && tail -n +2 $csv | sort -f -t ';' -k1,1 -k3n,3n >> $csv.sort";
my $clisort;
my $cliuniq;
if ($ref eq "chr") { # sort by chromosome number
	$clisort = "tail -n +2 $csv | sort -f -t ';' -V -k".($chri+1).",".($chri+1)." -k".($starti+1).",".($starti+1)."n > $csv.sort";
	$cliuniq = "cut -d ';' -f1 $csv.sort | uniq > $csv.uniqregion";
}
elsif ($ref eq "acc") { # sort by accession number
	$clisort = "tail -n +2 $csv | sort -f -t ';' -V -k".($acci+1).",".($acci+1)." -k".($starti+1)."n,".($starti+1)."n > $csv.sort";
	$cliuniq = "cut -d ';' -f2 $csv.sort | uniq > $csv.uniqregion";
}
#print "$clisort\n$cliuniq\n";
system ($clisort) == 0 || die "** Error on $clisort **\n";
system ($cliuniq) == 0 || die "** Error on $cliuniq **\n";

# GFF3 header
open(UNIQREGION, "$csv.uniqregion") || die "Error, can't open $csv.uniqregion file\n";
while (my $region = <UNIQREGION>) {
	chomp($region);
	if ($ref eq "acc") {
		print GFF "##sequence-region $region 1 ".$refsize{$region}."\n";
	}
	else { # ref = chr
		if ($region =~ /^[0-9]+$/ ) { # If chr number add Chr
        	print GFF "##sequence-region Chr$region 1 ".$refsize{$region}."\n";
        }
        elsif (lc($region) =~ /^mito/) {
        	print GFF "##sequence-region ChrM 1 ".$refsize{$region}."\n";
        }
        elsif (lc($region) =~ /^chloro/) {
            print GFF "##sequence-region ChrC 1 ".$refsize{$region}."\n";
        }
        else { # if scaffold or contig or other as /^chr/ || /^sc/ || /^contig/ ...
            print GFF "##sequence-region $region 1 ".$refsize{$region}."\n";
        }
    }
}
close(UNIQREGION);

# tRNA
open(CSV, "$csv.sort") || die "Error, can't open $csv.sort file\n";
# Traitement du fichier
while ($line = <CSV>) {
    chomp($line);
    my @tab = split (/;/,$line);

    my ($id, $name, $start, $stop,$strand);

    #Chr / source / seq_type
    $tab[$chri] =~ tr/\s//;

    # ID prefixe
    $tab[$chri] = "C" if ($tab[$chri] eq "chloroplast");
    $tab[$chri] = "M" if ($tab[$chri] eq "mitochondrion");
    if (($tab[$chri] !~ /^[Ss]caffold/) && ($tab[$chri] !~ /^[Cc]hr/) && ($tab[$chri] !~ /^[Cc]ontig/)){
        $id = "Chr".$tab[$chri]; # arath
    }
    else {
        $id = $tab[$chri]  ;
    }
    # Sequence reference
    if ($ref eq "chr") {
       $name=$id;

    }
    else {
        $tab[$acci] =~ tr/\s//;
        $name = $tab[$acci];
    }

    # chr positions
    if ($tab[$starti] < $tab[$stopi]) { # strand +
        $start = $tab[$starti];
        $stop = $tab[$stopi];
        $strand = "+";
    }
    else {# strand -
        $start = $tab[$stopi];
        $stop = $tab[$starti];
        $strand = "-";
    }
    #print "$id / $start / $stop / $strand \n";

    # ID / Type_name
    my $letter = "";
    $letter = "A" if ($tab[$typei] eq "Ala");
    $letter = "C" if ($tab[$typei] eq "Cys");
    $letter = "D" if ($tab[$typei] eq "Asp");
    $letter = "E" if ($tab[$typei] eq "Glu");
    $letter = "F" if ($tab[$typei] eq "Phe");
    $letter = "G" if ($tab[$typei] eq "Gly");
    $letter = "H" if ($tab[$typei] eq "His");
    $letter = "I" if ($tab[$typei] eq "Ile");
    $letter = "K" if ($tab[$typei] eq "Lys");
    $letter = "L" if ($tab[$typei] eq "Leu");
    $letter = "Me" if (($tab[$typei] =~ /^Met[_-]{0,1}[Ee]/) || ($tab[$typei] =~ /^eMet/));
    $letter = "Mi" if (($tab[$typei] =~ /^Met[_-]{0,1}[Ii]/) || ($tab[$typei] =~ /^iMet/));
    $letter = "M" if ($tab[$typei] eq "Met");
    $letter = "N" if ($tab[$typei] eq "Asn");
    $letter = "P" if ($tab[$typei] eq "Pro");
    $letter = "Q" if ($tab[$typei] eq "Gln");
    $letter = "R" if ($tab[$typei] eq "Arg");
    $letter = "S" if ($tab[$typei] eq "Ser");
    $letter = "U" if ($tab[$typei] eq "SeC");
    $letter = "T" if ($tab[$typei] eq "Thr");
    $letter = "V" if ($tab[$typei] eq "Val");
    $letter = "W" if ($tab[$typei] eq "Trp");
    $letter = "Y" if ($tab[$typei] eq "Tyr");
    $letter = "Sup" if ($tab[$typei] eq "Sup");
    $letter = "U" if ($tab[$typei] eq "Sec");
	$letter = "O" if ($tab[$typei] eq "Pyl");
    #print $letter."\n";


    if ($letter eq "") {
        print "Error : can't define letter for tRNA type : ".$tab[$typei]."\t".$tab[$aci]."\n" ;
    }
    else {
      	$typeNB{$letter} ++;
  		# print tDNA
      	print GFF "$name\tplantRNA\tgene\t".$start."\t".$stop."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.";Name=trn$letter-".$tab[$aci].";locus_type=pre_trna\n";
      	if ($ft eq "all"){ # print tRNA and exon / intros features
      		print GFF "$name\tplantRNA\ttRNA\t".$start."\t".$stop."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1;Parent=".$id.".trn$letter".$typeNB{$letter}.";Name=trn$letter-".$tab[$aci]."\n";
			if (($tab[$intronbi] != 0) && ($tab[$intronei] != 0)) { # intron
				if ($strand eq "+") {
                	print GFF "$name\tplantRNA\texon\t".$start."\t".($tab[$intronbi]-1)."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1.exon.1;Parent=".$id.".trn$letter".$typeNB{$letter}.".1\n";
                	print GFF "$name\tplantRNA\texon\t".($tab[$intronei]+1)."\t".$stop."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1.exon.2;Parent=".$id.".trn$letter".$typeNB{$letter}.".1\n";
            	}
				else {
                	print GFF "$name\tplantRNA\texon\t".$start."\t".($tab[$intronei]-1)."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1.exon.1;Parent=".$id.".trn$letter".$typeNB{$letter}.".1\n";
                	print GFF "$name\tplantRNA\texon\t".($tab[$intronbi]+1)."\t".$stop."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1.exon.2;Parent=".$id.".trn$letter".$typeNB{$letter}.".1\n";
            	}
			}
			else {
			    print GFF "$name\tplantRNA\texon\t".$start."\t".$stop."\t.\t$strand\t.\tID=".$id.".trn$letter".$typeNB{$letter}.".1.exon.1;Parent=".$id.".trn$letter".$typeNB{$letter}.".1\n";
			}
		}
    }
}

close GFF;
close CSV;

system("rm $csv.sort $csv.uniqregion") == 0 || die "** Error on rm $csv.sort **\n";

exit(0);
