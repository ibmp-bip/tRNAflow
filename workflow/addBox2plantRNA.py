#! /usr/bin/env python

##########################################################################################
# Name:		addBox2plantRNA.py
# Author:	Valerie Cognat valerie.cognat@ibmp-cnrs.unistra.fr
# Copyrigth :	IBMP - CNRS
# Description:	complete information in nucleus tRNA file : box A/B and pattern CAA/polyT
#
# Created : 2021, July the 21th
# Version 1.0 - 27 july 2021 - adapt plant2tRNA.py for csv
# Version 1.6 - 3 June 2021 - Remove CCA & polyT pattern in organelles
##########################################################################################


import os
import pandas as pd
import sys
import re
import argparse
import textwrap
from Bio import SeqUtils

## FUNCTIONS ##

# function to extract CAA motif
def get_caa(myupstream):
	# if nan value, specific property - return empty pattern
	if(myupstream != myupstream):
		return ""
	#print(myupstream)
	#match=re.search(r'CAA[a-z]{0,10}CAA[a-z]{0,10}CAA', myupstream) ## returns only the leftmost pattern - need the rightmost for CAA
	caamatch=re.findall(r'CAA', myupstream)
	if(len(caamatch)>0): # At least a CAA motif
		caaend = myupstream.rfind('CAA') + 3 # last CAA pos in the upstream # the biggest number
		#print("CAA: ", caaend)

		# Search differents patterns with 2, 3 or 4 CAA, keep the longest if last CAA in pattern is ended as the same position at the last CAA in upstream
		match=re.findall(r'CAA[a-z]{0,10}CAA[a-z]{0,10}CAA[a-z]{0,10}CAA', myupstream)
		if(len(match)>0): ## motif exists
			mypattern = match[-1]#rightmost pattern
			pos = myupstream.rfind(mypattern) + len(mypattern) # pos of the end of the pattern
			#print("BIG: ",mypattern, pos)
			if (pos >= caaend): # if the same as the last CAA
				return mypattern # return the largest pattern

		match2=re.findall(r'CAA[a-z]{0,10}CAA[a-z]{0,10}CAA', myupstream)
		if(len(match2)>0): ## motif exists
			mypattern = match2[-1]#rightmost pattern
			pos = myupstream.rfind(mypattern) + len(mypattern) # pos of the end of the pattern
			#print("MED: ",myupstream, mypattern, pos)
			if (pos >= caaend): # if the same as the last CAA
				return mypattern # return the largest pattern

		match3=re.findall(r'CAA[a-z]{0,10}CAA', myupstream)
		if(len(match3)>0):
			mypattern = match3[-1]
			pos = myupstream.rfind(mypattern) + len(mypattern) # pos of the end of the pattern
			#print("SMALL: ",myupstream, mypattern, pos)
			if (pos >= caaend): # if the same as the last CAA
				return mypattern # return the largest pattern

		return caamatch[-1] # return the rightmost CAA pattern if no longer motif found or longer motif before last CAA
	else: # pas de CAA
		return ""

# function to extract CAA motif position
def get_CAAposition(motif,utr):
	if (motif != ""):
		return(-len(utr) + utr.rfind("CAA")) # since version 1.2
		#return(-len(utr) + utr.rfind(motif)) # in version 1.0 & 1.1
	else:
		return ""

# funtion to extract polyT motif
def get_polyT(mydownstream):
	#print(mydownstream)
	# if nan value, specific property - return empty pattern
	if(mydownstream != mydownstream):
		return ""
	match=re.search(r'T{4,25}', mydownstream) ## returns the leftmost pattern
	if(match is None):
		return ""
	else:
		return match.group()

# function to extract polyT motif position
def get_polyTposition(motif,utr):
	if (motif != ""):
		return(utr.find(motif) + 1)
	else:
		return ""

#---
# funtion to extract boxA motif, its location and its length
def get_boxA(mytrna):
	#consensus = "TRGYNNANNNG" / "TRNYNNANNNG" / "TNNYNNANNNG"
	#print(mytrna)
	# if nan value, specific property - return empty pattern
	if(mytrna != mytrna):
		return ["","","",""]
	match=re.search(r'T.{2}[CT].{2}A.{2,3}GG', mytrna.upper())
	#print (match)
	if(match is None): # boxA not found
		return ["","","",""]
	else:
		if ((match.start() > 5) & (match.start() < 11)):
			return [match.group(0),match.start()+1,len(match.group(0)),"predicted"]
		else:
			return ["","","",""]



# funtion to extract polyT motif
def get_boxB(mytrna):
	#print(mytrna)
	#consensus = "GWTCRANNC"
	consensus = "NGHTYNRNNYN"
	length = len(consensus)
	# if nan value, specific property - return empty pattern
	if(mytrna != mytrna):
		return ["","","",""]
	match=SeqUtils.nt_search(mytrna.upper(),consensus)
	#print (match)
	if(len(match) == 1):  # consensus not found
		return ["","","",""]
	else:
		for i in range(1, len(match)):
			if (match[i] > 48): # return the first motif around position 50.
				return [mytrna[match[i]:(match[i]+11)], match[i]+1, length, "predicted"]
		# if nothing return :
		return ["","","",""]


# check args
def main(args):
	print ("Parameters: ")
	print ("- tRNA file: ", args.input)
	print ("- new tRNA file: ", args.output)
	print ("- box: ", args.nobox)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog = "addBox2plantRNA.py", formatter_class=argparse.RawDescriptionHelpFormatter, description = textwrap.dedent('''\
		Script to extract CAA and polyT motifs in UTRs and predict box A/B in tRNA file.
		  - CAA pattern will be the rightmost one in the upstream region.
		       [A CAA pattern is defined by CAA separated by less than 10 nucleotides.]
		  - polyT pattern will be the leftmost stretch of more 4 T in the downstream region.
		  - box A pattern looks like TNNYNNANN(N)GG near position 8 in tRNA
		  - box B pattern looks like  NGHTYNRNNYN at the end of the tRNA (> position 48)
		  '''))
	parser.add_argument('-i', '--input', type=str, help='csv tRNA file', required=True)
	parser.add_argument('-o', '--output', type=str, help='csv output tRNA file', required=True)
	parser.add_argument('--nobox', action='store_true', help='to add empty box and pattern columns in the output file - for organellar genomes')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()

	main(args)

	# read csv file and convert in dataframe
	df_tRNA = pd.read_csv(args.input, sep = ';')

	#Add CCA column
	df_tRNA['CCA'] = "CCA"

	if (args.nobox):
		df_tRNA[['boxA_seq','boxA_start','boxA_length','boxA_comment']] = ["","","",""]
		df_tRNA[['boxB_seq','boxB_start','boxB_length','boxB_comment']] = ["","","",""]
		df_tRNA[['CAA_seq','CAA_start','CAA_comment']] =  ["","",""]
		df_tRNA[['polyT_seq','polyT_start','polyT_comment']] =  ["","",""]

	else:
		# add boxA information
		df_tRNA[['boxA_seq','boxA_start','boxA_length','boxA_comment']] = [get_boxA(trna) for trna in df_tRNA['Sequence']]
		# add boxB information
		df_tRNA[['boxB_seq','boxB_start','boxB_length','boxB_comment']] = [get_boxB(trna) for trna in df_tRNA['Sequence']]
		# add CAA motif
		df_tRNA['CAA_seq'] =  [get_caa(up) for up in df_tRNA['upstream']]
		# add CAA motif position
		df_tRNA['CAA_start'] = df_tRNA.apply(lambda x: get_CAAposition(x.CAA_seq, x.upstream), axis=1)
		# add CAA comment
		df_tRNA['CAA_comment'] = ""
		# add polyT motif
		df_tRNA['polyT_seq'] =  [get_polyT(down) for down in df_tRNA['downstream']]
		# add polyT motif position
		df_tRNA['polyT_start'] = df_tRNA.apply(lambda x: get_polyTposition(x.polyT_seq, x.downstream), axis=1)
		# add polyT comment
		df_tRNA['polyT_comment'] = ""
		# add Import columns
		df_tRNA[['import_type','import_Pred/Exp','import_comment']] = ["","",""]

	#Add tRNA_comment column
	df_tRNA['tRNA_comment'] = ""

	# write output file
	df_tRNA.to_csv(args.output, sep = ';', index = False)
