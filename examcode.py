#!/usr/bin/env python3
# Programme for BLAST analysis
# Written by B203499, 10 Dec 2021
# Import the required modules
import sys, subprocess # for system interaction e.g. command line arguments
import os, shutil # for working with filesystems
print ("Imported sys, subprocess, os, shutil")
import numpy as np # for scientific computing
print ("Imported numpy as np")
import re # for regular expressions
print ("Imported re")
# for BioPython
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 's2248602@ed.ac.uk' # otherwise gives bug
print ("Imported BioPython modules: Seq, SeqIO, Entrez")
import pandas as pd
print ("Imported pandas")

####################### User info setup #####################
original_stdout = sys.stdout

# to evaluate the user input making sure either "y" OR "n" is correctly entered, it loops until returns boolean:
def eval(user_bool):
	# user_bool = input("\n>Please input y/n:\t").lower().strip()
	if user_bool == "y":
		return True
	if user_bool == "n":
		return False
	else:
		print("[Warning] CANNOT take input other than y OR n, please re-enter...")
		user_bool = input("\n>Please RE-ENTER [y/n] [default: y]:\t") or "y"
		return eval(user_bool)

Entrez_email = input("\n> Please enter your email address for NCBI access: \t") or "s2248602@ed.ac.uk"
Entrez_api_input = input("\n> Do you wish to use your NCBI API KEY? [y/n]:\t") or "y"
if eval(Entrez_api_input) == True:
	Entrez_default = subprocess.check_output("echo $NCBI_API_KEY", shell=True).rstrip().decode()
	Entrez_api_key = input("Please enter your NCBI API KEY: \t") or Entrez_default

#################################################################

####################### BLAST search #####################

while True:
	blastn_input = input("> Do you wish to perform BLASTN for nucleotide against nucleotide search? [y/n]:\t") or "n"
	if eval(blastn_input) == True:
		n_db = input("> Please enter the database for BLASTN") or "nt"
		n_query = input("> Please enter the query term for BLASTN") or "Mammalia COX1 complete"
		blastn_com = "blastn -db {0} -query {1} -out blastn.out -outfmt 6 -max_target_seqs 50 -num_threads 64".format(n_db, n_query)
		subprocess.call(blastn_com, shell=True)
		filename = "blastn.out"
		break
	elif:
		blastx_input = input("> Do you wish to perform BLASTX for translated nucleotide against protein search? [y/n]:\t") or "n"
		if eval(blastx_input) == True:
			x_db = input("> Please enter the database for BLASTX") or "nr"
			x_query = input("> Please enter the query term for BLASTX") or "Mammalia COX1 complete"
			blastx_com = "blastx -db {0} -query {1} -out blastx.out -outfmt 6 -max_target_seqs 50 -num_threads 64".format(x_db, x_query)
			subprocess.call(blastx_com, shell=True)
			filename = "blastx.out"
			break
	elif:
		tblastn_input = input("> Do you wish to perform TBLASTN for protein against nucleotide search? [y/n]:\t") or "n"
		if eval(tblastn_input) == True:
			tn_db = input("> Please enter the database for TBLASTN") or "nt"
			tn_query = input("> Please enter the query term for TBLASTN") or "Mammalia COX1 complete"
			tblastn_com = "tblastn -db {0} -query {1} -out tblastn.out -outfmt 6 -max_target_seqs 50 -num_threads 64".format(tn_db, tn_query)
			subprocess.call(tblastn_com, shell=True)
			filename = "tblastn.out"
			break
	else:
		blastp_input = input("> Do you wish to perform BLASTN for protein against protein search? [y/n]:\t") or "n"
		if eval(blastp_input) == True:
			p_db = input("> Please enter the database for BLASTP") or "nr"
			p_query = input("> Please enter the query term for BLASTP") or "Mammalia COX1 complete"
			blastp_com = "blastp -db {0} -query {1} -out blastp.out -outfmt 6 -max_target_seqs 50 -num_threads 64".format(p_db, p_query)
			subprocess.call(blastp_com, shell=True)
			filename = "blastp.out"
			break

## should also allow for fasta output for downstream analysis.

####################### BLAST search done #####################

# use PANDA to retrieve the list of accession number given by the blastp result :
dir_bpout = os.getcwd() + "/" + filename
df = pd.read_csv(dir_bpout, sep='\t', names=["Query acc.", "Subject acc.", "Percent ID", "Length", "Mismatch", "Gapopen", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore"])
acc_list = list(df["Subject acc."])
num_seq = len(acc_list)
print("The result file from previous BLAST <{0}> has {1} sequences.".format(filename, num_seq))
with open(filename+"_acclist.txt", 'w') as fi:
	sys.stdout = fi
	print(*acc_list, sep='\n')
	sys.stdout = original_stdout

## should also give option to store fasta files.

####################### Get SeqRecord objects for analysis #####################

## 1. Read-in from local files:
scratch_input = input("\n> Do you wish to read in your Seq for BLAST analysis from local files? [y/n]:\t") or "n"
# use variables to keep track of the user's choice:
scr_var = eval(scratch_input)
fasta_var = False
while scr_var == True:
	scr_fasta_input = input("\n> Do you wish to read in your Seq from local FASTA files? [y/n]:\t") or "n"
	if eval(scr_fasta_input) == True:
		fasta_var = True
		fastafile_default = "/localdisk/data/BPSM/Lecture19/sequences.fasta"
		fastafile_input = input("\n> Please enter the local path of the FASTA file :\t") or fastafile_default
		# make sure the input path exist:
		while os.path.isfile(fastafile_input) == False:
			print("[Warning] This given path does NOT seem to exist, please re-enter...")
			fastafile_input = input("\n> Please re-enter the path of the FASTA file:\t") or fastafile_default
		# read-in the FASTA sequence and generate records:
		records = SeqIO.parse( fastafile_input, "fasta" )
		ids_dict = SeqIO.to_dict(SeqIO.parse(fastafile_input, "fasta"))
		for seq_record in ids_dict.values():
			print(seq_record)
		break
	else:
		scr_genbank_input = input("\n> Do you wish to read in your Seq from local GenBank files? [y/n]:\t") or "n"
		if eval(scr_genbank_input) == True:
			genbank_default = "/localdisk/data/BPSM/Lecture19/sequence.gb"
			genbank_input = input("\n> Please enter the local path of the GenBank file :\t") or genbank_default
			# make sure the input path exist:
			while os.path.isfile(genbank_input) == False:
				print("[Warning] This given path does NOT seem to exist, please re-enter...")
				fastafile_input = input("\n> Please re-enter the path of the GenBank file:\t") or genbank_default
			# read-in the GenBank file and generate records:
			records = SeqIO.parse(genbank_input, "genbank")
			# prompt sequence features to the user:
			for seq_record in records:
				for feature in seq_record.features:
					print(feature.type)
					print(feature.location)
			break

## 2. Search with user-defined values for SeqRecord:
def efetch_id():
	db_input = ("\n> Please enter the name of the database to search with:\t") or "protein"
	id_input = ("\n> Please enter the ID of the record:\t") or "NP_008368.2"
	genbank = Entrez.efetch(db=db_input, id=id_input, rettype="gb", retmode="text")
	record = genbank.read()
	return result

def esearch_term():
	db_input = ("\n> Please enter the name of the database to search with:\t") or "protein"
	term_input = ("\n> Please enter the term to search with:\t") or "Mammalia COX1 complete"
	retmax_input = ("\n> Please enter the return maximum of results:\t") or "10"
	result = Entrez.read(Entrez.esearch(db_input, term_input, retmax_input))
	return result

while scr_var == False:
	# Retrieve a user-defined SeqRecord ID:
	esear_ID_input = input("\n> Do you wish to retrieve SeqRecord with a specific ID? [y/n]:\t") or "n"
	if eval(esear_ID_input) == True:
		record = efetch_id()
		break
	else: 
		# Search GenBank for Seqrecord with a user-defined term:
		esear_term_input = input("\n> Do you wish to read in your Seq from local GenBank files? [y/n]:\t") or "n"
		if eval(esear_term_input) == True:
			record = esearch_term()
			break

if scr_var == False:
	for seq_record in records:
		for feature in seq_record.features:
			print(feature.type)
			print(feature.location)

if fasta_var = True:
	id_list = []
	for seq_record in records:
		id_list.append(seq_record.id)
		records = []
	for i in id_list:
		genbank = Entrez.efetch(db="protein", id=i, rettype="gb", retmode="text")
		record = genbank.read()
		records.append(record)
	
####################### SeqRecord Genarated and Sequence features Prompted #####################

## use the sequence fetched from ID:
clustalo_msf = "clustalo -i {0} -o {1}.msf --threads=64 --outfmt=msf".format(fastafile_input, filename)
subprocess.call(clustalo_msf, shell=True)
print("Completed generating output files [{0}.msf].".format(filename))

