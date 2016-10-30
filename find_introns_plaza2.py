#!/usr/bin/python
#usage: python find_introns_plaza.py scaffoldfile.fas annotationfile.csv cdsfastasequences.fas
#created: Oct 22 of 2016
#Use CDS fasta sequences and Annotation(csv) file to get correspondent Genomic sequences
#in this version only passes to memory - dictionary - the genomic scaffolds belonging to respective CDS sequences

import sys
import subprocess

filescaf = sys.argv[1]	#Scaffold file in fasta format containing genomic sequences
fileanno = sys.argv[2]	#Annotation file in csv format
filelist = sys.argv[3]	#CDS Fasta file containing sequences to obtain the respective genomic sequences

#open files
try:
    fgen = open(filescaf)
    flist = open(filelist)
except IOError:
    print ("Scaffold and FASTA files doesn't exist!")

#open annotation file and creates an set for it
set_anno = set(line.strip() for line in open (fileanno, 'r'))



#Create a dictionary with CDS Fasta sequences
fastalist={}

for line in flist:
    line = line.rstrip()
    if line[0] == '>':
        words=line.split() 
        name=words[0][1:]
        fastalist[name]=''
    else:
        fastalist[name] = fastalist[name] + line



#Create a dictionary with annotation correspondent to all the CDS sequences in previous dictionary
cds_anno = dict()           #value in this dictionary will be a tupple
listscaffolds = []			#criates an list of scaffolds corresponding to cds's

for i in fastalist.keys():  #read keys from CDS fasta dictionary
	for line in set_anno:   #read lines in the annotation set
		if i in line:
			line = line.split('";"')
			start_cds = line[4]	#start coordinate
			stop_cds = line[5]	#stop coordinate
			scaffold = line[9]	#scaffold name
			cds_anno[i] = start_cds,stop_cds,scaffold	#add values to dictionary where the key is the CDS id
			listscaffolds += scaffold.splitlines()


#Create dictionary for all the Scaffolds in the genomic fasta file 
seqs={}
parsing = 0
for i in fgen:
	i = i.rstrip()
	
	if ">" in i:
		scafid = i.split()[0][1:]
	
	if parsing > 0 and ">" not in i:          #Second 2: if parsing is greater than 0, than this is the first line of fasta belonging to a scaffold selected in the first step. When it finds a line with ">", then parsing will be set to 0 and everything starts begin
		seqs[scafid] = seqs[scafid] + i
	else:
		parsing = 0

	if ">" in i and scafid in listscaffolds:   #First 1: if id scaffold in my list of scaffolds, then parsing will be set 1
		seqs[scafid] = ''                      #create a dictionary with the scaffold header
		parsing += 1

	
#Compare all three files based in dictionary cds_anno containing the coordinates of CDS correspondent genomic sequence in scaffolds
for i in cds_anno.items():

	cdsid = i[0]
	start_cds = list(i[1])[0]
	stop_cds = list(i[1])[1]
	scaffold = list(i[1])[2]
	
	if scaffold in seqs.keys():
		if cdsid in fastalist.keys():
			
			genomic = str(seqs[scaffold])[int(start_cds)+1:int(stop_cds)+1]
			output = open(cdsid+"_"+scaffold+".fasta", "w")       #open output for writing
			output.write(">"+scaffold+"\n"+str(genomic)+"\n"+">"+cdsid+"\n"+fastalist[cdsid]);   #write FASTA output containing the CDS and Genomic sequences
			output.close()
			doalignment = subprocess.call("muscle \-in"+" "+str(cdsid+"_"+scaffold+".fasta")+" "+"\-out"+" "+str(cdsid+"_"+scaffold+".aln")  ,shell = True)  #call shell subprocess using Muscle to align the CDS sequence to its correspondent Genomic sequence 
