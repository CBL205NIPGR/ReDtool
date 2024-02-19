#!/usr/bin/env python3

# Name:     ReDtool.py
# Sign:     Shiva & Vikash
# Version:  1.0
#created:   Thu Feb 16 2024  

###########################

__author__	= " Shiva & Vikash"
__copyright__	= "GNU GENERAL PUBLIC LICENSE"
__email__	= " shivabrahmam455@gmail.com , vikash.yadav@nipgr.ac.in "


# import the requried libraries
import argparse
import time 

# storing the start time of program 
start_time = time.time() 


# function to extract dna sequence from the fasta file
def sequence(file):
    fin = open(f"{file}",'r')
    seq = ""
    # extracting the nucleotide sequence
    for lines in fin:
        if lines.startswith(">"):
          continue
        lines=lines.strip()
        seq+=lines
    return seq
    
# funtion to make reverse compliment of the sequence
def reverse_comp(seq):
    comp=seq.maketrans("ATGCatgc","TACGtacg")
    comp=seq.translate(comp)
    rev_comp = comp[::-1]
    return rev_comp

# funtion to make complementry sequence
def complementry(seq):    
    comp=seq.maketrans("ATGCatgc","TACGtacg")
    comp_trans=seq.translate(comp)
    return comp_trans

# creating the arguments for command line 
parser = argparse.ArgumentParser(
                    prog='ReDtool.py',
                    description='Identifying restriction sites in Linear and Circular DNA',
                    epilog='For liner digestion last position will be "0"')

parser.add_argument("-i","--inputfile", dest="inputfile", metavar="Input_file",type = str , required=True,
                     help = "Enter the input file" )

parser.add_argument("-r", "--res_site",dest="enzyme", metavar="restriction_site", type = str , required=True,
                     help = "Enter the recognition site")

parser.add_argument("-o", "--outputfile", dest="outputfile", metavar="Output_file", type = str, required=True,
                     help = "Enter the output file name, output can be either tab-delimited text, XML, and ods")

parser.add_argument("-s","--sequence",dest="sequence", metavar="Sequence", type = str, default= "no",
                     help = "Give if the sequence is requried")

parser.add_argument("-rc","--rev_comp",dest="compliment",metavar="reverse Complimentary", type=str,default="no",
                     help= "Enter if need to check for complimentary")

parser.add_argument("-c","--plasmid",dest="circular",metavar="circular", type = str, default= "no",
                     help = "Enter the yes to check for plasmid")

parser.add_argument("-v","--version", action="version", version="ReDtool v1.0")


# create the variablles to store the input values
# taking the input from users through command line
args = parser.parse_args()

file = args.inputfile
out = args.outputfile
seg = args.sequence
recog_site = args.enzyme
comp = args.compliment
plas = args.circular

# extracting the restriction enzyme
enzyme = recog_site.replace("^", "")

# extracting the sequence from file
seq = sequence(file)

# preparing inputs for reverse compliment
enz =complementry(enzyme)
recog_enz = reverse_comp(recog_site)

# function of restriciton enzyme identifier
def restriction_digest(seq,output,res,recog,seg,comp,plas):
    seq = seq
    enzyme = res.upper()
    out =output
    recog_site = recog.upper()
    seg = seg
    comp = comp
    plas = plas

    # identifying the recognisiton site
    recog_pos = recog_site.find("^")
    trim_site = recog_site[recog_pos+1:]

    # reading the input file
    fout = open(f"{out}",'w')

    # to check if the sequence has the restriction enzyme
    # if present execute this block of code 
    if enzyme in seq:

        pos=[]
        # identifying the positions of the restriction site
        index=seq.find(enzyme)
        while index != -1:
            pos.append(index)
            index = seq.find(enzyme, index+1)

        # add the lasting fragment position
        last_pos = 0 - recog_pos
        pos.append(last_pos)

        # finding the length of fragments
        frag_len=[]
        for x in range (len(pos)-1):
            frag=(pos[x+1] - pos[x])
            frag_len.append(frag)

        pos1=pos[0]
        frag_len.insert(0,pos1+recog_pos)

        frag_len.pop()

        # finding the length of last fragment and inserting
        last_frag = len(seq[pos[-2]+recog_pos:])
        frag_len.append(last_frag)

        # sequences of fragments
        first_frag_seq = seq[:pos[0]+recog_pos]
        frag_seqs = [first_frag_seq]

        for frag in range (len(pos)-1):
            seq_of_frags = seq[pos[frag]+recog_pos:frag_len[frag+1]+pos[frag]+recog_pos]
            frag_seqs.append(seq_of_frags)

        # finding the length of the sequences in complimentry
        if comp == "yes":
            len_seq = len(seq)
            pos_rev = []
            for positions in pos:
                rev_pos = len_seq - positions
                pos_rev.append(rev_pos)

            pos_rev[-1] = 0-recog_pos
        
        # indentifying the position, length, and sequence for the plasmid
        if plas == "yes":
            # positions of plasmid
            pos.pop()

            # length of the fragments
            f_len = frag_len.pop(0)
            l_len = frag_len.pop()
            last_fraglen_plas = f_len + l_len
            frag_len.append(last_fraglen_plas)
            
            # sequence of last fragment
            f_seq = frag_seqs.pop(0)
            l_seq = frag_seqs.pop()
            last_fragseq_plas = l_seq+f_seq
            frag_seqs.append(last_fragseq_plas)


        # counting the total numberr of restriction sites
        n_res_sites=len(pos)
        n_frags = len(frag_len)

        # inserting the header infromation lines in the file
        if plas == "yes":
            print(f"The total number of restriction sites identified for {enzyme}:",n_res_sites,"\n",file = fout)

        elif plas == "no":
            print(f"The total number of restriction sites identified for {enzyme}:",n_res_sites-1,"\n",file = fout)

        print(f"The total number of fragments identified: {n_frags}","\n",file=fout)
        print(f"The recognition site: {recog_site}","\n","\n",file=fout)


        # print the output into new file
        # prints based on the condition
        if comp == "no":
            if seg == "no":
                print("S.No",'\t',"position",'\t',"Fragment length",file=fout)
                for i in range(len(pos)):
                    print(i+1,'\t',pos[i]+recog_pos,'\t',frag_len[i],file=fout)

            elif seg == "yes":
                print("S.No",'\t',"position",'\t',"Fragment length",'\t',"Fragment sequence",file=fout)
                for i in range(len(pos)):
                    print(i+1,'\t',pos[i]+recog_pos,'\t',frag_len[i],'\t',frag_seqs[i],file=fout)
        
        elif comp == "yes":
            if seg == "no":
                print("S.No",'\t',"position",'\t',"Fragment length",file=fout)
                for i in range(len(pos)):
                    print(i+1,'\t',pos_rev[i]+recog_pos,'\t',frag_len[i],file=fout)
            
            elif seg == "yes":
                print("S.No",'\t',"position",'\t',"Fragment length",'\t',"Fragment sequence",file=fout)
                for i in range(len(pos)):
                    print(i+1,'\t',pos_rev[i]+recog_pos,'\t',frag_len[i],'\t',frag_seqs[i],file=fout)
            
        # closing the file once the it is completed
        fout.close()
        print('\n','\n','Successfully generated the file','\n','\n')
        

    # if it doesnt have restriciton enzyme print this 
    else:
        print(f"\n Sequence dosen't have restriction site {enzyme} \n")

# to start the digestion 
if comp == "yes":
    comp_seq = reverse_comp(seq)
    restriction_digest(comp_seq,out,enzyme,recog_site,seg,comp,plas)

elif comp == "no":
    restriction_digest(seq,out,enzyme,recog_site,seg,comp,plas)

# storing the end time of program
end_time = time.time() 

# Calculating the runtime in seconds
runtime_seconds = end_time - start_time

# Converting seconds to hours, minutes, and seconds
hours = int(runtime_seconds // 3600)
minutes = int((runtime_seconds% 3600) // 60)
seconds = int(runtime_seconds% 60)

# total run time of program
print(f"Run time: {hours:02}hr:{minutes:02}min:{seconds:02}sec \n")