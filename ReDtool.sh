#!/bin/bash
#used for running multiple files at a time

#give the location of input and output directories
#use "pwd" to know the directory in terminal or cmd prompt
in_dir="/mnt/e/ASUS/Documents/python_files/homo"       #change according to input directory 
out_dir="/mnt/e/ASUS/Documents/python_files/homo"     #change according to output directory
res_site="a^agctt"				  #change the resriction site
sequence="no"					  #change to yes if sequence is requried, else no
reverse_compliment="no"				  #change to yes to perform reverse compliment, else no
circular="no"					  #change to yes to do circular digestion, else no for linear


#running the loop to run the ReDtool
for f_file in "$in_dir"/*.fa; do		  #based on the extension of file change the extension to fasta/txt/fa/fna
	echo "Input file:"$f_file
	#ectracting the basename of the file
	filename=$(basename "$f_file" .fa)	  #based on the extension of file change the extension to fasta/txt/fa/fna

	#output file name
	o_file="$out_dir/$filename.txt"		  #based on the requirement of the output, extension can be changed to txt/xls/ods
	echo "Output file:"$o_file

	#Running the ReDtool

	python3 ReDtool.py -i "$f_file" -r "$res_site" -o "$o_file" -rc "$reverse_compliment" -c "$circular" -s "$sequence"
done
