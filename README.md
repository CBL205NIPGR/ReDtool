# ReDtool
ReDtool is a command line tool for in-silico digestion

## Requirement's
ReDtool is based on python3 and it doesn't require any additional dependencies, it is recommended to use the latest version of python.

## Running the ReDtool 
### 1. Download the repository
```
https://github.com/CBL205NIPGR/ReDtool
```
###  2. Help section of the tool
```
python ReDtool.py -h
```
User should provide three mandatory inputs to the tool, these are:
1. -i: input file, either in fasta format or txt format ex: sequence.fasta
2. -o: output file, it can be either tab-separated txt, xls, or ods file ex: output.txt
3. -r: recognition site mark with ^ at restriciton site

#### Other optional inputs based on the requirements 
1. -s: 'yes' For the sequence of digested fragments
2. -rc: 'yes' If digestion of reverse complement required
3. -c: 'yes' If input DNA is circular

### 3. Usage
```
python ReDtool.py -i Test_input.fasta -r ^gatc -o Test_output_tetra.txt -s yes
```
With restriction enzyme ^ symbol indicates the recognisiton site.
For hexa cutter
```
python ReDtool.py -i Test_input.fasta -r c^catgg -o Test_output_hexa.txt -s yes 
```
### Output files
Test_output_tetra.txt is a basic output file which includes recognition site psoition, corresponding frangment length and sequence, Test_output_tetra.txt is a hexa cutter restriciton site output. 

## Usage of the ReDtool.sh
It is used to atuomate the running of the 'ReDtool.py', in ReDtool.sh user needs to open and change input and output directories, recognisiton site and it is metioned in the script where the changes can be done based on inputs. Windows users cannot use this with either installing ubuntu terminal or using virtual box.
### 1. open the ReDtool.sh
```
nano ReDtool.sh
```
Make the required changes as indicated with comments in accordance with the requirements. save the file with 'ctrl+s' and exit with 'ctrl+x'.
### 2. Making the shell script executable
```
chmod +x ReDtool.sh
```
### 3. Running the shell script
```
./ReDtool.sh
```
