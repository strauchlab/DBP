import sys
import operator
import os
import shutil,  subprocess
import re
#from itertools import permutations
#from itertools import combinations
from itertools import *

########## input ###############

starting_pdb = "ress.pdb"
starting_blueprint = open('starting.blue', 'r').readlines()

######## #what do we want to sample? ##############

topology = "H14p,L1,H14p,L1,H14p" #"E7,L1,H14p,L2,M,E1"
topology2 = "" #"L3,H14p,L4,E7"

L1 = ['BB','GBB','GB','GABB']

########### blocks ################## 
E5p = [5,6]
E5 = [5]
E7 = [7]
E7p = [7,8]
E1 = [0,1 ]
E2 = [2]
E3 = [3]
E8 = [8]
H10p = [10,11]
H14p = [14,15,16]
H18p = [17,18,19]
H1418 = [14,15]

loop_lengths = [3,4]
loop2_lengths = [1,2]
M = [len(starting_blueprint)]
motif_start_tracker = 1


######### command line options   ######################
def command_line(cdir,pdb,bp,res_start,res_stop):
    residues = str(res_start) 
    for i in range(res_start+1, res_stop+1):
	residues += ','+str(i)
    command = "cd " + cdir +"\n"
    command += "cp /work/strauch/denovoproteins/HHH/new/ress.pdb . \n"

    command += "mv " + pdb +" " + bp + ".pdb \n" 
    command += "cp /work/strauch/denovoproteins/HHH/new/desHHH.xml . \n"
    command += " /work/strauch/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -nstruct 1000 -jd2:ntrials 1 -parser:protocol"
    command += " desHHH.xml  -database /work/strauch/Rosetta/main/database  -s " + bp + ".pdb " + " -mute all -parser:script_vars blue=" + bp +".blue "
    command += " -holes::dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc "
    #command += " -parser:script_vars blue2=" + bp2 +".blue "
    #command += " -parser:script_vars begin=" + str(res_start) + "  -parser:script_vars residues=" + residues 
    #command += " -parser:script_vars bres=5 -parser:script_vars sres=6 "
    #command += " -parser:script_vars end=" + str(res_stop)
    
    command += " -no_his_his_pairE -nblist_autoupdate true -chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm"
    return command

#########   sheet paring  ##################
#header = "SSPAIR 1-3.A.0;1-4.A.0;2-3.A.0" #"SSPAIR 1-2.A.0;2-3.A.0"
headers = []

for i in range(0,4):
	for j in range(-2,2):
		headers.append('SSPAIR 1-3.A.' + str(i) +';2-3.A.' + str(j))
print "genearte headers: " ,headers

######### output directory ##############

#make a new folder for the organized decoys
outdirectory = 'sample_topos_two_a'# + topology
if not os.path.exists(outdirectory): os.makedirs(outdirectory)


########### stuff ###########
motif =""
for line in starting_blueprint:
    print line.rstrip()
    motif += line#.rstrip()
print "length " ,len(starting_blueprint)

def make_loopGB(outfile, slen):
    #for i in range(slen):
    outfile.write("0 V LG R\n")
    outfile.write("0 V LB R\n")

def make_sheet(outfile, slen):
    for i in range(slen):
        outfile.write("0 V EB R\n")

def make_helix(outfile, slen):
    for i in range(slen):
        outfile.write("0 V HA R\n")


def place_abegos(outfile, abegos):
    for i in abegos:
	outfile.write("0 V L" + i + " R\n")


curr_dir =  os.getcwd()
print curr_dir

def condor_com(cdir):
    com = "Executable  = "+ cdir + '/run.sh\n'
    com += "transfer_executable = false\n"
    com += "universe    = vanilla\n"
    com += "Error       = "+ cdir + '/stderr.log\n'
    com += "Output      = "+ cdir + '/stdout.log\n'
    com += "Log         = " + cdir + '/condor.log\n'
    com += 'queue 1\n'
    return com

dirlist = open('subfolders_topos.list', 'w')

tops = topology.split(',')
print "-----array with the split string -------" , tops

permutations = []
name = "p1"
elements = []

for item in tops:
	if 'E' in item:
		name += '-%sE'
		elements.append(eval(item))
	elif 'H' in item:
		name += '-%sH'
		elements.append(eval(item))

	elif 'L' in item:
		name += '-%sL'
		elements.append(eval(item))
	
	elif 'M' in item:
		name += '-%sM'
		elements.append(eval(item))

	else:
		print "something wrong in teh topology definition"

print elements
print "----------------- name" , name

permutations =  product(*elements) 
permuts = []
for i in permutations:
	permuts.append(i)
counter = 0
#counter_headers = 

#### need to create a method at some point instead of just repeating ....
tops2 = topology2.split(',')
print "-------------" , tops2

permutations2 = []
name2 = ""
elements2 = []

for item in tops2:
        if 'E' in item:
                name2 += '-%sE'
                elements2.append(eval(item))
        elif 'H' in item:
                name2 += '-%sH'
                elements2.append(eval(item))

        elif 'L' in item:
                name2 += '-%sL'
                elements2.append(eval(item))

        elif 'M' in item:
                name2 += '-%sM'
                elements2.append(eval(item))

        else:
                print "something wrong in teh topology definition"

print elements2
print "---- name 2  ------" , name2

res = 0
motif_adjustment = 5
res = motif_adjustment

permutations2 =  product(*elements2)
permuts2 = []
for i in permutations2:
        permuts2.append(i)

print "permuts2: --- " , permuts2


for topo in permuts:
    res = motif_adjustment
    counter += 1
    # first blue print
    bluename = name  % (topo)
    #bluename += '_' + replaced	

    if not os.path.exists(outdirectory + '/' + bluename): os.makedirs(outdirectory + '/' + bluename)
    print "bluename:", bluename

    name_items = name.split('-')
 
    dirlist.write(bluename + '\n')
    bpfile = open(outdirectory + '/' + bluename + '/' + bluename + '.blue', 'w')
    #bpfile.write(header + '\n') 
    
    print >> bpfile, "1 L LX ."
    firstsection = True #to include the silly 2nd residue from starting pdb

    for i in range(len(topo)):

	#ss = int(topo[i])

	if firstsection is True:
	   #ss  = ss -1 dd
	   firstsection = False
		#todo : change to generic
           print >> bpfile, "2 V HA R"

           make_helix(bpfile, int(topo[i]) -1 )
           res += int(topo[i]) + 2 # counting the very first residue here too


	elif 'H' in name_items[i+1] :
	   make_helix(bpfile, int(topo[i]) )
	   res += int(topo[i])
	elif 'E' in name_items[i+1]:
	   make_sheet(bpfile, int(topo[i]) )
	   res += int(topo[i])
	elif 'M' in name_items[i+1]:
	   bpfile.write(motif)
	   
	elif 'L' in name_items[i+1]:
	   place_abegos(bpfile, topo[i])
	   res += len(topo[i])

    #silly thing but apparently the last residues always seems to be a loop these days
    print >> bpfile, "0 V LX R"
    res += 1

    bpfile.close()
    print "residues " , res
    subprocess.Popen( "cp " + starting_pdb + " " + outdirectory + '/' + bluename ,shell=True)

    outcommand = open(outdirectory + '/' + bluename + '/run.sh', 'w')
    outcondor = open(outdirectory + '/' + bluename + '/job.cdr' , 'w')

    subprocess.Popen( "chmod +x " + outdirectory + '/' + bluename + '/run.sh' ,shell=True)

    outcommand.write("#!/bin/bash\n")
	## todo: counter for where the motif starts to set right what should or should not be rebuild ...
    print >> outcommand, command_line(curr_dir + '/' + outdirectory + '/' + bluename, starting_pdb , bluename,res+2,res+10 )
    print >> outcondor, condor_com(curr_dir + '/' + outdirectory + '/' + bluename )
    outcommand.close()

dirlist.close()

print "sampled topologies: " , counter
