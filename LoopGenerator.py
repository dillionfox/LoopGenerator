#!/usr/bin/env python
"""
REQUIREMENTS: Numpy, BioPython, CCD, and LoopTools.
This code was designed to take a protein structure, make a copy of it,
shift it, and then connect the N-terminus of the original protein
to the C-terminus of the copied protein with a LOOP. The loop
can be generated 2 different ways, but the first method (denoted
OPTION 1) is the highly preferred method. There are several other ways
of doing this, but OPTION 1 works well enough that I didn't see a
need to pursue any of them.

Dillion Fox, 3/2018

"""

import numpy as np; from numpy import linalg as npl
from numpy.random import random
import Bio.PDB
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
import CCD
import LoopTools

print "~~~ TODO ~~~"
print "1. I think the next thing to do is figure out how to port this"
print "code into AUTOMACS. From there I can figure out how to write"
print "pdb files containing the combined structures, and I can also"
print "have it throw away structures that have steric clashes."
print ""
print "2. make it work for coarse-grained models"
print ""

#########################################################
### load USER DEFINED variables into global namespace ###
#########################################################

global RMSD_threshold; RMSD_threshold = 0.2			# cut off for loop closing
global max_iterations; max_iterations = 500			# number of iterations for CCD algorithm
global max_success;    max_success = 10000			# maximum number of structures you want to make
#global dl;             dl = 2 					# variability in chain length. i.e. 15 +/- 1
global DIST_FACTOR;    DIST_FACTOR = 0.2			# N-to-C terminal distance variability for searching ArchDB for structures
								# pdbname is your input structure
global pdbname ;       pdbname = "/home/dillion/data/reflectin/add_linker/structures/small_mem_protein-inserted_frame-738.pdb"
global protein_shift;  protein_shift  = np.array([30,0,0])	# how much you want to shift your input structure


def mine_loop_data(loop_length,CC_DIST,LEN_TOL):
	"""
	This function calls the LoopMiner class and extracts dihedral data from
	ArchDB and writes it to a file.

	"""

	OUTFILE = "loop_dihedrals.out"				# write dihedrals to output file
	ArchDB = "ArchDB/archdb.mini.tab"			# ArchDB database
	DIST_TOL = np.floor(CC_DIST*DIST_FACTOR)		# tolerance for identifying geometrically similar loops
	numResLoop = range(int(loop_length-LEN_TOL),int(loop_length+LEN_TOL+1))	# desired number of residues in loop
	print "looking for loops with", loop_length, "+/-", LEN_TOL, ", and span a distance of", CC_DIST, "+/-", DIST_TOL
	miner = LoopTools.LoopMiner(numResLoop,CC_DIST,DIST_TOL,OUTFILE)
	[dihedral_list,native_structures] = miner.mineArchDB(ArchDB)
	return [dihedral_list,native_structures]

def extract_coors(structure):
	"""
	This function extracts coordinates from a BioPython structure and returns
	a numpy array

	"""

	import Bio.PDB
	v = []
	for model in structure:
		for chain in model:
			for residue in chain:
				for atom in residue:
					if isinstance(atom.get_coord(), (list,)):
						v.append(np.array(atom.get_coord()))
					else:
						l = []
						l.append(atom.get_coord()[0])
						l.append(atom.get_coord()[1])
						l.append(atom.get_coord()[2])
						v.append(np.array(l))
	return v

def compute_CC_DIST(first_res,target_coors):
	"""
	This function computes the C-terminal to N-terminal distance that
	the linker will have to connect. This data is readily accessible in 
	ArchDB so you can pre-select loops that have similar geometry to yours

	"""

	Bio.PDB.PDBIO().set_structure(first_res)
	l = []
	base_coors = []
	for model in first_res:
		for chain in model:
			for res in chain:
				for atom in res:
					l = []
					l.append(atom.get_coord()[0])
					l.append(atom.get_coord()[1])
					l.append(atom.get_coord()[2])
					base_coors.append(np.array(l))
	CC_DIST = base_coors[-1] - target_coors[0]
	return npl.norm(CC_DIST)

def run_CCD(chain,target_coors,success=0,n=1):
	"""
	This is the main function that calls the Cyclic Coordinate Descent 
	class functions

	"""

	init = CCD.CCD(chain, target_coors) 
	check, RMSD, it, V = init.run(n,max_it=max_iterations,threshold=RMSD_threshold)
	if check == 'success': 
		success += 1
	print("{0}\t{1}\t\t{2}\t{3}\t{4}".format(n, str(round(success*100/float(n+1),2))+'%', check, str(RMSD), it))
	return [success, V]

def make_first_res():
	"""
	This function loads your structure with BioPython and extracts the 
	LAST residue. This residue will be the FIRST residue of the loop, but
	will later be deleted. I made it this way to ensure a smooth transition
	from native structure to loop. It was the easiest way that I could think 
	of to check the dihedrals of all of the loop residues.

	"""

	loopmaker = LoopTools.LoopMaker()
	[base_res, target_coors] = loopmaker.define_initial_cond(pdbname,protein_shift)
	first_res = loopmaker.anchor_is_start(base_res)
	return [first_res,target_coors]

import mda
def run_loop(n):
	di = dihedral_list[n]
	[first_res,target_coors] = make_first_res()
	loop = loopmaker.finish_grafted_chain(first_res,di)
	[success,V] = run_CCD(extract_coors(loop),target_coors,success,n)
	if V != 0:
		mda.write_merged(pdbname, "linker"+str(n)+".pdb",n)
		print "[energy, identification number]", V, n 

################
### OPTION 1 ###
################
def run_graft(loop_length):
	"""
	This function mines a Loop Database to find naturally occurring loop
	structures with criteria that matches the loop you want to generate.
	It extracts the dihedrals from the backbone of the loop from the 
	database and builds a new backbone. It will then close the loop
	using the Cyclic Coordinate Descent algorithm. Each closed loop
	it generates is assigned an energetic score based on the CHARMM
	method of computing the Dihedral Angle Potential (V_DA)

	"""
	[first_res,target_coors] = make_first_res()
	CC_DIST = compute_CC_DIST(first_res,target_coors)
	global dihedral_list
	dihedral_list = [] ; LEN_TOL = 2
	start = first_res
	par = 'n'

	if LEN_TOL == 0:
		while len(dihedral_list) < 1:
			[dihedral_list, native_structures] = mine_loop_data(loop_length,CC_DIST,LEN_TOL)
			LEN_TOL += 1
	else:
		[dihedral_list, native_structures] = mine_loop_data(loop_length,CC_DIST,LEN_TOL)

	n_it = len(dihedral_list); print "number of loops to test:", n_it
	success = 0 ; n = 0 ; loopmaker = LoopTools.LoopMaker() ; last_success = 0; V_list = []

	if par == 'y' or par == 'yes':
		Parallel(n_jobs=12)(delayed(run_loop,has_shareable_memory)(n) for n in range(n_it))

	else:
		for di in dihedral_list:
			[first_res,target_coors] = make_first_res()
			loop = loopmaker.finish_grafted_chain(first_res,di)
			[success,V] = run_CCD(extract_coors(loop),target_coors,success,n)
			print "[energy, identification number]", V, n
			if V != 0:
				mda.write_merged(pdbname, "linker"+str(n)+".pdb",n)
				V_list.append([V,n])
			n += 1
			if len(V_list)>max_success or n == len(dihedral_list):
				print success, "structures found"
				for v in V_list:
					print "[energy, identification number, native structure]", v, n, native_structures[n]
				exit()
	return None

################
### OPTION 2 ###
################
def run_straight_aligned(loop_length):
	"""
	This is a much simpler, though significantly less effective, option 
	for generating loops. It builds a straight chain and uses CCD to try
	to close it. If CCD can't close it, you're screwed. 

	"""

	print "WARNING: this method usually doesn't give good results"
	[first_res,target_coors] = make_first_res()
	loopmaker = LoopTools.LoopMaker()
	loop_structure = loopmaker.finish_straight_chain(first_res,loop_length)
	chain = loopmaker.extract_coors(loop_structure)
	run_CCD(chain,target_coors)
	return 0

####################
### execute code ###
####################
if __name__ == "__main__":
	loop_length = 10
	#run_straight_aligned(loop_length)
	run_graft(loop_length)
