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
import os
from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory
from utils import CCD
from utils import LoopTools

"""
          ~~~ TODO ~~~

- make it work for coarse-grained models

"""

#########################################################
### load USER DEFINED variables into global namespace ###
#########################################################

global RMSD_threshold;    RMSD_threshold = 0.25			# cut off for loop closing
global max_iterations;    max_iterations = 1000			# number of iterations for CCD algorithm
global max_success;       max_success = 50			# maximum number of structures you want to make
global DIST_FACTOR;       DIST_FACTOR = 0.2			# N-to-C terminal distance variability for searching ArchDB for structures
global pdbname;           pdbname = "structures/merged8.pdb"
global protein_shift;     protein_shift  = np.array([60,0,0])	# how much you want to shift your input structure
global protein_rotation;  protein_rotation  = 0 		# how much you want to rotate the shifted structure (degrees)
global loop_length;       loop_length = 14
global clash_cutoff;	  clash_cutoff = 0.8

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
	target_coors = loopmaker.make_structure2(pdbname, protein_shift, protein_rotation)
	base_res = loopmaker.define_initial_cond(pdbname,protein_shift)
	first_res = loopmaker.anchor_is_start(base_res)
	return [first_res,target_coors]

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
	dihedral_list = [] ; LEN_TOL = 0
	start = first_res

	if LEN_TOL == 0:
		while len(dihedral_list) < 1:
			[dihedral_list, native_structures] = mine_loop_data(loop_length,CC_DIST,LEN_TOL)
			LEN_TOL += 1
	else:
		[dihedral_list, native_structures] = mine_loop_data(loop_length,CC_DIST,LEN_TOL)

	n_it = len(dihedral_list); print "number of loops to test:", n_it
	success = 0 ; n = 0 ; loopmaker = LoopTools.LoopMaker() ; last_success = 0; V_list = []

	for di in dihedral_list:
		[first_res,target_coors] = make_first_res()
		loop = loopmaker.finish_grafted_chain(first_res,di)
		[success,V] = run_CCD(extract_coors(loop),target_coors,success,n)
		linker_pdb = "linker"+str(n)+".pdb"
		if V != 0:
			mdatools = LoopTools.MDATools(pdbname, linker_pdb, protein_shift)
			overlap = mdatools.check_overlap()
			if overlap < clash_cutoff:
				print "STERIC CLASHES DETECTED. Closest contact =", overlap 
				os.remove(linker_pdb)
			else:
				print "SUCCESS: writing merged structure. Closest contact  =", overlap, ", dihedral energy =", V
				mdatools.write_merged(n)
				V_list.append([V,n])
		if len(V_list)>=max_success or n == len(dihedral_list):
			print success, "structures found"
			for v in V_list:
				print "[energy, identification number, native structure]", v[0], n, native_structures[n]
			exit()
		n += 1

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
	first_res = make_first_res()
	target_coors = loopmaker.make_structure2(pdbname, protein_shift, protein_rotation)
	loopmaker = LoopTools.LoopMaker()
	loop_structure = loopmaker.finish_straight_chain(first_res,loop_length)
	chain = loopmaker.extract_coors(loop_structure)
	run_CCD(chain,target_coors)
	return 0

####################
### execute code ###
####################
if __name__ == "__main__":
	#run_straight_aligned(loop_length)
	run_graft(loop_length)
