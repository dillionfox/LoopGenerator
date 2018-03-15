#!/usr/bin/env python
"""
Cyclic Coordinate Descent (CCD)
Parts of this class were borrowed from https://github.com/andre-wojtowicz/ccd-loop-closer/blob/master/CCDLoopCloser.py 
Based on article 'Cyclic coordinate descent: A robotics algorithm for protein loop closure.' - Canutescu AA, Dunbrack RL Jr.
"""
import numpy as np; from numpy import linalg as npl
from numpy.random import random
import Bio.PDB
import CCD
import LoopMaker

############################################
### load variables into global namespace ###
############################################
global max_iterations; max_iterations = 300
global dl; dl = 0 # variability in chain length. i.e. 15 +/- 1
global DIST_FACTOR; DIST_FACTOR = 0.05

########################################################
### Graft existing loops from the pdb onto structure ###
########################################################
def mine_loop_data(loop_length,CC_DIST,LEN_TOL):
	#import loop_miner
	OUTFILE = "loop_dihedrals.out"
	ArchDB = "ArchDB/archdb.mini.tab"
	DIST_TOL = np.floor(CC_DIST*DIST_FACTOR)
	numResLoop = range(int(loop_length-LEN_TOL),int(loop_length+LEN_TOL+1))						# desired number of residues in loop
	print "looking for loops with", loop_length, "+/-", LEN_TOL, ", and span a distance of", CC_DIST, "+/-", DIST_TOL

	miner = LoopMaker.LoopMiner(numResLoop,CC_DIST,DIST_TOL,OUTFILE)
	dihedral_list = miner.mineArchDB(ArchDB)							# mineArchDB is basically main function, it starts everything and ends there too
	return dihedral_list

def extract_coors(structure):
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
	### RUN CCD ###
	init = CCD.CCD(chain, target_coors) 
	check, RMSD, it = init.run(n,max_it=max_iterations)
	if check == 'success': success += 1
	print("{0}\t{1}\t\t{2}\t{3}\t{4}".format(n, str(round(success*100/float(n+1),2))+'%', check, str(RMSD), it))
	return success

def run_graft(loop_length):
	#import loop_miner
	[first_res,target_coors] = make_first_res()
	CC_DIST = compute_CC_DIST(first_res,target_coors)
	dihedral_list = [] ; LEN_TOL = 0
	start = first_res

	while len(dihedral_list) < 1:
		dihedral_list = mine_loop_data(loop_length,CC_DIST,LEN_TOL)
		LEN_TOL += 1

	success = 0 ; n = 0 ; loopmaker = LoopMaker.LoopMaker()
	for di in dihedral_list:
		[first_res,target_coors] = make_first_res()
		loop = loopmaker.finish_grafted_chain(first_res,di)
		success = run_CCD(extract_coors(loop),target_coors,success,n)
		n += 1
	return None

###############################################################################################
### this function starts the chain with the appropriate dihedrals, direction, and placement ###
###############################################################################################
def run_straight_aligned(loop_length):
	print "this usually doesn't give good results"
	[first_res,target_coors] = make_first_res()
	loopmaker = LoopMaker.LoopMaker()
	loop_structure = loopmaker.finish_straight_chain(first_res,loop_length)
	chain = loopmaker.extract_coors(loop_structure)

	run_CCD(chain,target_coors)
	return 0

def make_first_res():
	pdbname = "/home/dillion/data/reflectin/add_linker/structures/small_mem_protein-inserted_frame-738.pdb"
	protein_shift = np.array([30,0,0])
	loopmaker = LoopMaker.LoopMaker()
	[base_res, target_coors] = loopmaker.define_initial_cond(pdbname,protein_shift)
	first_res = loopmaker.anchor_is_start(base_res)
	return [first_res,target_coors]

####################
### execute code ###
####################
if __name__ == "__main__":
	loop_length = 10
	#run_straight_aligned(loop_length)
	run_graft(loop_length)
