##############################################################################
##############################################################################
###                                                                        ###  
### This class mines a loop database (ArchDB) and extracts the information ###
###     needed to reconstruct the backbone of loops that meet geometric    ###
###                              criteria                                  ###
###                                                                        ###
###                         Dillion Fox, 3/2018                            ###
###                                                                        ###  
##############################################################################
##############################################################################

import urllib
import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from pylab import *
import MDAnalysis as mda
import MDAnalysis.analysis.distances as mdad

url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" 	# link to RCSB
RTD = 57.2957795 											# radians to degrees

class LoopMiner:
	"""
	Collection of functions useful for mining ArchDB

	"""

	def __init__(self,numResLoop,CC_DIST,DIST_TOL,OUTFILE):
		self.numResLoop = numResLoop
		self.CC_DIST = CC_DIST
		self.DIST_TOL = DIST_TOL
		self.OUTFILE = OUTFILE

	### quick and dirty vector manipulation ###
	def dot(self, v1, v2): return sum(a * b for a, b in zip(v1, v2))
	def cross(self, v1, v2): return [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]
	def mag(self, v): return math.sqrt(sum(v[a]*v[a] for a in range(len(v))))
	def angle(self, v1, v2): return math.acos(dot(v1, v2)/(self.mag(v1) * self.mag(v2)) )
	def sub(self, v1, v2): return [v1[i] - v2[i] for i in range(3)]

	def calcDihedral(self, p0, p1, p2, p3):
		"""
		compute dihedral from 4 points using function that used to be listed on Wikipedia

		"""

		V01 = self.sub(p0, p1)
		V32 = self.sub(p3, p2)
		V12 = self.sub(p1, p2)
		vec1 = self.cross(V12, V01)
		vec2 = self.cross(V12, V32)
		a = self.angle(vec1, vec2)
		a = math.degrees(a)
		if self.dot(cross(vec1, vec2), V12) > 0: a = -a
		return a

	def computePhiPsi(self,p):
		"""
		separate function to loop through lists of lists containing all of the important info

		""" 

		new_p = []
		for c in p:
			for i in range(1, len(c)):
				c[i][6] = self.calcDihedral(c[i-1][5], c[i][3], c[i][4], c[i][5])
			for i in range(0, len(c) - 1):
				c[i][7] = self.calcDihedral(c[i][3], c[i][4], c[i][5], c[i + 1][3])
			new_p.append(c)
		return new_p
	
	@staticmethod
	def readPDB(f, start, stop, chain):
		"""
		read pdb, somewhat carefully. Not very robust.

		"""

		p = []												# stores information about all residues
		c = []												# stores all information about EACH residue 
		curr_Seq = None
		curr_Chain = None
		curr_Type = None	
		bb = {}												# stores backbone atoms
		saved = False											# new aa
		for line in open(f):
			if line[:3] == "TER": 									# if end of chain
				if c:
					p.append(c)
					c = []
				continue
			if line[:4] != "ATOM":
				continue
			atomName = line[12:16].strip()
			if atomName not in [ "N", "CA", "C" ]:							# Ignore non-backbone atoms
				continue
			resType = line[17:20]									#
			resChain = line[21]									# read data from pdb
			resSeq = int(line[22:26])								#
			if resSeq >= start and resSeq <= stop and resChain == chain: 				# if in the "loop region" and on the correct chain
				if curr_Seq != resSeq or curr_Chain != resChain or curr_Type != resType:	# if new residue, reset variables
					curr_Seq = resSeq
					curr_Type = resType
					curr_Chain = resChain
					bb = {}
					saved = False
				x = float(line[30:38])								#
				y = float(line[38:46])								# assign coordinates to variables
				z = float(line[46:54])								#
				bb[atomName] = [x, y, z]							# store vector in list
				if len(bb) == 3 and not saved:							# if new residue and all bb atoms were identified
					aa = [curr_Seq, curr_Chain, curr_Type, bb["N"], bb["CA"], bb["C"], None, None, f]	# store all information about residue
					c.append(aa)									# note: last two entries will store phi & psi
					saved = True
		if c:
			p.append(c)
		return p

	@staticmethod
	def download_pdb(pdbID):
		"""
		download pdb from RCSB	

		"""

		pdb = url+str(pdbID)
		open( pdbID+".pdb", "w" ).write( urllib.urlopen(pdb).read() )
		return 0

	### "Main Function" ###
	def mineArchDB(self,f):											# read file downloaded from ARCH DB
		"""
		This is the "main" function. It iterates through the ArchDB file (that the user supplies!),
		downloads the PDB's, extracts the dihedrals, and deletes the PDB file.

		"""

		pdb_list = []											# store PDB ID's
		dihedral_list = []										# store dihedral angles
		old_pdbID = 'null'
		total=0												# track number of loops (for debugging)
		minLength = 4
		maxLength = 20
		for line in open(f):
			flag = 0										# if information cannot be salvaged, this variable will break loop
			list = line.split()									# break lines delimited by spaces
			if (list[5] == "HH" and int(list[6]) in self.numResLoop):# and total < 3:			# if loop is between two helices and has appropriate number of residues
				pdbID = list[1]									# store information about loop
				chain = list[2]									#
				try:										# in case information cannot be read, use "try" statement
					SSstart = float(list[3])						# secondary structure ends here
					SSstop = float(list[4])							# secondary structure starts here again after loop (confusing ...)
					Nterm = float(list[7])							# this many residues between SS end and loop starting
					Cterm = float(list[8])							# this many residues between loop ending and SS starting
					dist = float(list[9])
				except:
					pass
					flag = 1
	
				if dist < (self.CC_DIST + self.DIST_TOL) and dist > (self.CC_DIST - self.DIST_TOL):
					total+=1
					loop_start = SSstart + Nterm							# start of loop (residue number)
					loop_stop = SSstop - Cterm							# end of loop (residue number)
					if (loop_stop - loop_start) > minLength and (loop_stop - loop_start) < maxLength:
						if pdbID not in pdb_list and flag != 1:						# if new pdb, download it
							if os.path.isfile(old_pdbID+'.pdb'): 					# delete old pdb before downloading new one
								os.remove(old_pdbID+'.pdb')
							pdb_list.append(pdbID)							# keep track of pdb's used
							self.download_pdb(pdbID)							# download pdb
							p = self.readPDB(pdbID+'.pdb', loop_start, loop_stop, chain)			# read pdb, which calls functions to calculate phi/psi
							curr_loop_dihedrals = self.computePhiPsi(p)					# save dihedral angles
							old_pdbID = pdbID
			
						elif pdbID in pdb_list and flag != 1:						# if pdb was already downloaded, use it
							p = self.readPDB(pdbID+'.pdb', loop_start, loop_stop, chain)			# read pdb, which calls functions to calculate phi/psi
							curr_loop_dihedrals = self.computePhiPsi(p)
			
						for List in curr_loop_dihedrals:						# take each pair of angles and add to total list
							loop_dihedrals = []
							type_error = 0
							for k in List:
								type_error = 0
								try:
									phi = float(k[6])
									psi = float(k[7])
								except TypeError:
									type_error = 1
								if type_error != 1:
									loop_dihedrals.append([phi, psi])
						dihedral_list.append(loop_dihedrals)					# comprehensive list of pairs of angles
		return [dihedral_list, pdb_list]

	### Output options ###
	def write_output(dihedral_list):									# save output to file, in format to go with plotting script
		"""
		Save the output in a format that is readable by the LoopMaker functions

		"""

		pdb_list = []
		extra = 0
		with open(self.OUTFILE, "a") as fout:
			for i in dihedral_list:
				for j in i:
					fout.write(str(j[0])+" "+str(j[1])+" ")
				fout.write("\n")

def run_LoopMiner():
	"""
	here are some functions I used while writing/debugging this class

	"""

	def print_output(dihedral_list):									# print angles, don't bother writing to file (debugging)
		"""
		self explanatory

		"""

		for i in dihedral_list:
			flag = 0
			try:
				phi = float(i[0])
				psi = float(i[1])
			except TypeError:
				flag = 1
			if flag != 1:
				print phi, psi
	
	def read_outfile(outfile):
		"""
		This function reads the file that is output by the class above

		"""

		dihedral_list = []
		for line in open(outfile):
			list = line.split()									# break lines delimited by spaces
			x = []
			for l in list:
				x.append(float(l))
			dihedral_list.append(x)
		return dihedral_list
	
	def rama_plot(dihedral_list):
		"""
		Generate Ramachandran heatmap

		"""

		xList = []
		yList = []
		for d in dihedral_list:
			for i in range(0,len(d),2):
				xList.append(d[i])
				yList.append(d[i+1])
		fig, ax = plt.subplots()
		plot = hist2d(xList, yList, bins = 30, norm=LogNorm(),cmap='hot')
		ax.set_xlabel("phi")
		ax.set_ylabel("psi")
		cbar = colorbar()
		cbar.set_label("number of occurances")
		show()
	if os.path.exists(self.OUTFILE) == True:
		dihedral_list = read_outfile(OUTFILE)	
		rama_plot(dihedral_list)
	else:
		dihedral_list = mineArchDB(filename)									# mineArchDB is basically main function, it starts everything and ends there too
		#print_output(dihedral_list)										# output options
		write_output(dihedral_list)
		dihedral_list = []
		dihedral_list = read_outfile(OUTFILE)	
		rama_plot(dihedral_list)

##############################################################################
##############################################################################
###                                                                        ###  
###  This class builds loops piggybacking off of the PeptideBuilder class  ###
###                     and the BioPython package                          ###
###                                                                        ###  
###                        Dillion Fox, 3/2018                             ###
###                                                                        ###  
##############################################################################
##############################################################################

import os
import numpy as np
import Bio.PDB

class LoopMaker:
	"""
	This is the main collection of functions I wrote to generate loops using 2 methods.
	The most useful method of making loops is to find an existing one in the PDB with
	similar geometry to what you want and graft it onto your structure. The match won't
	be perfect, so you can use the Cyclic Coordinate Descent algorithm in the CCD.py
	file to close it for you. Note that this code only produces GLYCINE chains.

	"""

        def __init__(self):
		return None

	def write_pdb(self,coors,pdbname):
		"""
		I find BioPython to be annoying, so this and many of the other functions
		in this file only serve to avoid using BioPython

		"""

		outfile = open(pdbname,"w")
		it = 0
		# loop contains extra residue at beginning and extra residue at end. Don't put those in PDB
		for i in coors:
			if it%4==0: atom = "N"
			if it%4==1: atom = "CA"
			if it%4==2: atom = "C"
			if it%4==3: atom = "O"
			t1 = "ATOM"					# ATOM
			t2 = 1						# INDEX
			t3 = atom					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = "GLY"					# RESIDUE NAME
			t6 = "A"					# CHAIN
			t7 = 0						# RESIDUE NUMBER
			t8 = ""						# INSERTION CODE
			t9 =  float(i[0])				# X
			t10 = float(i[1])				# Y
			t11 = float(i[2])				# Z
			t12 = 0.0					# OCCUPANCY
			t13 = 0.0					# TEMPERATURE FACTOR
			t14 = ""					# ELEMENT SYMBOL
			t15 = ""					# CHARGE ON ATOM
			outfile.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15))
			it+=1
		outfile.close()
		return 0

	def define_initial_cond(self,pdbname,protein_shift):
		"""
		The first thing you have to define is the ANCHOR (C-terminus of one protein)
		and the TARGET (N-terminus of the protein you're connecting to). This function
		returns the C-terminal residue as a BioPython Structure, and the coordinates
		of the target atom. The BioPython Structure is required because it serves as
		a "dummy" residue. It will be the first residue on the chain to ensure that the 
		chain connects to the protein properly. It will late be deleted.

		"""

		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		model = pdb_parser.get_structure("ref",pdbname)[0]
		#last_res =  self.model_selectres(model,self.model_numres(model))
		last_res = model['A'][self.model_numres(model)]
		return last_res

	def make_structure2(self, pdbname, protein_shift, protein_rotate):
		"""
		This is a substitute for the function above which doesn't work. It uses MDAnalysis...

		"""

		prot2 = mda.Universe(pdbname)
		prot2.atoms.translate(protein_shift)
		prot2.atoms.rotateby(protein_rotate,[0,0,1])
		prot2.atoms.write(pdbname+"_2.pdb")
		sel = prot2.select_atoms('resid 1 and (name N or name CA or name C)')
		return sel.positions
	
	def extract_coors(self,structure):
		"""
		Get the coordinates out of the stupid BioPython Structure

		"""

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
	
	def finish_straight_chain(self,loop,loop_length):
		"""
		Make a perfectly straight glycine chain building off of the C-terminus
		of your structure

		"""

		for i in range(1,loop_length):
			loop = add_residue(loop,"G")
		return loop

	def finish_grafted_chain(self,loop,di):
		"""
		Make a chain with dihedral angles specified by a structure found in the PDB

		"""

		for i in range(len(di)):
			#print di[i][0], di[i][1]
			loop = add_residue(loop,"G",di[i][0],di[i][1])
		return loop
	
	def model_numres(self,model):
		"""
		I hate BioPython

		"""

		length = 0
		for residue in model.get_residues():
			length+=1
		return length
	
	#def model_selectres(self,model,resnum):
	#	"""
	#	Because I sometimes forget

	#	"""
	#	import ipdb; ipdb.set_trace()
	#	return model['A'][resnum]
	
	def res_get_atom_coors(self,res):
		"""
		I hate BioPython

		"""
		v = []
		v.append(np.array([res['N' ].get_coord()[0],res['N' ].get_coord()[1],res['N' ].get_coord()[2]]))
		v.append(np.array([res['CA'].get_coord()[0],res['CA'].get_coord()[1],res['CA'].get_coord()[2]]))
		v.append(np.array([res['C' ].get_coord()[0],res['C' ].get_coord()[1],res['C' ].get_coord()[2]]))
		return v
	
	def get_bb_atoms(self,model):
		"""
		I hate BioPython

		"""
		atoms = []
		for chain in model:
			for res in chain:
				atoms.append(res['N'])
				atoms.append(res['CA'])
				atoms.append(res['C'])
		return np.array(atoms)
	
	def anchor_is_start(self,base_res):
		"""
		Dumb name, useful function. I couldn't figure out how to delete the class instance defining the 
		BioPython structure, so when I'd make many loops they all built off of eachother instead of
		starting from scratch. Therefore I re-use this function in a loop to keep generating new
		starting points in the loop. It's a memory hog, but it's a decent fix for now.

		"""

		# biopython is stupid so I have to extract coordinates, save in pdb, and reload to make coordinates into a "structure"
		out = Bio.PDB.PDBIO()
		pdbname = "anchor_res.pdb"
		coors = np.array([base_res['N'].get_coord(), base_res['CA'].get_coord(), base_res['C'].get_coord(), base_res['O'].get_coord()])
		self.write_pdb(coors,pdbname)
		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		loop_structure = pdb_parser.get_structure("ref",pdbname)
		os.remove( pdbname ) # delete unnecessary pdb

		return loop_structure

def run_LoopMaker():
	"""
	Suggested use case

	"""

	pdbname = "/home/dillion/data/reflectin/add_linker/structures/small_mem_protein-inserted_frame-738.pdb"
	protein_shift = np.array([30,0,0])
	loop_length = 15
	loopmaker = LoopMaker()
	[base_res, target_coors] = loopmaker.define_initial_cond(pdbname,protein_shift)
	loop_structure = loopmaker.anchor_is_start(base_res)
	loop_structure = loopmaker.finish_straight_chain(loop_structure)
	coors = loopmaker.extract_coors(loop_structure)
	print coors


##############################################################################
##############################################################################
###                                                                        ###
###  This class employs MDAnalysis functions to manipulate the linker and  ###
###   protein structures. I would have used BioPython for these features   ###
###          but I found MDAnalysis was much better suited for it.         ###
###                                                                        ###
###                                                                        ###
###                         Dillion Fox, 3/2018                            ###
###                                                                        ###  
##############################################################################
##############################################################################

class MDATools:
	def __init__(self, protein_pdb, linker_pdb, shift):
		self.protein_pdb = protein_pdb
		self.linker_pdb = linker_pdb
		self.shift = shift
		self.rotate = 30
		self.protein_pdb2 = protein_pdb+"_2.pdb"

	def merge_structures(self):
		"""
		This function takes the original protein structure and the newly created linker
		and merges them into one MDAnalysis universe.
	
		"""
		
		prot1 = mda.Universe(self.protein_pdb) ; pl = len(prot1.atoms.residues)
		prot2 = mda.Universe(self.protein_pdb2)
		link = mda.Universe(self.linker_pdb) ; ll = len(link.atoms.residues)
		
		# protein 1
		sel1 = 'resid '
		for i in range(pl):
			sel1+=str(i)+' '
		
		# linker
		sel2 = 'resid '
		for i in range(ll):
			link.atoms.residues[i].resid += (pl+1)
			sel2+=str(link.atoms.residues[i].resid)+' '
		
		# protein 2
		sel3 = 'resid '
		for i in range(pl):
			prot2.atoms.residues[i].resid += (pl+ll)
			sel3+=str(prot2.atoms.residues[i].resid)+' '
		
		u = mda.Merge(prot1.atoms,link.atoms,prot2.atoms)
		sel = [sel1, sel2, sel3]
		return [u, sel]
	
	def check_overlap(self):
		"""
		This function checks for steric clashes between the two proteins and the linker.
		If clashes are found, it returns False, otherwise it returns True.
	
		"""
	
		def find_min(arr):
			epsilon = 0.01
			el = np.where(np.array(arr) < arr.min()+epsilon)
			return [el, arr.min()]

		[u, sel] = self.merge_structures()

		newsel = 'resid '
		for res in sel[1].split()[2:-1]:
			newsel += str(res)+' '
		
		protein1 = u.select_atoms(sel[0])
		#linker2 = u.select_atoms(sel[1])
		linker2 = u.select_atoms(newsel)
		protein3 = u.select_atoms(sel[2])
		
		dist12 = mdad.distance_array(protein1.positions,linker2.positions)
		dist23 = mdad.distance_array(protein3.positions,linker2.positions)
		dist13 = mdad.distance_array(protein1.positions,protein3.positions)
		
		return min(dist12.min(), dist23.min(), dist13.min()) 
	
	def write_merged(self, n):
		u = self.merge_structures()[0]
		u.atoms.write("merged"+str(n)+".pdb")
		return None


##############################################################################
##############################################################################
###                                                                        ###
###     !!! THIS IS A MODIFIED VERSION OF : PeptideBuilder library !!!     ###
###            written by Matthew Z. Tien, Dariya K. Sydykova,             ###
###                 Austin G. Meyer, and Claus O. Wilke.                   ###
###                   (Protected under GNU license)                        ###
###                                                                        ### 
###         This class contains a bunch of useful functions that           ###
###       piggyback on Biopython functions that generate residues.         ###
###          The Biopython classes are very annoying to use and            ###
###               manipulate so this just makes it easier.                 ###
###                                                                        ###
###                         Dillion Fox, 3/2018                            ###
###                                                                        ###  
##############################################################################
##############################################################################

from Bio.PDB import *
from Bio.PDB.Atom import *
from Bio.PDB.Residue import *
from Bio.PDB.Chain import *
from Bio.PDB.Model import *
from Bio.PDB.Structure import *
from Bio.PDB.Vector import *
from Bio.PDB.Entity import*
import math, warnings

class Geo():
	def __repr__(self):
		repr = ""
		for var in dir(self):
			if var in self.__dict__: # exclude member functions, only print member variables
				repr += "%s = %s\n" % ( var, self.__dict__[var] )
		return repr

class GlyGeo(Geo):
	def __init__(self):
		self.CA_N_length=1.46
		self.CA_C_length=1.52
		self.N_CA_C_angle=110.8914
		self.C_O_length=1.23
		self.CA_C_O_angle=120.5117
		self.N_CA_C_O_diangle= 180.0
		self.phi=-120
		self.psi_im1=140
		self.omega=180.0
		self.peptide_bond=1.33
		self.CA_C_N_angle =116.642992978143
		self.C_N_CA_angle= 121.382215820277
		self.residue_name= 'G'

def geometry():
	return GlyGeo()

def calculateCoordinates(refA, refB, refC, L, ang, di):
	AV=refA.get_vector(); BV=refB.get_vector(); CV=refC.get_vector()
	CA=AV-CV; CB=BV-CV
	##CA vector
	AX=CA[0]; AY=CA[1]; AZ=CA[2]
	##CB vector
	BX=CB[0]; BY=CB[1]; BZ=CB[2]
	##Plane Parameters
	A=(AY*BZ)-(AZ*BY); B=(AZ*BX)-(AX*BZ); G=(AX*BY)-(AY*BX)
	##Dot Product Constant
	F= math.sqrt(BX*BX + BY*BY + BZ*BZ) * L * math.cos(ang*(math.pi/180.0))
	##Constants
	const=math.sqrt( math.pow((B*BZ-BY*G),2) *(-(F*F)*(A*A+B*B+G*G)+(B*B*(BX*BX+BZ*BZ) + A*A*(BY*BY+BZ*BZ)- (2*A*BX*BZ*G) + (BX*BX+ BY*BY)*G*G - (2*B*BY)*(A*BX+BZ*G))*L*L))
	denom= (B*B)*(BX*BX+BZ*BZ)+ (A*A)*(BY*BY+BZ*BZ) - (2*A*BX*BZ*G) + (BX*BX+BY*BY)*(G*G) - (2*B*BY)*(A*BX+BZ*G)
	X= ((B*B*BX*F)-(A*B*BY*F)+(F*G)*(-A*BZ+BX*G)+const)/denom
	if((B==0 or BZ==0) and (BY==0 or G==0)):
		const1=math.sqrt( G*G*(-A*A*X*X+(B*B+G*G)*(L-X)*(L+X)))
		Y= ((-A*B*X)+const1)/(B*B+G*G)
		Z= -(A*G*G*X+B*const1)/(G*(B*B+G*G))
	else:
		Y= ((A*A*BY*F)*(B*BZ-BY*G)+ G*( -F*math.pow(B*BZ-BY*G,2) + BX*const) - A*( B*B*BX*BZ*F- B*BX*BY*F*G + BZ*const)) / ((B*BZ-BY*G)*denom)
		Z= ((A*A*BZ*F)*(B*BZ-BY*G) + (B*F)*math.pow(B*BZ-BY*G,2) + (A*BX*F*G)*(-B*BZ+BY*G) - B*BX*const + A*BY*const) / ((B*BZ-BY*G)*denom)
	#GET THE NEW VECTOR from the orgin
	D=Vector(X, Y, Z) + CV
	with warnings.catch_warnings():
		# ignore inconsequential warning
		warnings.simplefilter("ignore")
		temp=calc_dihedral(AV, BV, CV, D)*(180.0/math.pi)
	di=di-temp
	rot= rotaxis(math.pi*(di/180.0), CV-BV)
	D=(D-BV).left_multiply(rot)+BV
	return D

def makeGly(segID, N, CA, C, O, geo):
	'''Creates a Glycine residue'''
	##Create Residue Data Structure
	res= Residue((' ', segID, ' '), "GLY", '    ')
	res.add(N)
	res.add(CA)
	res.add(C)
	res.add(O)
	return res

def getReferenceResidue(structure):
	resRef = structure[0]['A'].child_list[-1]
	assert is_aa(resRef)
	return resRef

def add_residue_from_geo(structure, geo):
	resRef= getReferenceResidue(structure)
	AA=geo.residue_name
	segID= resRef.get_id()[1]
	segID+=1
	##geometry to bring together residue
	peptide_bond=geo.peptide_bond
	CA_C_N_angle=geo.CA_C_N_angle
	C_N_CA_angle=geo.C_N_CA_angle
	##Backbone Coordinages
	N_CA_C_angle=geo.N_CA_C_angle
	CA_N_length=geo.CA_N_length
	CA_C_length=geo.CA_C_length
	phi= geo.phi
	psi_im1=geo.psi_im1
	omega=geo.omega
	N_coord=calculateCoordinates(resRef['N'], resRef['CA'], resRef['C'], peptide_bond, CA_C_N_angle, psi_im1)
	N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")
	CA_coord=calculateCoordinates(resRef['CA'], resRef['C'], N, CA_N_length, C_N_CA_angle, omega)
	CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")
	C_coord=calculateCoordinates(resRef['C'], N, CA, CA_C_length, N_CA_C_angle, phi)
	C= Atom("C", C_coord, 0.0, 1.0, " ", " C",0,"C")
	##Create Carbonyl atom (to be moved later)
	C_O_length=geo.C_O_length
	CA_C_O_angle=geo.CA_C_O_angle
	N_CA_C_O_diangle=geo.N_CA_C_O_diangle
	carbonyl=calculateCoordinates(N, CA, C, C_O_length, CA_C_O_angle, N_CA_C_O_diangle)
	O= Atom("O",carbonyl , 0.0 , 1.0, " "," O", 0, "O")
	if(AA=='G'):
		res=makeGly(segID, N, CA, C, O, geo)
	else:
		res=makeGly(segID, N, CA, C, O, geo)
	resRef['O'].set_coord(calculateCoordinates(res['N'], resRef['CA'], resRef['C'], C_O_length, CA_C_O_angle, 180.0))
	ghost= Atom("N", calculateCoordinates(res['N'], res['CA'], res['C'], peptide_bond, CA_C_N_angle, psi_im1), 0.0 , 0.0, " ","N", 0, "N")
	res['O'].set_coord(calculateCoordinates( res['N'], res['CA'], res['C'], C_O_length, CA_C_O_angle, 180.0))
	structure[0]['A'].add(res)
	return structure
    
def add_residue(structure, residue, phi=-120, psi_im1=140, omega=-370):
	if isinstance( residue, Geo ):
		geo = residue
	else:
		geo=geometry() 
		geo.phi=phi
		geo.psi_im1=psi_im1
		if omega>-361:
			geo.omega=omega
	add_residue_from_geo(structure, geo)
	return structure
