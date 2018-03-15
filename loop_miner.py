# This code reads a file downloaded from arch DB (http://sbi.imim.es/archdb/), a database containing
# lots of information about loops. Download whichever file you are interested in, and specify the path of the file:

import urllib
import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from pylab import *

url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=" 	# link to RCSB
RTD = 57.2957795 											# radians to degrees

class LoopMiner:
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
		V01 = self.sub(p0, p1)
		V32 = self.sub(p3, p2)
		V12 = self.sub(p1, p2)
		vec1 = self.cross(V12, V01)
		vec2 = self.cross(V12, V32)
		a = self.angle(vec1, vec2)
		a = math.degrees(a)
		if self.dot(cross(vec1, vec2), V12) > 0: a = -a
		return a
	### separate function to loop through lists of lists containing all of the important info ###
	def computePhiPsi(self,p):
		new_p = []
		for c in p:
			for i in range(1, len(c)):
				c[i][6] = self.calcDihedral(c[i-1][5], c[i][3], c[i][4], c[i][5])
			for i in range(0, len(c) - 1):
				c[i][7] = self.calcDihedral(c[i][3], c[i][4], c[i][5], c[i + 1][3])
			new_p.append(c)
		return new_p
	
	### read pdb, somewhat carefully. Not very robust ###
	@staticmethod
	def readPDB(f, start, stop, chain):
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
	def download_pdb(pdbID):										# download pdb from RCSB	
		pdb = url+str(pdbID)
		open( pdbID+".pdb", "w" ).write( urllib.urlopen(pdb).read() )
		return 0

	### "Main Function" ###
	def mineArchDB(self,f):											# read file downloaded from ARCH DB
		pdb_list = []											# store PDB ID's
		dihedral_list = []										# store dihedral angles
		old_pdbID = 'null'
		total=0												# track number of loops (for debugging)
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
					if (loop_stop - loop_start) > 6 and (loop_stop - loop_start) < 10:
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
		return dihedral_list

	### Output options ###
	def write_output(dihedral_list):									# save output to file, in format to go with plotting script
		pdb_list = []
		extra = 0
		with open(self.OUTFILE, "a") as fout:
			for i in dihedral_list:
				for j in i:
					fout.write(str(j[0])+" "+str(j[1])+" ")
				fout.write("\n")

if __name__ == "__main__":
	def print_output(dihedral_list):									# print angles, don't bother writing to file (debugging)
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
		dihedral_list = []
		for line in open(outfile):
			list = line.split()									# break lines delimited by spaces
			x = []
			for l in list:
				x.append(float(l))
			dihedral_list.append(x)
		return dihedral_list
	
	def rama_plot(dihedral_list):
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
		print "MOVING ON"
		rama_plot(dihedral_list)
	else:
		dihedral_list = mineArchDB(filename)									# mineArchDB is basically main function, it starts everything and ends there too
		#print_output(dihedral_list)										# output options
		write_output(dihedral_list)
		dihedral_list = []
		dihedral_list = read_outfile(OUTFILE)	
		rama_plot(dihedral_list)
