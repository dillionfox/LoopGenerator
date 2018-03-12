import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from pylab import *

RTD = 57.2957795 											# radians to degrees
infile = sys.argv[1]
chains = []
dihedral_file = 'loop_dihedrals.out'

for i in range(2,len(sys.argv)):
	chains.append(sys.argv[i])

### quick and dirty vector manipulation ###

def dot(v1, v2):
	return sum(a * b for a, b in zip(v1, v2))

def cross(v1, v2):
	return [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]

def mag(v):
	return math.sqrt(sum(v[a]*v[a] for a in range(len(v))))

def angle(v1, v2):
	return math.acos(dot(v1, v2)/(mag(v1) * mag(v2)) )

def sub(v1, v2):
	return [v1[i] - v2[i] for i in range(3)]

### calculate dihedrals - math procedure described in link below ###
### www.math.fsu.edu/~quine/MB_10/6_torsion.pdf ###

def calcDihedral(p0, p1, p2, p3):
	V01 = sub(p0, p1)
	V32 = sub(p3, p2)
	V12 = sub(p1, p2)
	vec1 = cross(V12, V01)
	vec2 = cross(V12, V32)
	a = angle(vec1, vec2)
	a = math.degrees(a)
	if dot(cross(vec1, vec2), V12) > 0:
		a = -a
	return a

### separate function to loop through lists of lists containing all of the important info ###

def computePhiPsi(p):
	new_p = []
	for c in p:
		for i in range(1, len(c)):
			c[i][6] = calcDihedral(c[i - 1][5], c[i][3], c[i][4], c[i][5])
		for i in range(0, len(c) - 1):
			c[i][7] = calcDihedral(c[i][3], c[i][4], c[i][5], c[i + 1][3])
		new_p.append(c)
	return new_p

### read pdb, somewhat carefully. Not very robust ###

def readPDB(f, chains):
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
		if resChain in chains: 									# if in the "loop region" and on the correct chain
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
				aa = [curr_Seq, curr_Chain, curr_Type, bb["N"], bb["CA"], bb["C"], None, None]	# store all information about residue
				c.append(aa)									# note: last two entries will store phi & psi
				saved = True
	if c:
		p.append(c)
	return p

def readDihedralFile(dihedral_file):
	phi = []
	psi = []
	for line in open(dihedral_file):
		list = line.split()
		for i in range(len(list)-1):
			phi.append(float(list[i]))
			psi.append(float(list[i+1]))
			i+=2
	return [phi, psi]

def print_output(dihedral_list):									# print angles, don't bother writing to file (debugging)
	for chain in dihedral_list:
		for res in chain:
			#print res[6], res[7]
			flag = 0
			try:
				phi = float(res[6])
				psi = float(res[7])
			except TypeError:
				flag = 1
			if flag != 1:
				print phi, psi

def write_output(dihedral_list):									# save output to file, in format to go with plotting script
	pdb_list = []
	extra = 0
	with open(OUTFILE, "a") as fout:
		for i in dihedral_list:
			fout.write(str(i[0][0]) + " " + str(i[0][1]) + " " + str(i[1][0]) + " " + str(i[1][1]) + " " + str(i[2][0]) + " " + str(i[2][1]) + " " + str(i[3][0]) + " " + str(i[3][1]) + " " + str(i[4][0]) + " " + str(i[4][1]) + "\n")
#				fout.write(str(phi) + " " + str(psi) + "\n")

def rama_plot(dihedral_list, dihedral_file):
	xList = []
	yList = []
	for chain in dihedral_list:
		for res in chain:
			flag = 0
			try:
				phi = float(res[6])
				psi = float(res[7])
			except TypeError:
				flag = 1
			if flag != 1:
				xList.append(phi)
				yList.append(psi)

	xends = [xList[0], xList[1], xList[len(xList)-1],xList[len(xList)-2]]
	yends = [yList[0], yList[1], yList[len(yList)-1],yList[len(yList)-2]]

	phi_file = dihedral_file[0]
	psi_file = dihedral_file[1]

	fig, ax = plt.subplots()
	plot = hist2d(phi_file, psi_file, bins = 40, norm=LogNorm())
	#plot = hist2d(xends, yends, bins = 40, norm=LogNorm())
	#ax.set_xlabel("phi")
	#ax.set_ylabel("psi")
	#cbar = colorbar()
	#cbar.set_label("number of occurances")
	#show()

	#fig = plt.figure()
	#ax = fig.add_subplot(111)
	ax.scatter(xList, yList, s=100, marker='o', c='black')
	ax.scatter(xends, yends, s=100, marker='o', c='white')
	#ax.scatter(phi_file, psi_file, s=0.2, marker='o', c='0.5')
	#ax.scatter(Ou2[0], Ou2[1], Ou2[2], s=30, marker='o', c='black')
	plt.show()


atominfo = readPDB(infile, chains)
dihedral_list = computePhiPsi(atominfo)

dihedral_file_in = readDihedralFile(dihedral_file)

#print_output(dihedral_list)
#write_output(dihedral_list)
rama_plot(dihedral_list, dihedral_file_in)
