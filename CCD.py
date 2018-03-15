#!/usr/bin/env python
"""
Cyclic Coordinate Descent (CCD)
Parts of this class were borrowed from https://github.com/andre-wojtowicz/ccd-loop-closer/blob/master/CCDLoopCloser.py 
Based on article 'Cyclic coordinate descent: A robotics algorithm for protein loop closure.' - Canutescu AA, Dunbrack RL Jr.

This class is an implementation of the cyclic coordinate descent algorithm.
"""
import numpy as np; from numpy import linalg as npl
from numpy.random import random

####################################################################
### This class creates chain objects that are adjusted using CCD ###
####################################################################
class CCDError(Exception):
	pass
class CCD:
	def __init__(self, chain, target):
		self.chain_length = len(chain[:])
		if self.chain_length < 3: return "ERR", "chain length must equal at least 3"
		self.chain = chain[:]
		if len(target) != 3: return "ERR", "target length must equal to 3"
		self.target = target[:]
	
	def write_pdb(self,n):
		outfile = open("linker"+str(n)+".pdb","w")
		it = 0
		for i in self.chain:
			if it%3==0: atom = "N"
			if it%3==1: atom = "CA"
			if it%3==2: atom = "C"
			t1 = "ATOM"					# ATOM
			t2 = 1						# INDEX
			t3 = atom					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = "GLY"					# RESIDUE NAME
			t6 = "X"					# CHAIN
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

	def calc_rmsd(self):
		rmsd = 0.0
		for i in range(1, 4): # iterate through LAST THREE atoms (i.e. [-3, -2, -1] bc of pbc's)
			vec = self.chain[-i] - self.target[-i]
			rmsd += np.dot(vec, vec)
		return np.sqrt(rmsd/3.0)
	
	@classmethod
	def rotate_points(self, M0, O, rot_angle, vec_norm): # EQUATION 3
		R = M0 - O   
		try: R_norm = self.normalize(R)
		except CCDError: return M0 # don't rotate when the point is on rotation axis
		S = np.cross(R_norm, vec_norm)
		return (npl.norm(R) * np.cos(rot_angle) * R_norm + npl.norm(R) * np.sin(rot_angle) * S) + O
	
	@staticmethod
	def get_rotation_axis(P, Q):
		params = {}
		params['x_1']=float(P[0]); params['y_1']=float(P[1]); params['z_1']=float(P[2])
		params['a']=float(Q[0]-P[0]); params['b']=float(Q[1]-P[1]); params['c']=float(Q[2]-P[2])
		return params
	
	@staticmethod
	def get_rotation_central_point(line, P, Q, M_Ox):
		QP = P - Q
		line_xyz = np.array([line['x_1'], line['y_1'], line['z_1']])
		line_abc = np.array([line['a'], line['b'], line['c']])
		t = (np.dot(QP, M_Ox) - np.dot(QP, line_xyz))/float(np.dot(QP, line_abc))
		O_x = np.array([line['x_1'] + line['a'] * t, line['y_1'] + line['b'] * t, line['z_1'] + line['c'] * t])
		return O_x
	
	@staticmethod
	def is_on_rotation_axis(pnt, axis):
		if round(axis['a'], 4) != 0.0: t = (pnt[0] - axis['x_1']) / axis['a']
		elif round(axis['b'], 4) != 0.0: t = (pnt[1] - axis['y_1']) / axis['b']
		elif round(axis['c'], 4) != 0.0: t = (pnt[2] - axis['z_1']) / axis['c']
		else: t = 0.0
		if (pnt[0] == (axis['x_1']+axis['a']*t) and pnt[1] == (axis['y_1']+axis['b']*t) and pnt[2] == (axis['z_1']+axis['c']*t)): return True
		else: return False
	
	@staticmethod
	def normalize(vec):
		"""find unit vector"""
		nrm = npl.norm(vec)
		if nrm == 0.0: raise CCDError('Normalisation error; vector length eq 0')
		else: return vec / nrm

	def check_dihedrals(self):
		for i in range(0,len(self.chain)-8,4): # iterate by RESIDUE
			# OXYGEN ATOMS NOT INCLUDED IN DIHEDRAL CALCULATIONS
			p1 = self.chain[i] 	# N  
			p2 = self.chain[i+1]	# CA
			p3 = self.chain[i+2]	# C
			p5 = self.chain[i+4]	# N
			p6 = self.chain[i+5]	# CA

			def compute_phi(p0,p1,p2,p3):
				b0 = -1.0*(p1 - p0) ; b1 = p2 - p1 ; b2 = p3 - p2
				b0xb1 = np.cross(b0, b1)
				b1xb2 = np.cross(b2, b1)
				b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
				y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
				x = np.dot(b0xb1, b1xb2)
				return np.degrees(np.arctan2(y, x))
			
			#print "phi:", compute_phi(p1,p2,p3,p5), ", psi:", compute_phi(p2,p3,p5,p6)

		return True

	def run(self, n, threshold=0.3, max_it=5000):
		it = 0
		while True:
			rmsd = self.calc_rmsd()
			if rmsd < threshold: 
				write = self.write_pdb(n)
				return "success", rmsd, it
			if it == max_it: return "MAX_IT", rmsd, it

			# iterate through all residues EXCEPT the target ones
			for i in range(0,self.chain_length-2,4):# for almost each edge find best rotation angle
				N_coors = np.array(self.chain[i])
				Ca_coors = np.array(self.chain[i+1])
				rot_axis = self.get_rotation_axis(N_coors, Ca_coors)
				Ca_N_vec_norm = self.normalize(Ca_coors - N_coors)
				b = 0.0; c = 0.0 

				# Now we have a rotation vector, let's see what that does for the RMSD of the target
				for j in range(1, 4):	# these last 3 atoms are the target
					initial_coors = self.chain[-j]; target_coors = self.target[-j]
					rot_point = self.get_rotation_central_point(rot_axis, N_coors, Ca_coors, initial_coors)
					R_j = initial_coors - rot_point
					target_vec = target_coors - rot_point
					if self.is_on_rotation_axis(initial_coors, rot_axis): continue
					# now we have everything we need to compute the rotations (given by b&c)
					S_j = np.cross(self.normalize(R_j), Ca_N_vec_norm) # see Fig 1b
					b += 2 * npl.norm(R_j) * np.dot(target_vec, self.normalize(R_j))
					c += 2 * npl.norm(R_j) * np.dot(target_vec, S_j)
				rot_angle = np.arctan2(c/np.sqrt(b**2+c**2), b/np.sqrt(b**2+c**2)) # EQUATION 10 
				# APPLY ROTATION TO NEXT POINTS
				for j in range(i+2, self.chain_length):
					pntO = self.get_rotation_central_point(rot_axis, N_coors, Ca_coors, self.chain[j])
					chain_copy = self.chain
					self.chain[j] = self.rotate_points(self.chain[j], pntO, rot_angle, Ca_N_vec_norm)
				if self.check_dihedrals() == False:
					print "this isn't going to work"
					return "failed"
			it += 1

