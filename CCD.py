#!/usr/bin/env python
"""
Cyclic Coordinate Descent (CCD)
Parts of this class were borrowed from https://github.com/andre-wojtowicz/ccd-loop-closer/blob/master/CCDLoopCloser.py 
Based on article 'Cyclic coordinate descent: A robotics algorithm for protein loop closure.' - Canutescu AA, Dunbrack RL Jr.
"""
import numpy as np; from numpy import linalg as npl
from numpy.random import random

global max_iterations; max_iterations = 300

class CCDError(Exception):
	pass

class CCD:
	def __init__(self, chain, fixed):
		"""
		chain - list of 3+ atoms that are moved.
		fixed - list of target three atoms (C-Anchor).
		"""
		if len(chain) < 3: # initial check of input corectness
			return "ERR", "chain length must equal at least 3"
		if len(fixed) != 3:
			return "ERR", "fixed length must equal to 3"
		self.chain = chain[:]
		self.fixed = fixed[:]
	
	@classmethod
	def arbitrary_rotate(self, pntM, pntO, rot_angle, uvecTheta): # EQUATION 3
		"""
		3D rotation of a point.
		pntM - moving point.
		pntO - central point of rotation.
		rot_angle - rotation angle (radian).
		uvecTheta - unit vector of y in local coordinate system.
		"""
		vecR = pntM - pntO   
		try:
			uvecR = self.normalize(vecR)
		except CCDError: # don't rotate when the point is on rotation axis
			return pntM
		uvecS = np.cross(uvecR, uvecTheta)
		return (npl.norm(vecR) * np.cos(rot_angle) * uvecR + npl.norm(vecR) * np.sin(rot_angle) * uvecS) + pntO

	def write_pdb(self,n):
		outfile = open("linker"+str(n)+".pdb","w")
		for i in self.chain:
			t1 = "ATOM"					# ATOM
			t2 = 1						# INDEX
			t3 = "C"					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = "AAA"					# RESIDUE NAME
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
		outfile.close()
		return 0
	
	def calc_rmsd(self):
		rmsd = 0.0
		for i in range(1, 4): # iterate through LAST THREE atoms (i.e. [-3, -2, -1] bc of pbc's)
			vec = self.chain[-i] - self.fixed[-i]
			rmsd += np.dot(vec, vec)
		return np.sqrt(rmsd/3.0)
	
	@staticmethod
	def get_rotation_axis(pntP, pntQ):
		params = {}
		params['x_1']=float(pntP[0]); params['y_1']=float(pntP[1]); params['z_1']=float(pntP[2])
		params['a']=float(pntQ[0]-pntP[0]); params['b']=float(pntQ[1]-pntP[1]); params['c']=float(pntQ[2]-pntP[2])
		return params
	
	@staticmethod
	def get_rotation_central_point(line, pntP, pntQ, pntM_Ox):
		vecQP = pntP - pntQ
		line_xyz = np.array([line['x_1'], line['y_1'], line['z_1']])
		line_abc = np.array([line['a'], line['b'], line['c']])
		t_numerator = np.dot(vecQP, pntM_Ox) - np.dot(vecQP, line_xyz)
		t_denominator = np.dot(vecQP, line_abc)
		t = t_numerator / float(t_denominator)
		pntO_x = np.array([line['x_1'] + line['a'] * t, line['y_1'] + line['b'] * t, line['z_1'] + line['c'] * t])
		return pntO_x
	
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
		nrm = npl.norm(vec)
		if nrm == 0.0: raise CCDError('Normalisation error; vector length eq 0')
		else: return vec / nrm

	def run(self, plot, threshold=1.0, max_it=5000):
		if plot=='y' or plot=='yes' or plot=='on':
			from mpl_toolkits.mplot3d import Axes3D
			import matplotlib.pyplot as plt
			plt.ion()
			fig = plt.figure()
			ax = fig.add_subplot(111)
			def reshape(C):
				X=[];Y=[];Z=[]
				for q in C:
					X.append(q[0]);Y.append(q[1]);Z.append(q[2])
				return [X,Y,Z]
			[X,Y,Z]=reshape(self.fixed); plt.scatter(X,Y,s=100,marker=(5,1),label="TARGET")
			[X,Y,Z]=reshape(self.chain)
			line1, = ax.plot(X, Y, 'b-')
		it = 0
		while True:
			rmsd = self.calc_rmsd()
			if rmsd < threshold: 
				write = self.write_pdb(it)
				return "success", rmsd, it
			if it == max_it: return "MAX_IT", rmsd, it
			for i in range(len(self.chain)-2):# for almost each edge find best rotation angle
				N_coors = self.chain[i]
				Ca_coors = self.chain[i+1]
				rot_axis = self.get_rotation_axis(N_coors, Ca_coors)
				Ca_N_vec_norm = self.normalize(Ca_coors - N_coors)
				b = 0.0; c = 0.0 # S = (a,b,c)
				for j in range(1, 4):
					initial_coors = self.chain[-j]; target_coors = self.fixed[-j]
					rot_point = self.get_rotation_central_point(rot_axis, N_coors, Ca_coors, initial_coors)
					R_j = initial_coors - rot_point
					target_vec = target_coors - rot_point
					if self.is_on_rotation_axis(initial_coors, rot_axis): continue
					S_j = np.cross(self.normalize(R_j), Ca_N_vec_norm) #see Fig 1b
					b += 2 * npl.norm(R_j) * np.dot(target_vec, self.normalize(R_j))
					c += 2 * npl.norm(R_j) * np.dot(target_vec, S_j)
				rot_angle = np.arctan2(c/np.sqrt(b**2+c**2), b/np.sqrt(b**2+c**2)) # EQUATION 10 
				# APPLY ROTATION TO NEXT POINTS
				for j in range(i+2, len(self.chain)):
					pntO = self.get_rotation_central_point(rot_axis, N_coors, Ca_coors, self.chain[j])
					self.chain[j] = self.arbitrary_rotate(self.chain[j], pntO, rot_angle, Ca_N_vec_norm)
			if plot=='y' or plot=='yes' or plot=='on':
				[X,Y,Z]=reshape(self.chain)
				line1.set_ydata(Y)
				fig.canvas.draw()
				it += 1

def plot(C0,CF,FIXED):
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	def reshape(C):
		X=[];Y=[];Z=[]
		for q in C:
			X.append(q[0]);Y.append(q[1]);Z.append(q[2])
		return [X,Y,Z]
	[X,Y,Z]=reshape(C0); ax.scatter(X,Y,Z,s=200,marker='+',label="INITIAL")
	[X,Y,Z]=reshape(CF); ax.scatter(X,Y,Z,s=100,marker=(5,0),label="FINAL")
	[X,Y,Z]=reshape(FIXED); ax.scatter(X,Y,Z,s=100,marker=(5,1),label="TARGET")
	legend = ax.legend(loc='upper left')
	plt.show()

def make_random_chain_and_run(chain_len,chain_base,target,plot,number_of_loops=10):
	global max_iterations
	success = 0
	av_dist = 3.1 # computed from CG structure
	for n in range(number_of_loops):
		chain = [chain_base]
		current_point = chain_base
		for j in range(chain_len):
			next_point = random(3)
			next_point = next_point/npl.norm(next_point)
			next_point = current_point + next_point*av_dist
			chain.append(next_point)
			current_point = next_point
		for j in range(-5, -2):
			phi = 2*np.pi*random()
			rot_axis = CCD.get_rotation_axis(chain[j], chain[j+1])
			vec_PQ = chain[j+1] - chain[j]
			vec_PQ_norm = vec_PQ / npl.norm(vec_PQ)			
			for k in range(j+2, 0):
				center_point = CCD.get_rotation_central_point(rot_axis, chain[j], chain[j+1], chain[k])
				chain[k] = CCD.arbitrary_rotate(chain[k], center_point, phi, vec_PQ_norm)
		init = CCD(chain, target)
		check, RMSD, it = init.run(plot,max_it=max_iterations)
		if check == 'success': success += 1
		print("{0}\t{1}\t\t{2}\t{3}\t{4}".format(n, str(round(success*100/float(n+1),2))+'%', check, str(RMSD), it))

def run_CCD():
	print "need to constrain the first 3 beads to overlap with the C terminus"
	print "also need to not print the first 3 and last 3 beads since they are"
	print "supposed to overlap"
	chain_base = np.array([9.840000,16.639999,25.469999])
	target = [np.array([27.95,28.18,31.82]),np.array([29.68,28.79,34.75]),np.array([30.23,29.99,37.54])]
	chain_length = 12
	number_of_loops = 10
	plot='off'
	make_random_chain_and_run(chain_length,chain_base,target,plot,number_of_loops)

if __name__ == "__main__":
	run_CCD()
