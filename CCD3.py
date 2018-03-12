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
	def __init__(self, chain, target):
		if len(chain) < 3: return "ERR", "chain length must equal at least 3"
		self.chain = chain[:]
		if len(target) != 3: return "ERR", "target length must equal to 3"
		self.target = target[:]
	
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

	@staticmethod
	def compute_phi(p):
		"""formula from Wikipedia article on "Dihedral angle"; formula was removed"""
		p0 = p[0] ; p1 = p[1] ; p2 = p[2] ; p3 = p[3]
		b0 = -1.0*(p1 - p0) ; b1 = p2 - p1 ; b2 = p3 - p2
		b0xb1 = np.cross(b0, b1)
		b1xb2 = np.cross(b2, b1)
		b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
		y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
		x = np.dot(b0xb1, b1xb2)
		return np.degrees(np.arctan2(y, x))

	def compute_dihedrals(self, chain):
		"""iterate through chain and compute phi and psi"""
		for i in range(len(chain)-4):
			p = [chain[i],chain[i+1],chain[i+2],chain[i+3]]
			if i%3==0:
				print "psi: ", self.compute_phi(p)
			if i%3==1:
				print "phi: ", self.compute_phi(p)
			if i%3==2:
				print "omega: ", self.compute_phi(p)
		return 0

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
			[X,Y,Z]=reshape(self.target); plt.scatter(X,Y,s=100,marker=(5,1),label="TARGET")
			[X,Y,Z]=reshape(self.chain)
			line1, = ax.plot(X, Y, 'b-')
		it = 0
		while True:
			rmsd = self.calc_rmsd()
			if rmsd < threshold: 
				write = self.write_pdb(it)
				print self.compute_dihedrals(self.chain)
				return "success", rmsd, it
			if it == max_it: return "MAX_IT", rmsd, it
			for i in range(len(self.chain)-2):# for almost each edge find best rotation angle
				N_coors = self.chain[i]
				Ca_coors = self.chain[i+1]
				rot_axis = self.get_rotation_axis(N_coors, Ca_coors)
				Ca_N_vec_norm = self.normalize(Ca_coors - N_coors)
				b = 0.0; c = 0.0 # S = (a,b,c)
				for j in range(1, 4):
					initial_coors = self.chain[-j]; target_coors = self.target[-j]
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
					self.chain[j] = self.rotate_points(self.chain[j], pntO, rot_angle, Ca_N_vec_norm)
			if plot=='y' or plot=='yes' or plot=='on':
				[X,Y,Z]=reshape(self.chain)
				line1.set_ydata(Y)
				fig.canvas.draw()
				it += 1

def make_random_chain_and_run(chain_len,chain_base,target_coors,plot,number_of_loops=10):
	global max_iterations
	success = 0
	av_dist = 3.1 # computed from CG structure
	for n in range(number_of_loops):
		chain = chain_base
		current_point = chain_base[-1]
		for j in range(chain_len): # actual chain will be 3+chain_len+3 b/c of target and anchor
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
				chain[k] = CCD.rotate_points(chain[k], center_point, phi, vec_PQ_norm)
		init = CCD(chain, target_coors)
		check, RMSD, it = init.run(plot,max_it=max_iterations)
		if check == 'success': success += 1
		print("{0}\t{1}\t\t{2}\t{3}\t{4}".format(n, str(round(success*100/float(n+1),2))+'%', check, str(RMSD), it))

def run_CCD():
	chain_base =   [np.array([7.81,19.43,28.83]),np.array([10.17,16.96,28.01]),np.array([9.84,16.63,25.46])]
	target_coors = [np.array([27.95,28.18,31.82]),np.array([29.68,28.79,34.75]),np.array([30.23,29.99,37.54])]
	chain_length = 20
	number_of_loops = 1
	plot='on'
	make_random_chain_and_run(chain_length,chain_base,target_coors,plot,number_of_loops)

if __name__ == "__main__":
	run_CCD()
