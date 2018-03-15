''' !!! MODIFIED VERSION OF !!! : PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

This class contains a bunch of useful functions that
piggyback on Biopython functions that generate residues.
The Biopython classes are very annoying to use and
manipulate so this just makes it easier.

This file is provided to you under the GNU General Public
License, version 2.0 or later.'''

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

def initialize_res(residue):
	if isinstance( residue, Geo ):
		geo = residue
	else:
		geo=geometry(residue) 
	segID=1
	AA= geo.residue_name
	CA_N_length=geo.CA_N_length
	CA_C_length=geo.CA_C_length
	N_CA_C_angle=geo.N_CA_C_angle
	CA_coord= [0.,0.,0.]
	C_coord= [CA_C_length,0,0]
	N_coord = [CA_N_length*math.cos(N_CA_C_angle*(math.pi/180.0)),CA_N_length*math.sin(N_CA_C_angle*(math.pi/180.0)),0]
	N= Atom("N", N_coord, 0.0 , 1.0, " "," N", 0, "N")
	CA=Atom("CA", CA_coord, 0.0 , 1.0, " "," CA", 0,"C")
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
	cha= Chain('A')
	cha.add(res)
	mod= Model(0)
	mod.add(cha)
	struc= Structure('X')
	struc.add(mod)
	return struc

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

def make_extended_structure(AA_chain):
	geo = geometry(AA_chain[0])
	struc=initialize_res(geo)
	for i in range(1,len(AA_chain)): 
		AA = AA_chain[i]
		geo = geometry(AA)
		add_residue(struc, geo)
	return struc
    
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
 
def make_structure(AA_chain,phi,psi_im1,omega=[]):
	geo = geometry(AA_chain[0])
	struc=initialize_res(geo)
	if len(omega)==0:
		for i in range(1,len(AA_chain)): 
			AA = AA_chain[i]
			add_residue(struc, AA, phi[i-1], psi_im1[i-1])
	else:
		for i in range(1,len(AA_chain)): 
			AA = AA_chain[i]
			add_residue(struc, AA, phi[i-1], psi_im1[i-1], omega[i-1])
	return struc
    
def make_structure_from_geos(geos):
	model_structure=initialize_res(geos[0])
	for i in range(1,len(geos)):
		model_structure=add_residue(model_structure, geos[i])
	return model_structure
