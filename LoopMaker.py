import os
import numpy as np
import Bio.PDB

class LoopMaker:
        def __init__(self):
		return None

	def write_pdb(self,coors,pdbname):
		outfile = open(pdbname,"w")
		it = 0
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
		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		model = pdb_parser.get_structure("ref",pdbname)[0]
		first_res =  self.model_selectres(model,1)
		last_res =  self.model_selectres(model,self.model_numres(model))
		base_coors = self.res_get_atom_coors(last_res) # C-terminus is base
		target_coors = self.res_get_atom_coors(first_res) # need to connect to N-terminus
		for v in range(3):
			target_coors[v] += protein_shift
		return [last_res, target_coors]

	def write_protein_copy(self,pdbname,protein_shift):
		print "this doesn't work for some reason"
		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		structure = pdb_parser.get_structure("ref",pdbname)
		for model in structure:
			for chain in model:
				for res in chain:
					for atom in res:
						#atom.set_coord = np.array(atom.get_coord()) + np.array(self.protein_shift)
						atom.set_coord = [0,0,0]

		out = Bio.PDB.PDBIO()
		out.set_structure(structure)
		out.save( "structure2.pdb" )
	
	def extract_coors(self,structure):
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
		import PeptideBuilderMod
		for i in range(1,loop_length):
			loop = PeptideBuilderMod.add_residue(loop,"G")
		return loop

	def finish_grafted_chain(self,loop,di):
		import PeptideBuilderMod
		for i in range(len(di)):
			print di[i][0], di[i][1]
			loop = PeptideBuilderMod.add_residue(loop,"G",di[i][0],di[i][1])
		return loop
	
	def model_numres(self,model):
		length = 0
		for residue in model.get_residues():
			length+=1
		return length
	
	def model_selectres(self,model,resnum):
		return model['A'][resnum]
	
	def res_get_atom_coors(self,res):
		v = []
		v.append(np.array([res['N' ].get_coord()[0],res['N' ].get_coord()[1],res['N' ].get_coord()[2]]))
		v.append(np.array([res['CA'].get_coord()[0],res['CA'].get_coord()[1],res['CA'].get_coord()[2]]))
		v.append(np.array([res['C' ].get_coord()[0],res['C' ].get_coord()[1],res['C' ].get_coord()[2]]))
		return v
	
	def get_bb_atoms(self,model):
		atoms = []
		for chain in model:
			for res in chain:
				atoms.append(res['N'])
				atoms.append(res['CA'])
				atoms.append(res['C'])
		return np.array(atoms)
	
	def align_with_pdb(self,base_res):
		# make first residue in loop
		import PeptideBuilderMod
		geo = PeptideBuilderMod.geometry()
		loop = PeptideBuilderMod.initialize_res(geo)
	
		# PeptideBuilder fucks something up, so you need to write to pdb and read in before using SuperImposer
		out = Bio.PDB.PDBIO()
		out.set_structure(loop)
		out.save( "loop1.pdb" )
		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		loop_model = pdb_parser.get_structure("ref","loop1.pdb")[0]
		os.remove( "loop1.pdb" ) # delete unnecessary pdb
		
		# define atoms to align
		base_atoms = np.array([base_res['N'], base_res['CA'], base_res['C']])
		loop_atoms = self.get_bb_atoms(loop_model)
		super_imposer = Bio.PDB.Superimposer()
		super_imposer.set_atoms(base_atoms, loop_atoms)
		super_imposer.apply(loop_model.get_atoms())
		return loop

	def anchor_is_start(self,base_res):
		# biopython is stupid so I have to extract coordinates, save in pdb, and reload to make coordinates into a "structure"
		out = Bio.PDB.PDBIO()
		pdbname = "anchor_res.pdb"
		test = np.array([base_res['N'].get_coord(), base_res['CA'].get_coord(), base_res['C'].get_coord(), base_res['O'].get_coord()])
		self.write_pdb(test,pdbname)
		pdb_parser = Bio.PDB.PDBParser(QUIET = True)
		loop_structure = pdb_parser.get_structure("ref",pdbname)
		os.remove( pdbname ) # delete unnecessary pdb
		return loop_structure

if __name__ == "__main__":
	pdbname = "/home/dillion/data/reflectin/add_linker/structures/small_mem_protein-inserted_frame-738.pdb"
	protein_shift = np.array([30,0,0])
	loop_length = 15
	loopmaker = LoopMaker()
	[base_res, target_coors] = loopmaker.define_initial_cond(pdbname,protein_shift)
	loop_structure = loopmaker.anchor_is_start(base_res)
	loop_structure = loopmaker.finish_straight_chain(loop_structure)
	coors = loopmaker.extract_coors(loop_structure)
	print coors

