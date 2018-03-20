import MDAnalysis as mda
import MDAnalysis.analysis.distances as mdad
import sys

def write_merged(protein_pdb, linker_pdb, n):

	#protein_pdb = "/home/dillion/data/reflectin/add_linker/structures/small_mem_protein-inserted_frame-738.pdb"
	#linker_pdb = sys.argv[1]
	
	prot1 = mda.Universe(protein_pdb) ; pl = len(prot1.atoms.residues)
	prot2 = mda.Universe(protein_pdb)
	link = mda.Universe(linker_pdb) ; ll = len(link.atoms.residues)
	
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
	prot2.atoms.translate([30,0,0])
	for i in range(pl):
		prot2.atoms.residues[i].resid += (pl+ll)
		sel3+=str(prot2.atoms.residues[i].resid)+' '
	
	u = mda.Merge(prot1.atoms,link.atoms,prot2.atoms)
	u.atoms.write("merged"+str(n)+".pdb")
	
	#protein = u.select_atoms(sel1)
	#linker = u.select_atoms(sel2)
	
	#dist = mdad.distance_array(protein.positions,linker.positions)
