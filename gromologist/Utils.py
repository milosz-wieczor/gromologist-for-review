import gromologist as gml
from typing import Optional, Iterable

# TODO make top always optional between str/path and gml.Top


def generate_dftb3_aa(top: "gml.Top", pdb: "gml.Pdb", selection: str):
    special_atoms = {'N': -0.43, 'H': 0.35, 'HN': 0.35, 'C': 0.55, 'O': -0.47}
    atoms = top.get_atoms(selection)
    print("The following atoms were found:")
    for at in atoms:
        print(str(at))
    out = input("Proceed? (y/n)\n")
    if out.strip().lower() != 'y':
        return
    top.parameters.add_dummy_def('LA')
    mols = list(set(at.molname for at in atoms))
    for mol in mols:
        molecule = top.get_molecule(mol)
        current_atoms = [at for at in molecule.atoms if at in atoms]
        atom_indices = [at.num for at in current_atoms]
        current_bonds = molecule.get_subsection('bonds').entries_bonded
        for bond in current_bonds:
            if bond.atom_numbers[0] in atom_indices and bond.atom_numbers[1] in atom_indices:
                bond.interaction_type = '5'
                bond.params_state_a = []
        for atom in current_atoms:
            if atom.atomname not in special_atoms.keys():
                atom.charge = 0.0
            else:
                atom.charge = special_atoms[atom.atomname]
        cas = [at for at in current_atoms if at.atomname == 'CA']
        cbs = [at for at in current_atoms if at.atomname == 'CB']
        assert len(cas) == len(cbs)
        for ca, cb in zip(cas, cbs):
            molecule.add_vs2(ca.num, cb.num, 0.72, 'LIN', 'LA')
            molecule.add_constraint(ca.num, cb.num, 0.155)
        # TODO add vs2 to PDB for each chain that is affected
        cas_all, cbs_all = [at for at in atoms if at.atomname == 'CA'], [at for at in atoms if at.atomname == 'CB']
        for ca, cb in zip(cas_all, cbs_all):
            mol = top.get_molecule(ca.molname)
            pdb_num_ca = mol._match_pdb_to_top(ca.num)[0]
            pdb_num_cb = mol._match_pdb_to_top(cb.num)[0]
            last_atom = mol._match_pdb_to_top(len(mol.atoms))[0]
            pdb.add_vs2(ca.resid, pdb_num_ca, pdb_num_cb, 'LIN', fraction=0.72, serial=last_atom)