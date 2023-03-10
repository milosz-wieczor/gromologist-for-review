import gromologist as gml
from typing import Optional, Iterable, Union

# TODO make top always optional between str/path and gml.Top


def generate_dftb3_aa(top: "gml.Top", selection: str, pdb: Optional[Union[str, "gml.Pdb"]] = None):
    """
    Prepares a DFT3B-compatible topology and structure, setting up amino acids
    for QM/MM calculations (as defined by the selection)
    :param top: gml.Top, a Topology object
    :param selection: str, a selection defining the residues to be modified
    :param pdb: gml.Pdb, a Pdb object (optional, alternatively can be an attribute of top)
    :return: None
    """
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
        if pdb is not None and top.pdb is None:
            top.add_pdb(pdb)

        for ca, cb in zip(cas_all, cbs_all):
            mol = top.get_molecule(ca.molname)
            for pdb_num_ca, last_atom in zip(mol._match_pdb_to_top(ca.num), mol._match_pdb_to_top(len(mol.atoms))):
                resid = top.pdb.atoms[pdb_num_ca].resnum
                chain = top.pdb.atoms[pdb_num_ca].chain
                top.pdb.add_vs2(resid, 'CA', 'CB', 'LIN', fraction=0.72, serial=last_atom, chain=chain)


# TODO move REST2 preparation here