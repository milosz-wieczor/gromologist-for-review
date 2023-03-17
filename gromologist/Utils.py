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

def parse_frcmod(filename):
    content = open(filename).readlines()
    atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded = {}, {}, {}, {}, {}, {}
    headers = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB', 'LJED']
    current = None
    for line in content:
        if any([line.strip().startswith(i) for i in headers]):
            current = line.strip()[:4]
            continue
        if current is None or not line.strip() or line.strip().startswith('#'):
            continue
        if current == 'BOND':
            types = tuple(x.strip() for x in line[:5].split('-'))
            vals = tuple(float(x) for x in line[5:].split()[:2])
            bondtypes[types] = [vals[1]/10, vals[0] * 200 * 4.184]
        elif current == 'ANGL':
            types = tuple(x.strip() for x in line[:8].split('-'))
            vals = tuple(float(x) for x in line[8:].split()[:2])
            angletypes[types] = [vals[1], vals[0] * 2 * 4.184]
        elif current == 'MASS':
            types = line.split()[0]
            mass = float(line.split()[1])
            atomtypes[types] = [mass]
        elif current == 'NONB':
            types = line.split()[0]
            rmin = float(line.split()[1])
            eps = float(line.split()[2])
            atomtypes[types].extend([rmin * 0.2 * 2**(-1/6), eps * 4.184])
        elif current == 'LJED':
            types = tuple(line.split()[:2])
            vals = tuple(line.split()[2:])
            assert vals[0] == vals[2] and vals[1] == vals[3]
            nonbonded[types] = [float(vals[0]) * 0.2 * 2**(-1/6), float(vals[1]) * 4.184]
        elif current == 'DIHE':
            types = tuple(x.strip() for x in line[:12].split('-'))
            vals = tuple(float(x) for x in line[12:].split()[:4])
            entry = [vals[2], 4.184 * vals[1] / vals[0], int((vals[3]**2)**0.5)]
            if types in dihedraltypes.keys():
                dihedraltypes[types].extend(entry)
            else:
                dihedraltypes[types] = entry
        elif current == 'IMPR':
            types = tuple(x.strip() for x in line[:12].split('-'))
            vals = tuple(float(x) for x in line[12:].split()[:3])
            entry = [vals[1], 4.184 * vals[0], int((vals[2]**2)**0.5)]
            impropertypes[types] = entry
    assert(all([len(val) == 3 for val in atomtypes.values()]))
    return atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded


def load_frcmod(top: "gml.Top", filename: str):
    atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded = parse_frcmod(filename)
    params = top.parameters
    for at in atomtypes.keys():
        params.add_atomtype(at, *atomtypes[at], action_default='r')
    for b in bondtypes.keys():
        params.add_bonded_param(b, bondtypes[b], 1, action_default='r')
    for a in angletypes.keys():
        params.add_bonded_param(a, angletypes[a], 1, action_default='r')
    for d in dihedraltypes.keys():
        # TODO add wildcards at the end?
        params.add_bonded_param(d, dihedraltypes[d], 9, action_default='r')
    for i in impropertypes.keys():
        params.add_bonded_param(i, impropertypes[i], 4, action_default='r')
    for n in nonbonded.keys():
        try:
            params.add_nbfix(*n, new_sigma=nonbonded[n][0], new_epsilon=nonbonded[n][1])
        except KeyError:
            print(f"Skipping NBFIX {n} as at least one of the types is not defined; if you want to keep it, "
                  "create/load the type and run this command again.")
