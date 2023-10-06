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
            if types in atomtypes.keys():
                atomtypes[types][0] = mass
            else:
                atomtypes[types] = [mass]
        elif current == 'NONB':
            types = line.split()[0]
            rmin = float(line.split()[1])
            eps = float(line.split()[2])
            if types in atomtypes.keys() and len(atomtypes[types]) == 1:
                atomtypes[types].extend([rmin * 0.2 * 2**(-1/6), eps * 4.184])
            else:
                atomtypes[types] = [0, rmin * 0.2 * 2 ** (-1 / 6), eps * 4.184]
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


def lib_to_rtp(lib: str, outfile: str = "new.rtp"):
    curr_resname = None
    atoms = {}
    bonds = {}
    connector = {}
    reading_atoms = False
    reading_bonds = False
    content = [line for line in open(lib) if line.strip()]
    for n, ln in enumerate(content):
        if not ln.startswith('!'):
            if reading_atoms:
                atoms[curr_resname].append((ln.strip().split()[0].strip('"'), ln.strip().split()[1].strip('"'),
                                            float(ln.strip().split()[7]), int(ln.strip().split()[5])))
            elif reading_bonds:
                bonds[curr_resname].append((int(ln.strip().split()[0]), int(ln.strip().split()[1])))
        if ln.startswith('!'):
            if len(ln.strip('!').split()[0].split('.')) < 3:
                continue
            else:
                reading_bonds, reading_atoms = False, False
                if ln.strip('!').split()[0].split('.')[3] == 'atoms':
                    reading_atoms = True
                    curr_resname = ln.strip('!').split()[0].split('.')[1]
                    atoms[curr_resname] = []
                    bonds[curr_resname] = []
                    connector[curr_resname] = []
                elif ln.strip('!').split()[0].split('.')[3] == 'connectivity':
                    reading_bonds = True
                elif ln.strip('!').split()[0].split('.')[3] == 'connect':
                    connector[curr_resname].append(int(content[n+1].strip()))

    with open(outfile, 'w') as out:
        for res in atoms.keys():
            out.write(f"[ {res} ]\n [ atoms ]\n")
            for at in atoms[res]:
                out.write(f"  {at[0]:4s}   {at[1]:4s}          {at[2]:8.5f}     {at[3]:3d}\n")
            out.write(f" [ bonds ]\n")
            for bd in bonds[res]:
                out.write(f"  {atoms[res][bd[0] - 1][0]:4s}   {atoms[res][bd[1] - 1][0]:4s}\n")
            if len(connector[res]) > 0 and connector[res][0] > 0:
                atomlist = [at[0] for at in atoms[res]]
                is_prot = True if 'CA' in atomlist else False
                is_na = True if "O4'" in atomlist else False
                if is_prot:
                    out.write(f"  -C  {atoms[res][connector[res][0] - 1][0]}\n")
                elif is_na:
                    out.write(f"  -O3'  {atoms[res][connector[res][0] - 1][0]}\n")
            out.write("\n\n")
            # TODO no idea how this is encoded in AMBER files


def generate_gaussian_input(pdb: Union["gml.Pdb", str], directive_file: str, outfile: str = 'inp.gau', charge: int = 0,
                            multiplicity: int = 1, group_a: Optional[str] = None, group_b: Optional[str] = None):
    """
    From a .pdb file and an existing Gaussian input, produces a new .gau input
    with correct atom names, coordinates, and possibly fragment assignment
    :param pdb: gml.Pdb or str, the structure object/file containing the desired coordinates
    :param directive_file: str, an existing Gaussian input from which the %- and #-prefixed lines will be taken
    :param outfile: str, a file to which the new input will be written
    :param charge: int, charge of the system (by default 0)
    :param multiplicity: int, multiplicity of the system (by default 1)
    :param group_a: str, selection to define 1st fragment if the counterpoise correction is used
    :param group_b: str, selection to define 2nd fragment if the counterpoise correction is used
    :return: None
    """
    gau_content = [line for line in open(directive_file)]
    pdb = gml.Pdb(pdb) if isinstance(pdb, str) else pdb
    pdb.add_elements()
    with open(outfile, 'w') as outf:
        for line in [ln for ln in gau_content if ln.strip().startswith('%')]:
            outf.write(line)
        for line in [ln for ln in gau_content if ln.strip().startswith('#')]:
            outf.write(line)
        outf.write(f"\ngromologist input to gaussian\n\n{charge} {multiplicity}\n")
        if group_a is None and group_b is None:
            for atom in pdb.atoms:
                outf.write(f" {atom.element}   {atom.x}  {atom.y}  {atom.z}\n")
        elif group_a is not None and group_b is not None:
            for atom in pdb.get_atoms(group_a):
                outf.write(f" {atom.element}(Fragment=1)   {atom.x:8.3f}  {atom.y:8.3f}  {atom.z:8.3f}\n")
            for atom in pdb.get_atoms(group_b):
                outf.write(f" {atom.element}(Fragment=2)   {atom.x:8.3f}  {atom.y:8.3f}  {atom.z:8.3f}\n")
        else:
            raise RuntimeError('Specify either both group_a and group_b, or neither')
        outf.write("\n")
