import os

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
    dat = True if filename.endswith('dat') else False
    content = open(filename).readlines()
    if any(['MOD4' in l and 'AC' in l for l in content]):
        raise RuntimeError("LJ type A/C not supported, terminating")
    content = content[1:] if dat else content
    atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded = {}, {}, {}, {}, {}, {}
    headers = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'HBON', 'NONB', 'LJED']
    identical_nonbonded = {}
    iterator = 0
    current = headers[iterator] if dat else None
    for line in content:
        if not dat:
            if any([line.strip().startswith(i) for i in headers]):
                current = line.strip()[:4]
                continue
            if current is None or not line.strip() or line.strip().startswith('#'):
                continue
        else:
            if not line.strip() and iterator < len(headers) - 2:
                iterator += 1
                current = headers[iterator]
                continue
            if line.strip() == "END":
                current = headers[-1]
        if current == 'BOND':
            if dat and '-' not in line[:5]:
                continue
            types = tuple(x.strip() for x in line[:5].split('-'))
            vals = tuple(float(x) for x in line[5:].split()[:2])
            bondtypes[types] = [vals[1] / 10, vals[0] * 200 * 4.184]
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
            if dat:
                try:
                   _ = float(line.split()[1])
                except:
                    if len(line.split()) >= 2 and all([t in atomtypes.keys() for t in line.split()]):
                        identical_nonbonded[line.split()[0]] = line.split()[1:]
                    continue
            tps = line.split()[0]
            rmin = float(line.split()[1])
            eps = float(line.split()[2])
            types = [tps] if tps not in identical_nonbonded.keys() else [tps] + identical_nonbonded[tps]
            for atype in types:
                if atype in atomtypes.keys() and len(atomtypes[atype]) == 1:
                    atomtypes[atype].extend([rmin * 0.2 * 2 ** (-1 / 6), eps * 4.184])
                else:
                    atomtypes[atype] = [0, rmin * 0.2 * 2 ** (-1 / 6), eps * 4.184]
        elif current == 'LJED':
            if dat:
                if not(len(line.split()) > 1 and line.split()[0] in atomtypes.keys() and line.split()[1] in atomtypes.keys()):
                    continue
            types = tuple(line.split()[:2])
            vals = tuple(line.split()[2:])
            assert vals[0] == vals[2] and vals[1] == vals[3]
            nonbonded[types] = [float(vals[0]) * 0.2 * 2 ** (-1 / 6), float(vals[1]) * 4.184]
        elif current == 'DIHE':
            types = tuple(x.strip() for x in line[:12].split('-'))
            vals = tuple(float(x) for x in line[12:].split()[:4])
            entry = [vals[2], 4.184 * vals[1] / vals[0], int((vals[3] ** 2) ** 0.5)]
            if types in dihedraltypes.keys():
                dihedraltypes[types].extend(entry)
            else:
                dihedraltypes[types] = entry
        elif current == 'IMPR':
            types = tuple(x.strip() for x in line[:12].split('-'))
            vals = tuple(float(x) for x in line[12:].split()[:3])
            entry = [vals[1], 4.184 * vals[0], int((vals[2] ** 2) ** 0.5)]
            impropertypes[types] = entry
    #assert (all([len(val) == 3 for val in atomtypes.values()]))
    atomtypes = {k: v for k, v in atomtypes.items() if len(v) == 3} # TODO that's temporary
    natomtypes = {k: v for k, v in atomtypes.items() if len(v) != 3}
    print natomtypes
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


def read_lib(lib: str) -> (dict, dict, dict):
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
                    connector[curr_resname].append(int(content[n + 1].strip()))
    return atoms, bonds, connector


def write_rtp(atoms, bonds, connector, outfile: str = "new.rtp"):
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


def amber2gmxFF(leaprc: str, outdir: str, amber_dir: Optional[str] = None):
    """
    Files that should be copied manually: watermodels.dat and tip*itp, .hdb, .tdb and .arn
    :param leaprc:
    :param outdir:
    :param amber_dir:
    :return:
    """
    content = [line.strip() for line in open(leaprc)]
    orig_dir = os.path.sep.join(leaprc.split(os.path.sep)[:-1]) + os.path.sep if os.path.sep in leaprc else ''
    if amber_dir is not None:
        amb = amber_dir
    else:
        amb = f'{orig_dir}../'
    libs = [amb + '/lib/' + line.split()[1] for line in content if len(line.split()) >= 2 and
            line.split()[0] == "loadOff"]
    dats = [amb + '/parm/' + line.split()[-1] for line in content if len(line.split()) >= 2 and
            line.split()[-2] == "loadamberparams"]
    pro_atoms, pro_bonds, pro_connectors = {}, {}, {}
    dna_atoms, dna_bonds, dna_connectors = {}, {}, {}
    rna_atoms, rna_bonds, rna_connectors = {}, {}, {}
    for lib in libs:
        print(f"Adding residues from {lib}")
        nucres = gml.Pdb.nucl_map.keys()

        def dict_filter(dict, restype):
            if restype == "DNA":
                return {k: v for k, v in dict.items() if k in nucres and 'D' in nucres}
            elif restype == "RNA":
                return {k: v for k, v in dict.items() if k in nucres and 'D' not in nucres}
            else:
                return {k: v for k, v in dict.items() if k not in nucres}

        a, b, c = gml.read_lib(lib)
        pro_atoms.update(dict_filter(a, 'protein'))
        pro_bonds.update(dict_filter(b, 'protein'))
        pro_connectors.update(dict_filter(c, 'protein'))
        dna_atoms.update(dict_filter(a, 'DNA'))
        dna_bonds.update(dict_filter(b, 'DNA'))
        dna_connectors.update(dict_filter(c, 'DNA'))
        rna_atoms.update(dict_filter(a, 'RNA'))
        rna_bonds.update(dict_filter(b, 'RNA'))
        rna_connectors.update(dict_filter(c, 'RNA'))
    new_top = gml.Top(amber=True)
    for dat in dats:
        print(f"Adding parameters from {dat}")
        load_frcmod(new_top, dat)
    os.mkdir(outdir)
    os.chdir(outdir)
    new_top.save_top('forcefield.itp', split=True)
    gml.write_rtp(pro_atoms, pro_bonds, pro_connectors, 'aminoacids.rtp')
    if dna_atoms:
        gml.write_rtp(dna_atoms, dna_bonds, dna_connectors, 'dna.rtp')
    if rna_atoms:
        gml.write_rtp(rna_atoms, rna_bonds, rna_connectors, 'rna.rtp')


def read_addAtomTypes(text: list) -> dict:
    reading, brack = False, 0
    types = {}
    for line in text:
        if line.strip().startswith('addAtomTypes'):
            reading = True
        if reading:
            brack += line.count('{')
            brack -= line.count('}')
        else:
            continue
        if brack == 0:
            reading = False
        data = line.strip().strip('{}').strip().split()
        if len(data) == 3:
            types[data[0].strip('"')] = data[1].strip('"')
    return types
