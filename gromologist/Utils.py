import os
from glob import glob
from shutil import copy2
from itertools import permutations
from copy import deepcopy

import gromologist as gml
from typing import Optional, Iterable, Union


# TODO make top always optional between str/path and gml.Top


def generate_dftb3_aa(top: Union[str, "gml.Top"], selection: str, pdb: Optional[Union[str, "gml.Pdb"]] = None):
    """
    Prepares a DFT3B-compatible topology and structure, setting up amino acids
    for QM/MM calculations (as defined by the selection)
    :param top: gml.Top, a Topology object
    :param selection: str, a selection defining the residues to be modified
    :param pdb: gml.Pdb, a Pdb object (optional, alternatively can be an attribute of top)
    :return: None
    """
    top = gml.obj_or_str(top=top)
    pdb = gml.obj_or_str(pdb=pdb)
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
    pdb = gml.obj_or_str(pdb=pdb)
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


# TODO move REST2 preparation here

def parse_frcmod(filename: str) -> (dict, dict, dict, dict, dict, dict, dict):
    """
    Parses either an frcmod file with extra parameters, or a parm file with
    core FF parameters (both have slightly different formats for opening/closing sections)
    :param filename: str, name of the file to be read
    :return: tuple of dict, FF parameres in their respective formats
    """
    dat = True if filename.endswith('dat') else False
    content = open(filename).readlines()
    if any(['MOD4' in l and 'AC' in l for l in content]):
        raise RuntimeError("LJ type A/C not supported, terminating")
    content = content[1:] if dat else content
    atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded, cmaptypes = {}, {}, {}, {}, {}, {}, {}
    cmapresol, cmapres, cmapvals, cmapread = None, [], [], False
    headers = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'HBON', 'NONB', 'LJED', 'CMAP']
    identical_nonbonded = {}
    iterator = 0
    current = headers[iterator] if dat else None
    for nl, line in enumerate(content):
        if not dat:
            if any([line.strip().startswith(i) for i in headers]):
                current = line.strip()[:4]
                continue
            if current is None or not line.strip() or line.strip().startswith('#'):
                continue
        else:
            if not line.strip() and iterator < len(headers) - 3:
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
            else:
                if len(line.split()) < 3:
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
        elif current == 'CMAP':
            types = tuple('C N CT C N'.split()) # in case the format ever changes in the future
            if line.startswith('%FLAG'):
                if line.split()[1] == "CMAP_RESLIST":
                    cmapres = content[nl+1].split()
                elif line.split()[1] == "CMAP_RESOLUTION":
                    cmapresol = line.split()[2]
                elif line.split()[1] == "CMAP_PARAMETER":
                    cmapread = True
            elif cmapread:
                if not line.strip():
                    cmapread = False
                else:
                    cmapvals.extend(line.strip().split())
            if len(cmapvals) > 0 and len(cmapvals) == int(cmapresol) ** 2:
                cmapvals = [str(round(4.184 * float(i), 10)) for i in cmapvals]
                for res in cmapres:
                    cmaptypes[(types, res)] = [cmapresol, cmapvals]
                cmapresol, cmapres, cmapvals, cmapread = None, [], [], False

    # assert (all([len(val) == 3 for val in atomtypes.values()]))
    atomtypes = {key: val for key, val in atomtypes.items() if len(val) == 3}
    non_atomtypes = [key for key, val in atomtypes.items() if len(val) != 3]
    if non_atomtypes:
        print(f"skipping atomtypes {non_atomtypes}, missing LJ parameters")
    return atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded, cmaptypes


def load_frcmod(top: Union[str, "gml.Top"], filename: str, return_cmap: bool = False, pro_atoms: Optional[dict] = None)\
        -> Optional[dict]:
    """
    Loads an .frcmod file into an existing topology. Can be also launched as
    gml.Top().load_frcmod(...)
    :param top: str or gml.Top, existing gmx topology
    :param filename: str, name of the frcmod file to load
    :param return_cmap: bool, if set to True will return cmaptypes
    :return: None or dict, depending on return_cmap
    """
    top = gml.obj_or_str(top)
    atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes, nonbonded, cmaptypes = parse_frcmod(filename)
    for k in cmaptypes.keys():
        try:
            cc = [adata[1] for adata in pro_atoms[k[1]] if adata[0] == "C"][0]
            na = [adata[1] for adata in pro_atoms[k[1]] if adata[0] == "N"][0]
            ca = [adata[1] for adata in pro_atoms[k[1]] if adata[0] == "CA"][0]
        except KeyError:
            cc, na, ca = "C", "N", "CA"
        cmaptypes[k].append(f'{cc} {na} {ca} {cc} {na}'.split())
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
    for c in cmaptypes.keys():
        params.add_bonded_param(cmaptypes[c][-1], [c[1], cmaptypes[c][:-1]], 1, action_default='a')
    if return_cmap:
        return cmaptypes
    else:
        return None


def read_lib(lib: str) -> (dict, dict, dict):
    """
    Reads a .lib file with residue definitions
    :param lib: str, name of the .lib file
    :return: tuple of dict, dictionary with atoms, bonds, and inter-residue connectors
    """
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
    atoms = {k: v for k, v in atoms.items() if 'BOX' not in k}
    bonds = {k: v for k, v in bonds.items() if 'BOX' not in k}
    connector = {k: v for k, v in connector.items() if 'BOX' not in k}
    return atoms, bonds, connector


def write_rtp(atoms: dict, bonds: dict, connector: dict, outfile: str = "new.rtp", ff='amber',
              impropers: Optional[dict] = None, cmap: Optional[dict] = None):
    """
    Writes an .rtp file given all per-residue dictionary with topology, extra dihedrals/impropers, CMAP etc.
    :param atoms: dict of tuple, atom names/types and their ID/charge
    :param bonds: dict of tuples, 1-based indices of pairs defining bonds
    :param connector: dict of tuples, connecting atoms from -1 or +1 residue in polymers
    :param outfile: str, to which file output should be written
    :param ff: str, 'amber' or 'charmm' to set correct [ bondedtypes ] (default interaction types)
    :param impropers: dict of tuples, atom names that should be involved in improper dihedrals
    :param cmap: dict of tuples, atom names that should be involved in the cmap correction
    :return: None
    """
    if ff.lower() not in ['amber', 'charmm']:
        raise RuntimeError("Only Amber and CHARMM are currently supported")
    btypes = '11941310' if ff == 'amber' else '15921310'
    print(f"Setting [ bondedtypes ] in {outfile} file for the {ff}-type force field, please make sure this is right")
    with open(outfile, 'w') as out:
        out.write(f"[ bondedtypes ]\n{' '.join(btypes)}\n\n")
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
            if impropers is not None:
                if res in impropers.keys():
                    out.write(f" [ impropers ]\n")
                    for imp in impropers[res]:
                        out.write(f" {imp[0]:5s} {imp[1]:5s} {imp[2]:5s} {imp[3]:5s}\n")
            if cmap is not None:
                if res in cmap.keys():
                    out.write(f" [ cmap ]\n")
                    for cmp in cmap[res]:
                        out.write(f" {cmp[0]:5s} {cmp[1]:5s} {cmp[2]:5s} {cmp[3]:5s} {cmp[4]:5s}\n")
            out.write("\n\n")


def dict_filter(indict: dict, restype: str) -> dict:
    """
    Filters dictionaries based on whether they belong to DNA, RNA or protein
    :param indict: dict, dictionary with keys that are residue names
    :param restype: str, type of residues (DNA, RNA, or anything else for protein)
    :return: dict, the filtered dictionary
    """
    nucres = gml.Pdb.nucl_map.keys()
    if restype == "DNA":
        return {k: v for k, v in indict.items() if k in nucres and 'D' in nucres}
    elif restype == "RNA":
        return {k: v for k, v in indict.items() if k in nucres and 'D' not in nucres}
    else:
        return {k: v for k, v in indict.items() if k not in nucres}


def fix_rtp(rtp_dict: dict, impr: bool = False, rna: bool = False) -> dict:
    """
    Makes Gromacs-specific changes in .rtp data, e.g. adjusts atom names, terminal atom naming,
    hydrogen numbering, RNA residue names, improper type order etc.
    :param rtp_dict: dict, dictionary with .rtp data (atoms or impropers)
    :param impr: bool, whether the entry contains improper dihedral data
    :param rna: bool, whether the entry represents RNA residues
    :return: dict, the dictionary with necessary modifications
    """
    to_copy = {}
    for res in rtp_dict.keys():
        if not rna:
            if "ILE" in res:
                for n, ent in enumerate(rtp_dict[res]):
                    if ent[0] == "CD1":
                        rtp_dict[res][n] = ('CD', *ent[1:])
                    elif ent[0].startswith("HD1"):
                        rtp_dict[res][n] = (ent[0].replace("HD1", "HD"), *ent[1:])
            if res.startswith("C") and res[1:] in gml.Pdb.prot_map.keys():
                if not impr:
                    for n, ent in enumerate(rtp_dict[res]):
                        if ent[0] == "O":
                            rtp_dict[res][n] = ('OC1', *ent[1:])
                        elif ent[0] == "OXT":
                            rtp_dict[res][n] = ('OC2', *ent[1:])
                else:
                    for n, ent in enumerate(rtp_dict[res]):
                        for m, atname in enumerate(ent):
                            if atname == "O":
                                rtp_dict[res][n][m] = 'OC2'
                            elif atname == "OXT":
                                rtp_dict[res][n][m] = 'OC1'
            for mid in ["A", "B", "G", "D", "E", "G1"]:
                if f"H{mid}3" in [e[0] for e in rtp_dict[res]] and f"H{mid}1" not in [e[0] for e in rtp_dict[res]]:
                    for n, ent in enumerate(rtp_dict[res]):
                        if ent[0] == f"H{mid}3":
                            rtp_dict[res][n] = (f'H{mid}1', *ent[1:])
            if impr:
                for n, ent in enumerate(rtp_dict[res]):
                    for m, atname in enumerate(ent):
                        if rtp_dict[res][n][m] == '-M':
                            rtp_dict[res][n][m] = '-C'
                        elif rtp_dict[res][n][m] == '+M':
                            rtp_dict[res][n][m] = '+N'
        if rna:
            for res in rtp_dict.keys():
                if res[0] in 'CAGU':
                    to_copy[res] = 'R' + res
                if impr:
                    continue
                for n, ent in enumerate(rtp_dict[res]):
                    try:
                        _ = ent[0]
                    except TypeError:
                        continue
                    if ent[0] == "OP1":
                        rtp_dict[res][n] = ('O1P', *ent[1:])
                    elif ent[0] == "OP2":
                        rtp_dict[res][n] = ('O2P', *ent[1:])
                    elif ent[0] == "H5''":
                        rtp_dict[res][n] = ("H5'2", *ent[1:])
                    elif ent[0] == "H5'":
                        rtp_dict[res][n] = ("H5'1", *ent[1:])
                    elif ent[0] == "H2'":
                        rtp_dict[res][n] = ("H2'1", *ent[1:])
                    elif ent[0] == "HO2'":
                        rtp_dict[res][n] = ("HO'2", *ent[1:])
                    elif ent[0] == "HO3'":
                        rtp_dict[res][n] = ("H3T", *ent[1:])
                    elif ent[0] == "HO5'":
                        rtp_dict[res][n] = ("H5T", *ent[1:])
    for cop in to_copy.keys():
        rtp_dict[to_copy[cop]] = rtp_dict[cop]
    for clean in to_copy.keys():
        del rtp_dict[clean]
    return rtp_dict


def find_leap_files(ambdir: str, filetype: str, leaprc: str, extras: Optional[list] = None):
    content = [line.strip() for line in open(leaprc)]
    extras = extras if extras is not None else []
    pref = os.path.sep if leaprc.startswith(os.path.sep) else ''
    reldir = pref + os.path.sep.join(leaprc.split(os.path.sep)[:-1]) + pref
    if filetype == 'lib':
        files = [line.split()[-1] for line in content if
                 len(line.split()) >= 2 and "loadoff" in [w.lower() for w in line.split()]] + extras
    elif filetype == 'parm':
        files = [line.split()[-1] for line in content if
                 len(line.split()) >= 2 and "loadamberparams" in [w.lower() for w in line.split()]] + extras
    else:
        raise RuntimeError('select "parm" or "lib"')
    files_in_ambdir = [ambdir + f'/{filetype}/' + f for f in files]
    files_in_reldir = [reldir + '/' + f for f in files]
    print(files_in_reldir)
    found_files = [f for f in files_in_ambdir + files_in_reldir if os.path.exists(f)]
    return found_files


def amber2gmxFF(leaprc: str, outdir: str, amber_dir: Optional[str] = None, base_ff: Optional[str] = None,
                extra_files: Optional[list] = None):
    """
    Reads a .leaprc file and all parameter dependencies from Amber to convert into a Gromacs .ff dir
    Files that should be copied manually: watermodels.dat and tip*itp, .hdb, .tdb and .arn
    :param leaprc: str, a file that sources dependencies from which the .ff will be created
    :param outdir: str, a new .ff directory that will contain the Gromacs-compatible files
    :param amber_dir: str, Abs path to the dir containing Amber prep, parm, lib directories if `leaprc` is a local file
    :param base_ff: str, name of the source .ff directory that will be used for auxilliary files
    :param extra_files: list of str, additional files that should be included (e.g. frcmod.tip3p)
    :return: None
    """
    orig_dir = os.path.sep.join(leaprc.split(os.path.sep)[:-1]) + os.path.sep if os.path.sep in leaprc else ''
    if amber_dir is not None:
        amb = amber_dir
    else:
        amb = f'{orig_dir}../'
    extra_files = extra_files if extra_files is not None else []
    extra_libs = [pth for pth in extra_files if pth.endswith('lib')]
    extra_dats = [pth for pth in extra_files if pth.endswith('dat') or 'frcmod' in pth.split(os.path.sep)[-1]]
    extra_prep = [pth for pth in extra_files if pth.endswith('in') or pth.endswith('prep')]
    libs = find_leap_files(amb, 'lib', leaprc, extra_libs)
    dats = find_leap_files(amb, 'parm', leaprc, extra_dats)
    pro_atoms, pro_bonds, pro_connectors = {}, {}, {}
    dna_atoms, dna_bonds, dna_connectors = {}, {}, {}
    rna_atoms, rna_bonds, rna_connectors = {}, {}, {}
    impropers = {}
    prepsel = "/prep/all*.in"
    for prep in glob(amb + prepsel) + glob(amb + "/prep/nuc*.in") + extra_prep:
        impropers.update(read_prep_impropers(prep))
    print('\n')
    for lib in libs:
        print(f"***** Adding residues from {lib} *****")
        a, b, c = gml.read_lib(lib)
        pro_atoms.update(fix_rtp(dict_filter(a, 'protein')))
        pro_bonds.update(dict_filter(b, 'protein'))
        pro_connectors.update(dict_filter(c, 'protein'))
        dna_atoms.update(dict_filter(a, 'DNA'))
        dna_bonds.update(dict_filter(b, 'DNA'))
        dna_connectors.update(dict_filter(c, 'DNA'))
        rna_atoms.update(fix_rtp(dict_filter(a, 'RNA'), rna=True))
        rna_bonds.update(fix_rtp(dict_filter(b, 'RNA'), rna=True))
        rna_connectors.update(fix_rtp(dict_filter(c, 'RNA'), rna=True))
    print('\n')
    pro_impropers = fix_rtp(dict_filter(impropers, 'protein'), impr=True)
    dna_impropers = dict_filter(impropers, 'DNA')
    rna_impropers = fix_rtp(dict_filter(impropers, 'RNA'), impr=True, rna=True)
    imp_types_all = set()
    for atm, imp in zip([pro_atoms, dna_atoms, rna_atoms], [pro_impropers, dna_impropers, rna_impropers]):
        for key in atm.keys():
            typedict = {ent[0]: ent[1] for ent in atm[key]}
            try:
                imp_names = imp[key]
            except:
                continue
            for impn in imp_names:
                try:
                    imp_types = tuple(typedict[i.lstrip('-+')] for i in impn)
                except:
                    continue
                if imp_types not in imp_types_all and imp_types[::-1] not in imp_types_all:
                    imp_types_all.add(imp_types)
    new_top = gml.Top(amber=True)
    cmaptypes, rtp_cmap = {}, {}
    print('\n')
    for dat in dats:
        print(f"***** Adding parameters from {dat} *****")
        cmaptypes.update(load_frcmod(new_top, dat, return_cmap=True, pro_atoms=pro_atoms))
    print('\n')
    for k in cmaptypes.keys():
        rtp_cmap[k[1]] = [f'-C N CA C +N'.split()]
    new_top = reorder_amber_impropers(new_top, imp_types_all)
    for typical_dummy in ['EP', 'DR']:
        new_top.parameters.add_dummy_def(typical_dummy)
    outdir = outdir + '.ff' if not outdir.endswith('.ff') else outdir
    atomtypes = read_addAtomTypes(leaprc)
    os.mkdir(outdir)
    os.chdir(outdir)
    all_types_rtp = set()
    for atoms_dict in [pro_atoms, dna_atoms, rna_atoms]:
        all_types_rtp.update({at[1] for res in atoms_dict.keys() for at in atoms_dict[res]})
    if all_types_rtp.difference(new_top.defined_atomtypes):
        print(f"WARNING: Type(s) {all_types_rtp.difference(new_top.defined_atomtypes)} defined in .rtp but not in "
              f"ffparams.itp, consider removing the entry from .rtp or adding it to ffparams.itp")
    for atype in new_top.parameters.atomtypes.entries_param:
        if atype.modifiers[0] == '0' and atype.types[0] in atomtypes.keys():
            atype.modifiers[0] = str(atomtypes[atype.types[0]][0])
        if atype.modifiers[1] == '0' and atype.types[0] in atomtypes.keys():
            atype.modifiers[1] = str(atomtypes[atype.types[0]][1])
    new_top.parameters.sort_dihedrals()
    new_top.save_top('forcefield.itp', split=True)
    gml.write_rtp(pro_atoms, pro_bonds, pro_connectors, 'aminoacids.rtp', impropers=pro_impropers, cmap=rtp_cmap)
    if dna_atoms:
        gml.write_rtp(dna_atoms, dna_bonds, dna_connectors, 'dna.rtp', impropers=dna_impropers)
    if rna_atoms:
        gml.write_rtp(rna_atoms, rna_bonds, rna_connectors, 'rna.rtp', impropers=rna_impropers)
    gmx_dir = new_top.gromacs_dir
    if not gmx_dir:
        print("Gromacs directory not found. Please move the additional files (.hdb, solvent .itp etc.) manually")
    else:
        if base_ff is None:
            base_ff = gmx_dir + '/amber99.ff'
        else:
            if not base_ff.startswith('/'):
                base_ff = gmx_dir + f'/{base_ff}'
        wildcards_to_copy = ['tip*', 'spc*', '*r2b', '*tdb', '*arn', '*hdb']
        files_to_copy = []
        for wildcard in wildcards_to_copy:
            files_to_copy.extend(glob(f'{base_ff}/{wildcard}'))
        print(f"Copying accessory files (.hdb, .atp, solvent .itp etc.) from {base_ff}\n")
        for file_to_copy in files_to_copy:
            copy2(file_to_copy, '.')
    with open('atomtypes.atp', 'w') as outfile:
        for atype in new_top.parameters.atomtypes.entries_param:
            outfile.write(f"{atype.types[0]:12s} {float(atype.modifiers[1]):12.5f}\n")
    new_top.find_undefined_types()
    # TODO set up watermodels.dat


def reorder_amber_impropers(new_top: "gml.Top", rtp_dihedrals: set) -> "gml.Top":
    """
    Modifying improper dihedral order, empirically checked against GMX FFs
    when errors come up
    :param new_top: gml.Top, a topology to process
    :return: gml.Top
    """
    try:
        _ = new_top.parameters.dihedraltypes
    except KeyError:
        pass
    else:
        dih_sub = new_top.parameters.dihedraltypes
        defined_imps = [dih.types for dih in dih_sub.entries_param if str(dih.interaction_type) in '24']
        rtp_imps = list(rtp_dihedrals)
        reorders = {}
        for rtp_imp in rtp_imps:
            if rtp_imp not in defined_imps and rtp_imp[::-1] not in defined_imps:
                pers = permutations(rtp_imp)
                match = False
                for per in pers:
                    if match:
                        break
                    for ref in defined_imps:
                        if dih_match(per, ref):
                            if 'X' not in ref:
                                match = per
                            break
                if match:
                    reorders[match] = rtp_imp
        for types_to_dup in reorders.keys():
            imp = [e for e in dih_sub.get_entries_by_types(*types_to_dup) if str(e.interaction_type) in '24'][0]
            imp.types = reorders[types_to_dup]
        # previous attempt at listing problematic dihedrals, but this gets out of hand too fast
        # new_top.parameters.dihedraltypes.reorder_improper(('CB', 'CT', 'C*', 'CW'), '1203')
        # new_top.parameters.dihedraltypes.reorder_improper(('CT', 'CW', 'CC', 'NB'), '0213')
        # new_top.parameters.dihedraltypes.reorder_improper(('CB', 'N2', 'CA', 'NC'), '0321')
        # new_top.parameters.dihedraltypes.reorder_improper(('CB', 'C5', 'N*', 'CT'), '3201')
        # new_top.parameters.dihedraltypes.reorder_improper(('CB', 'CP', 'N*', 'CT'), '3201')
        # new_top.parameters.dihedraltypes.reorder_improper(('C', 'C4', 'N*', 'CT'), '3201')
        # new_top.parameters.dihedraltypes.reorder_improper(('C', 'CS', 'N*', 'CT'), '3201')
        # new_top.parameters.dihedraltypes.reorder_improper(('C4', 'N2', 'CA', 'NC'), '1203')
        # new_top.parameters.dihedraltypes.reorder_improper(('N2', 'NA', 'CA', 'NC'), '1320')
    return new_top


def dih_match(query, ref):
    if all([q == r for q, r in zip(query, ref) if r != "X"]) or all([q == r for q, r in zip(query, ref[::-1]) if r != "X"]):
        return True
    else:
        return False

def read_addAtomTypes(textfile: str) -> dict:
    text = open(textfile).readlines()
    reading, brack = False, 0
    element_properties = gml.guess_element_properties()
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
    for k, v in types.items():
        if v and v in element_properties.keys():
            types[k] = element_properties[v]
        else:
            types[k] = (0, 0.0)
    return types

def guess_element_properties():
    element_properties = {
        "H": (1, 1.008), "O": (8, 15.999), "C": (6, 12.011), "N": (7, 14.007), "S": (16, 32.06), "P": (15, 30.974),
        "F": (9, 18.998), "Cl": (17, 35.45), "Br": (35, 79.904), "I": (53, 126.90), "Mg": (12, 24.305),
        "Ca": (20, 40.078)
    }
    return element_properties

def read_prep_impropers(prepfile: str) -> dict:
    """
    Reads improper dihedrals from a specified file in the leap/prep directory
    :param prepfile: str, input file for reading
    :return: dict of lists, each dict entry corresponds to a residue, the list contains 4-lists of atomtypes
    """
    impropers = {}
    current = None
    reading = False
    content = [line for line in open(prepfile)]
    for line in content:
        if len(line.split()) > 2:
            if line.strip().split()[1] == "INT":
                current = line.strip().split()[0]
        if line.strip() == "IMPROPER":
            reading = True
        if line.strip() == "DONE":
            reading = False
        if reading and current is not None and line.strip():
            if current not in impropers.keys():
                impropers[current] = []
            if len(line.split()) == 4:
                types = line.strip().split()
                impropers[current].append(types)
    return impropers


def prep_to_rtp(prepfile, outfile = "new.rtp", ff = 'amber'):
    content = [l.strip() for l in open(prepfile)]
    atoms = {}
    bonds = {}
    impropers = {}
    connectors = {}
    current_res = ""
    reading_loop_closing = False
    reading_impropers = False
    for ln in content:
        if len(ln.split()) > 2 and ln.split()[1].upper() == "INT":
            current_res = ln.split()[0]
            atoms[current_res] = []
            bonds[current_res] = []
            impropers[current_res] = []
            connectors[current_res] = []
        elif len(ln.split()) > 9 and ln.split()[3].upper() in "MSB3456E":
            if ln.split()[2] != "DU":
                atoms[current_res].append((ln.split()[1], ln.split()[2], float(ln.split()[-1]), len(atoms[current_res]) + 1))
                if len(atoms[current_res]) > 1:
                    bonds[current_res].append((int(ln.split()[4]) - 3, len(atoms[current_res])))
        elif ln.strip() == "IMPROPER":
            reading_impropers = True
        elif ln.strip().split() and ln.strip().split()[0] == "LOOP":
            reading_loop_closing = True
            reading_impropers = False
        elif ln.strip() == "DONE":
            reading_loop_closing = False
        elif reading_loop_closing and len(ln.split()) > 1:
            at1 = ln.split()[0]
            at2 = ln.split()[1]
            n1 = [ent[3] for ent in atoms[current_res] if ent[0] == at1][0]
            n2 = [ent[3] for ent in atoms[current_res] if ent[0] == at2][0]
            bonds[current_res].append((int(n1), int(n2)))
        elif reading_impropers and len(ln.split()) > 3:
            impropers[current_res].append(ln.split())
    bonds[current_res] = sorted(bonds[current_res], key=lambda x: x[0])
    write_rtp(atoms, bonds, connector=connectors, outfile=outfile, impropers=impropers, ff=ff)



def calc_Coulomb_force(top: Union[str, "gml.Top"], pdb: Union[str, "gml.Pdb"], force_on: str, force_from: str):
    """
    Calculates Coulomb forces exerted by one subset of atoms on each atom in another subset. Does NOT take into accout
    scaled 1-4 interactions.
    :param top: gml.Top or str, topology
    :param pdb: gml.Pdb or str, structure
    :param force_on: selection for which the forces will be calculated
    :param force_from: selection of atoms that are exerting the forces
    :return: Nx3 np.array, where N is the length of the force_on selection
    """
    import numpy as np
    k_coul = 138.9118  # in kJ/mol nm e²
    pdb = gml.obj_or_str(pdb=pdb)
    top = gml.obj_or_str(top=top)
    fon_indices = pdb.get_atom_indices(force_on)
    ffr_indices = pdb.get_atom_indices(force_from)
    topat = top.atoms
    fon_charges = [a.charge for a in [topat[i] for i in fon_indices]]
    ffr_charges = np.array([a.charge for a in [topat[i] for i in ffr_indices]])
    fon_coords = np.array([pdb.get_coords(force_on)])[0]/10
    ffr_coords = np.array([pdb.get_coords(force_from)])[0]/10
    vec_matrix = np.zeros((len(fon_indices), len(ffr_indices), 3))
    for i in range(len(fon_indices)):
        for j in range(len(ffr_indices)):
            vec_matrix[i, j, :] = fon_coords[i] - ffr_coords[j]
    dist_matrix = np.linalg.norm(vec_matrix, axis=2)[..., None]
    geom_matrix = vec_matrix / dist_matrix**3
    force_matrix = np.zeros((len(fon_indices), 3))
    for i in range(len(fon_indices)):
        force_matrix[i, :] = k_coul * fon_charges[i] * np.sum(ffr_charges[..., None] * geom_matrix[i, :, :], axis=0)
    return force_matrix


def calc_LJ_force(top: Union[str, "gml.Top"], pdb: Union[str, "gml.Pdb"], force_on: str, force_from: str):
    """
    Calculates Coulomb forces exerted by one subset of atoms on each atom in another subset. Does NOT take into accout
    scaled 1-4 interactions.
    :param top: gml.Top or str, topology
    :param pdb: gml.Pdb or str, structure
    :param force_on: selection for which the forces will be calculated
    :param force_from: selection of atoms that are exerting the forces
    :return: Nx3 np.array, where N is the length of the force_on selection
    """
    import numpy as np
    pdb = gml.obj_or_str(pdb=pdb)
    top = gml.obj_or_str(top=top)
    fon_indices = pdb.get_atom_indices(force_on)
    ffr_indices = pdb.get_atom_indices(force_from)
    topat = top.atoms
    sigma_matrix = np.array([[top.parameters.sigma_ij(topat[i], topat[j]) for i in ffr_indices] for j in fon_indices])
    epsilon_matrix = np.array(
        [[top.parameters.epsilon_ij(topat[i], topat[j]) for i in ffr_indices] for j in fon_indices])
    fon_coords = np.array([pdb.get_coords(force_on)])[0]/10
    ffr_coords = np.array([pdb.get_coords(force_from)])[0]/10
    vec_matrix = np.zeros((len(fon_indices), len(ffr_indices), 3))
    for i in range(len(fon_indices)):
        for j in range(len(ffr_indices)):
            vec_matrix[i, j, :] = fon_coords[i] - ffr_coords[j]
    dist_matrix = np.linalg.norm(vec_matrix, axis=2)[..., None]
    geom_matrix = vec_matrix / dist_matrix ** 2
    s6_matrix = (sigma_matrix / dist_matrix[:, :, 0]) ** 6
    force_matrix = np.zeros((len(fon_indices), 3))
    for i in range(len(fon_indices)):
        force_matrix[i, :] = 24 * np.sum(epsilon_matrix[i, :][..., None] * geom_matrix[i, :, :] * (2 * s6_matrix**2 - s6_matrix)[i, :][..., None], axis=0)
    return force_matrix
