from subprocess import run, PIPE
import os
import gromologist as gml
from typing import Optional, Iterable


def gmx_command(gmx_exe: str, command: str = 'grompp', answer: bool = False, pass_values: Optional[Iterable] = None,
                quiet: bool = False, **params) -> str:
    """
    Runs the specified gmx command, optionally passing keyworded or stdin arguments
    :param gmx_exe: str, a gmx executable
    :param command: str, the gmx command to launch
    :param answer: bool, whether to read & return the stderr + stdout of the command
    :param pass_values: iterable, optional values to pass to the command (like group selections in gmx trjconv)
    :param quiet: bool, whether to show gmx output
    :param params: dict, for any "-key value" option to be included pass entry formatted as {"key": value}
    :return: str, stdout/stderr output from the command (if answer=True)
    """
    if pass_values is not None:
        pv = (' '.join([str(x) for x in pass_values]) + '\n').encode()
    else:
        pv = None
    qui = ' &> /dev/null' if quiet else ''
    call_command = f'{gmx_exe} {command} ' + ' '.join([f'-{k} {v}' for k, v in params.items()]) + qui
    result = run(call_command.split(), input=pv, stderr=PIPE, stdout=PIPE)
    # result = call(call_command, shell=True)
    if answer:
        return result.stdout.decode() + result.stderr.decode()


def gen_mdp(fname: str, runtype: str = 'md', **extra_args):
    """
    Produces a default .mdp file for the rerun
    :param fname: str, name of the output file
    :param runtype: str, "mini" for minimization or anything else for dynamics
    :param extra_args: dict, optional extra parameter: value pairs (will overwrite defaults); use __ for -
    :return: None
    """
    mdp_defaults = {"integrator": "sd", "nstcomm": 100, "nstenergy": 5000, "nstlog": 5000, "nstcalcenergy": 100,
                    "nstxout-compressed": 5000, "compressed-x-grps": "System",
                    "compressed-x-precision": 2500, "dt": 0.002, "constraints": 'hbonds', "coulombtype": "Cut-off",
                    "ref-t": 300, "tau-t": 1.0, "ref-p": 1.0,
                    "rlist": 1.2, "rcoulomb": 1.2, "vdw-type": "Cut-off", "rvdw_switch": 0.8, "rvdw": 1.2,
                    "ld_seed": -1, "compressibility": "4.5e-5", "tau-p": 1.0,
                    "tc-grps": "System", "gen-vel": "yes", "gen-temp": 300, "pcoupl": "Berendsen",
                    "separate-dhdl-file": "no", "nsteps": 1000, "nstxout": 10000, "nstvout": 10000}
    mini_defaults = {"integrator": "steep", "nsteps": 1000, "emtol": 200, "emstep": 0.001, "nstlist": 10,
                     "pbc": "xyz", "coulombtype": "PME", "vdw-type": "Cut-off"}
    mdp_defaults.update(extra_args)
    for key in list(mdp_defaults.keys()):
        if '__' in key:
            mdp_defaults[key.replace('__', '-')] = mdp_defaults[key]
            del mdp_defaults[key]
    default = mini_defaults if runtype == 'mini' else mdp_defaults
    mdp = '\n'.join([f"{param} = {value}" for param, value in default.items()])
    with open(fname, 'w') as outfile:
        outfile.write(mdp)


def read_xvg(fname: str, cols: Optional[list] = None) -> list:
    """
    Reads an .xvg file into a 2D list
    :param fname: str, .xvg file to read
    :param cols: list of int, columns to select
    :return: list of lists, numeric data from the .xvg file
    """
    content = [[float(x) for x in line.split()[1:]] for line in open(fname) if not line.startswith(('#', '@'))]
    if cols is not None:
        if len(cols) == 1:
            content = [line[cols[0]] for line in content]
        else:
            content = [[line[x] for x in cols] for line in content]
    return content


def get_legend(gmx: str, fname: str) -> dict:
    """
    Performs a dummy run of gmx energy to read the matching between terms and numbers
    :param gmx: str, path to the gmx executable
    :param fname: str, path to the .edr file
    :return: dict, matches between the terms' names and their consecutive numbers
    """
    pp = run([gmx, 'energy', '-f', fname], input=b'0\n', stderr=PIPE, stdout=PIPE)
    output = pp.stderr.decode().split()
    return {output[i + 1].lower(): int(output[i]) for i in range(output.index('1'), len(output), 2)
            if output[i].isnumeric()}


def ndx(struct: gml.Pdb, selections: list, fname: str = 'gml.ndx') -> list:
    """
    Writes a .ndx file with groups g1, g2, ... defined by the
    list of selections passed as input
    :param struct: gml.Pdb, a structure file
    :param selections: list of str, selections compatible with `struct`
    :param fname: str, name of the resulting .ndx file
    :return: list of str, names of the group
    """
    groups = []
    group_names = []
    for n, sel in enumerate(selections, 1):
        groups.append([x + 1 for x in struct.get_atom_indices(sel)])
        group_names.append(f'g{n}')
    with open(fname, 'w') as out:
        for gname, gat in zip(group_names, groups):
            out.write(f'{gname}\n')
            for n, at in enumerate(gat):
                out.write(f'{at:8d}')
                if n % 15 == 14:
                    out.write('\n')
    return group_names


def calc_gmx_energy(struct: str, topfile: str, gmx: str = '', quiet: bool = False, traj: Optional[str] = None,
                    terms: str = 'potential', cleanup: bool = True, group_a: Optional[str] = None,
                    group_b: Optional[str] = None) -> dict:
    """
    Calculates selected energy terms given a structure/topology pair or structure/topology/trajectory set.
    :param struct: str, path to the structure file
    :param topfile: str, path to the topology file
    :param gmx: str, path to the gmx executable (if not found in the $PATH)
    :param quiet: bool, whether to print gmx output to the screen
    :param traj: str, path to the trajectory (optional)
    :param terms: str or list, terms which will be calculated according to gmx energy naming (can also be "all")
    :param cleanup: bool, whether to remove intermediate files (useful for debugging)
    :param group_a: str, selection defining group A to calculate interactions between group A and B
    :param group_b: str, selection defining group B to calculate interactions between group A and B
    :return: dict of lists, one list of per-frame values per each selected term
    """
    if not gmx:
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
    if group_a and group_b:
        group_names = ndx(gml.Pdb(struct), [group_a, group_b])
        gen_mdp('rerun.mdp', energygrps=f"{group_names[0]} {group_names[1]} ")
        gmx_command(gmx, 'grompp', quiet=quiet, f='rerun.mdp', p=topfile, c=struct, o='rerun', maxwarn=5, n='gml.ndx')
    else:
        gen_mdp('rerun.mdp')
        gmx_command(gmx, 'grompp', quiet=quiet, f='rerun.mdp', p=topfile, c=struct, o='rerun', maxwarn=5)
    gmx_command(gmx, 'mdrun', quiet=quiet, deffnm='rerun', rerun=struct if traj is None else traj)
    legend = get_legend(gmx, 'rerun.edr')
    if terms == 'all':
        terms = list(legend.keys())
    if isinstance(terms, str):
        terms = [terms]
    try:
        passv = [legend[i.lower()] for i in terms]
    except KeyError:
        raise RuntimeError(f'Could not process query {terms}; available keywords are: {legend.keys()}')
    gmx_command(gmx, 'energy', quiet=quiet, pass_values=passv, f='rerun')
    out = read_xvg('energy.xvg')
    if cleanup:
        to_remove = ['rerun.mdp', 'mdout.mdp', 'rerun.tpr', 'rerun.trr', 'rerun.edr', 'rerun.log', 'energy.xvg']
        if group_a and group_b:
            to_remove.append('gml.ndx')
        for filename in to_remove:
            os.remove(filename)
    return {term: [o[onum] for o in out] for term, onum in zip(terms, range(len(out[0])))}


def calc_gmx_dhdl(struct: str, topfile: str, traj: str, gmx: str = '', quiet: bool = False,
                  cleanup: bool = True) -> list:
    """
    Calculates selected energy terms given a structure/topology pair or structure/topology/trajectory set.
    :param struct: str, path to the structure file
    :param topfile: str, path to the topology file
    :param gmx: str, path to the gmx executable (if not found in the $PATH)
    :param quiet: bool, whether to print gmx output to the screen
    :param traj: str, path to the trajectory (optional)
    :param cleanup: bool, whether to remove intermediate files (useful for debugging)
    :return: dict of lists, one list of per-frame values per each selected term
    """
    if not gmx:
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
    gen_mdp('rerun.mdp', free__energy="yes", fep__lambdas="0 1", nstdhdl="1", separate__dhdl__file="yes",
            dhdl__derivatives="yes", init__lambda__state="0")
    gmx_command(gmx, 'grompp', quiet=quiet, f='rerun.mdp', p=topfile, c=struct, o='rerun', maxwarn=5)
    gmx_command(gmx, 'mdrun', quiet=quiet, deffnm='rerun', rerun=struct if traj is None else traj)
    out = read_xvg('rerun.xvg', cols=[0])
    if cleanup:
        to_remove = ['rerun.mdp', 'mdout.mdp', 'rerun.tpr', 'rerun.trr', 'rerun.edr', 'rerun.log']
        for filename in to_remove:
            os.remove(filename)
    return out
