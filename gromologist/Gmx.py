from subprocess import run, PIPE
import os


def gmx_command(gmx_exe, command='grompp', answer=False, pass_values=None, quiet=False, **params):
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


def gen_mdp(fname, runtype='md', **extra_args):
    """
    Produces a default .mdp file for the rerun
    :param fname: str, name of the output file
    :param runtype: str, "mini" for minimization or anything else for dynamics
    :param extra_args: dict, optional extra parameter: value pairs (will overwrite defaults)
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
    default = mini_defaults if runtype == 'mini' else mdp_defaults
    mdp = '\n'.join([f"{param} = {value}" for param, value in default.items()])
    with open(fname, 'w') as outfile:
        outfile.write(mdp)


def read_xvg(fname):
    """
    Reads an .xvg file into a 2D list
    :param fname: str, .xvg file to read
    :return: list of lists, numeric data from the .xvg file
    """
    content = [[float(x) for x in line.split()[1:]] for line in open(fname) if not line.startswith(('#', '@'))]
    return content


def get_legend(gmx, fname):
    """
    Performs a dummy run of gmx energy to read the matching between terms and numbers
    :param gmx: str, path to the gmx executable
    :param fname: str, path to the .edr file
    :return: dict, matches between the terms' names and their consecutive numbers
    """
    pp = run([gmx, 'energy', '-f', fname], input=b'0\n', stderr=PIPE, stdout=PIPE)
    output = pp.stderr.decode().split()
    return {output[i+1].lower(): int(output[i]) for i in range(output.index('1'), len(output), 2)
            if output[i].isnumeric()}


def calc_gmx_energy(struct, topfile, gmx='', quiet=False, traj=None, terms='potential', cleanup=True):
    """
    Calculates selected energy terms given a structure/topology pair or structure/topology/trajectory set.
    :param struct: str, path to the structure file
    :param topfile: str, path to the topology file
    :param gmx: str, path to the gmx executable (if not found in the $PATH)
    :param quiet: bool, whether to print gmx output to the screen
    :param traj: str, path to the trajectory (optional)
    :param terms: str or list, terms which will be calculated according to gmx energy naming (can also be "all")
    :param cleanup: bool, whether to remove intermediate files (useful for debugging)
    :return: dict of lists, one list of per-frame values per each selected term
    """
    if not gmx:
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
    gen_mdp('rerun.mdp')
    gmx_command(gmx, 'grompp', quiet=quiet,  f='rerun.mdp', p=topfile, c=struct, o='rerun', maxwarn=5)
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
        for filename in ['rerun.mdp', 'mdout.mdp', 'rerun.tpr', 'rerun.trr', 'rerun.edr', 'rerun.log', 'energy.xvg']:
            os.remove(filename)
    return {term: [o[onum] for o in out] for term, onum in zip(terms, range(len(out[0])))}
