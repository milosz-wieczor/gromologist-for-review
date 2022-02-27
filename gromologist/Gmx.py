from subprocess import call
import os


def gmx_command(gmx_exe, command='grompp', answer=False, pass_values=None, quiet=False, **params):
    if pass_values is not None:
        pv = "echo " + ' '.join([str(x) for x in pass_values]) + ' | '
    else:
        pv = ''
    qui = ' &> /dev/null' if quiet else ''
    call_command = f'{pv}{gmx_exe} {command} ' + ' '.join([f'-{k} {v}' for k, v in params.items()]) + qui
    result = call(call_command, shell=True)
    if answer:
        return result


def gen_mdp(fname, runtype='md', **extra_args):
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
    content = [[float(x) for x in l.split()[1:]] for l in open(fname) if not l.startswith(('#', '@'))]
    return content


def calc_gmx_energy(struct, topfile, gmx='', quiet=False):
    # TODO read all terms + hide output
    if not gmx:
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
    if not gmx:
        gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
    gen_mdp('rerun.mdp')
    gmx_command(gmx, 'grompp', quiet=quiet,  f='rerun.mdp', p=topfile, c=struct, o='rerun', maxwarn=5)
    gmx_command(gmx, 'mdrun', quiet=quiet, deffnm='rerun', rerun=struct)
    gmx_command(gmx, 'energy', quiet=quiet, pass_values=[9], f='rerun')
    out = read_xvg('energy.xvg')
    for filename in ['rerun.mdp', 'mdout.mdp', 'rerun.tpr', 'rerun.trr', 'rerun.edr', 'rerun.log', 'energy.xvg']:
        os.remove(filename)
    return out[-1][-1]
