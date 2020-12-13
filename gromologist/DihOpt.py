import gromologist as gml
from copy import deepcopy
import os
from subprocess import call

try:
    from scipy.optimize import dual_annealing, shgo
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    pass
try:
    from multiprocessing import Pool
except ImportError:
    pass


class DihOpt:
    def __init__(self, top, qm_ref, cutoff=None, qm_unit='kj/mol', processes=1, tmpi=False):
        self.orig_top = top
        self.orig_vals = None
        self.optimizable_params = self.orig_top.parameters.get_opt_dih()
        self.traj = None
        self.frame = None
        self.opt_vals = None
        self.energy_term = 0
        # TODO so far only works with Gaussian outputs
        qm_unit = qm_unit.lower()
        if qm_unit not in ['kj/mol', 'kcal/mol', 'hartree']:
            raise RuntimeError("Unit can be kj/mol, kcal/mol or hartree")
        # we're converting everything to kJ/mol
        scaling = 1.0 if qm_unit == 'kj/mol' else 4.184 if 'kcal/mol' else 2625.5
        self.cutoff = cutoff * scaling if cutoff is not None else 10000000
        if isinstance(qm_ref, str) and qm_ref.endswith('.dat'):  # a file with energy values
            self.qm_ref = [float(x.split()[0])*scaling for x in open(qm_ref)]
        elif isinstance(qm_ref, str) and qm_ref.endswith('.log'):  # a .log file from Gaussian
            self.qm_ref = self.read_gaussian_energies(qm_ref)
        elif isinstance(qm_ref, list) and type(qm_ref[0]) in [int, float]:  # a list with energy values
            self.qm_ref = [x*scaling for x in qm_ref]
        elif isinstance(qm_ref, list) and type(qm_ref[0]) == str:  # a list of Gaussian logs
            self.qm_ref = [e for eners in [self.read_gaussian_energies(ref_file) for ref_file in qm_ref] for e in eners]
        else:
            raise RuntimeError("For qm_ref, pass a Gaussian .log file, a list of Gaussian .logs, a .dat file with "
                               "per-frame energies, or a list with energy values.")
        self.qm_ref = [x - min(self.qm_ref) for x in self.qm_ref]
        self.mod_tops = [deepcopy(self.orig_top) for _ in range(processes)]
        self.processes = processes
        self.run_count = 0
        self.tmpi = tmpi
        self.gmx = self.find_gmx()
        self.gen_dirs()
        self.write_mdp()

    def calc_energy(self, ff_values=None, sys=0, cleanup=True):
        ff_values = self.orig_top.parameters.get_opt_dih() if ff_values is None else ff_values
        self.mod_tops[sys].parameters.set_opt_dih(ff_values)
        self.mod_tops[sys].save_top('opt{}/mod.top'.format(sys))
        self.gmx_grompp(sys)
        self.gmx_mdrun(sys)
        energy = self.gmx_energy(sys)
        if cleanup:
            self.cleanup(sys)
        return energy

    def calc_rmse(self, ff_values, sys=0, cleanup=True):
        energies = self.calc_energy(ff_values, sys, cleanup)
        rmse = [(x-y)**2 for x, y in zip(energies, self.qm_ref) if y < self.cutoff]
        rmse = (sum(rmse)/len(rmse))**0.5
        call('echo {r} >> opt{s}/rmse'.format(s=sys, r=rmse), shell=True)
        return rmse

    def get_bounds(self, sys=0):
        vals = self.mod_tops[sys].parameters.get_opt_dih()
        return [(0, 360), (0, 15)] * (len(vals) // 2)

    def optimize(self):
        # check for imports first
        self.run_count += 1
        self.orig_vals = self.calc_energy()
        if self.processes == 1:
            self.opt_vals = dual_annealing(self.calc_rmse, bounds=self.get_bounds(0), maxiter=200,
                                           x0=self.mod_tops[0].parameters.get_opt_dih())
            # opt_vals = shgo(self.calc_rmse, args=(0,), bounds=self.get_bounds(0), options={'maxev': 1000})
            self.mod_tops[0].save_top('opt{}_{}'.format(self.run_count, self.mod_tops[0].top))
            print(self.opt_vals)
        else:
            p = Pool()
            self.opt_vals = p.map(mappable, [(i, self) for i in range(self.processes)])

    def gmx_grompp(self, sys=0):
        call('{gmx} grompp -f run.mdp -p opt{sys}/mod.top -c {frame} -o opt{sys}/mod.tpr -po opt{sys}/mdout.mdp '
             '-maxwarn 5 >> gmp.log 2>&1'.format(gmx=self.gmx, sys=sys, frame=self.frame), shell=True)

    def gmx_mdrun(self, sys=0):
        if self.tmpi:
            tmpi = ' -ntmpi 1 '
        else:
            tmpi = ''
        call('{gmx} mdrun -s opt{sys}/mod.tpr -o opt{sys}/mod.trr -x opt{sys}/mod.xtc -e opt{sys}/mod.edr '
             '-g opt{sys}/mod.log -c opt{sys}/mod.gro -v -ntomp 1 {tmpi} -rerun {traj} >> mdr.log '
             '2>&1'.format(gmx=self.gmx, sys=sys, tmpi=tmpi, traj=self.traj), shell=True)

    def gmx_energy(self, sys=0):
        if not self.energy_term:
            self.get_energy_term(sys)
        call('echo "{et} 0" | {gmx} energy -f opt{sys}/mod.edr -o opt{sys}/mod.xvg >> edr.log '
             '2>&1'.format(gmx=self.gmx, sys=sys, et=self.energy_term), shell=True)
        vals = np.loadtxt('opt{sys}/mod.xvg'.format(sys=sys), comments=['#', '@'])[:, 1]
        return vals - np.min(vals)

    def get_energy_term(self, sys):
        call('echo "0" | {gmx} energy -f opt{sys}/mod.edr > legend.log 2>&1'.format(gmx=self.gmx, sys=sys), shell=True)
        self.energy_term = os.popen('cat legend.log | grep Potential | awk \'{for (i=1;i<=NF;i++){if ($i=="Potential")'
                                    ' {print $(i-1)}}}\'').read().strip()

    def cleanup(self, sys=0):
        call('rm opt{s}/mdout.mdp opt{s}/mod.xvg opt{s}/mod.edr opt{s}/mod.tpr opt{s}/mod.log'.format(s=sys),
             shell=True)

    @staticmethod
    def write_mdp():
        mdp_defaults = {"nstenergy": 1, "nstlog": 0, "nstcalcenergy": 1, "nstxout-compressed": 0, "pbc": "xyz",
                        "coulombtype": "Cut-off", "vdw-type": "Cut-off", "rlist": 2.0, "rcoulomb": 2.0, "rvdw": 2.0}
        mdp = '\n'.join(["{} = {}".format(param, value) for param, value in mdp_defaults.items()])
        with open('run.mdp', 'w') as outfile:
            outfile.write(mdp)

    def gen_dirs(self):
        for i in range(self.processes):
            try:
                os.mkdir('opt{}'.format(i))
            except FileExistsError:
                pass

    @staticmethod
    def find_gmx():
        """
        Attempts to find Gromacs internal files to fall back to
        when default .itp files are included using the
        #include statement
        :return: str, path to share/gromacs/top directory
        """
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
        return gmx

    def read_gaussian_energies(self, log, traj='opt_traj.pdb', structfile='struct.pdb'):
        log_contents = [x.strip().split() for x in open(log)]
        natoms = [int(x[1]) for x in log_contents if x and x[0] == "NAtoms="][0]
        energies = []
        read_geo = False
        geoms = []
        energy_tmp = None
        for n, line in enumerate(log_contents):
            if len(line) > 2 and line[0] == 'SCF' and line[1] == 'Done:':
                energy_tmp = float(line[4])
            if len(line) > 3 and line[1] == 'Stationary' and line[2] == 'point' and line[3] == 'found.':
                if energy_tmp is not None:
                    energies.append(energy_tmp)
                else:
                    raise RuntimeError("Keyword 'Stationary point found' before first 'SCF done' in Gaussian log, check"
                                       "your QM run")
                read_geo = True
            if read_geo and len(line) > 1 and line[0] == 'Standard' and line[1] == 'orientation:':
                struct = np.array([[float(line[x]) for x in [3, 4, 5]] for line in log_contents[n+5:n+5+natoms]])
                geoms.append(struct)
                read_geo = False
        opt_traj = gml.Traj(top=self.orig_top, array=geoms)
        opt_traj.save_traj_as_pdb(traj)
        opt_traj.save_pdb(structfile)
        self.traj = traj
        self.frame = structfile
        return [2625.5 * x for x in energies]


def mappable(arg):
    sys, dihopt = arg
    return dual_annealing(dihopt.calc_rmse, args=(sys, True), bounds=dihopt.get_bounds(0), maxiter=20, seed=sys,
                          x0=dihopt.mod_tops[0].parameters.get_opt_dih())