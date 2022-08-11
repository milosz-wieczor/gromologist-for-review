from concurrent.futures import ProcessPoolExecutor
import argparse
import os
import sys
from shutil import copy
import numpy as np
from scipy.integrate import simps
from scipy.optimize import minimize, fmin
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from subprocess import call


class CrooksPool:
    def __init__(self, struct, top, xtc0, xtc1, nruns=100, length=500, gmx=None, maxwarn=1, convergence=False,
                 alias='free', debug=False, offset=0, mpi=True, nmax=None, init_eq=50, temperature=300, stride=20,
                 plumed='', top2='', tmpi=True, mutate='', pmxff='', wt_only=False, mini=False, struct2='',
                 plumed2='', hbond=False, dt=2, weights=None, random=False, **kwargs):
        """
        A pool of workers that will initialize all simulations
        and collect the results afterwards
        :param struct: str, a valid .gro or .pdb structure file for the system
        :param top: str, a .top topology file for the system
        :param xtc0: str, a trajectory capturing the equilibrium ensemble for state A (can also be .dcd or any other
        format accepted by MDTraj)
        :param xtc1: str, a trajectory capturing the equilibrium ensemble for state B
        :param nruns: int, number of A->B runs (overall, twice this many runs will be performed)
        :param length: int, length of a single trajectory (in ps)
        :param gmx: str, Gromacs executable (gmx by default)
        :param maxwarn: int, max number of warnings accepted by grompp
        :param convergence: bool, whether to do convergence analysis
        :param alias: str, a base name for all output files
        :param debug: bool, whether to print detailed convergence info
        :param offset: int, start numbering workers from this number (useful when more jobs are added post-factum)
        :param kwargs: dict of str:str pairs, will be converted into additional .mdp options
        """
        self.struct = struct
        self.struct2 = struct2 if struct2 else self.struct
        self.top = top
        self.top2 = top2 if top2 else self.top
        self.name = alias
        self.bandwidth = 0
        self.debug = debug
        self.convergence = convergence
        self.xtc = {0: xtc0, 1: xtc1}
        self.random = random
        self.nruns = nruns
        self.stride = stride
        self.mpi = mpi
        self.mini = mini
        self.plumed = plumed
        self.plumed2 = plumed2 if plumed2 else self.plumed
        if self.mpi:
            from mpi4py.futures import MPIPoolExecutor
        self.init_length = init_eq
        self.sim_length = length
        self.extra_args = kwargs
        self.lincs = 'all-bonds' if not hbond else 'h-bonds'
        self.maxwarn = maxwarn
        self.offset = offset
        self.temperature = temperature
        self.mutate = mutate
        self.pmxff = pmxff
        self.wt_only = wt_only
        self.dt = dt / 1000
        self.nmax = nmax if nmax is not None else self.nruns
        self.init_nst = self.init_length * int(1 / self.dt)
        self.nst = self.sim_length * int(1 / self.dt)
        self.gmx = 'gmx' if gmx is None else gmx
        self.tmpi = tmpi
        self.workers = [Crooks(self, n, lam) for n in range(offset, self.nruns + offset) for lam in [0, 1]]
        self.analysis_only = [Crooks(self, n, lam) for n in range(offset) for lam in [0, 1]]
        if weights is None:
            self.weights = [np.ones(self.nruns) / self.nruns, np.ones(self.nruns) / self.nruns]
        else:
            self.weights = [np.loadtxt(weights.replace('X', '0')), np.loadtxt(weights.replace('X', '1'))]
            for w in self.weights:
                w /= np.sum(w)

    def run(self):
        """
        The main routine: extracts initial frames,
        runs grompp and mdrun in parallel
        :return: None
        """
        self.drop_frames()
        if self.mutate:
            print('mutating structures...'.format())
            with ProcessPoolExecutor() as executor:
                executor.map(self.apply_pmx, [(self.mutate, self.wt_only, w.id, w.initlambda, self.pmxff)
                                              for w in self.workers])
        with ProcessPoolExecutor() as executor:
            executor.map(self.run_worker, [(w, 'mini', self.stride, self.temperature, self.init_nst, self.nst,
                                            self.lincs, self.extra_args, self.mini, self.dt) for w in self.workers])

    @staticmethod
    def run_worker(data):
        worker, runtype, stride, temperature, init_nst, nst, lincs, extra_args, mini, dt = data
        worker.prep_runs('eq', stride, temperature, init_nst, nst, lincs, dt, extra_args)
        if mini:
            worker.prep_runs('mini', stride, temperature, init_nst, nst, lincs, dt, extra_args)
            # We can minimize if required
            worker.log('running energy minimization...'.format())
            worker.grompp_me('mini')
            worker.mdrun_me('mini')
        else:
            worker.cp('run{id}_l{lam}/frame{id}_l{lam}.gro'.format(id=worker.id, lam=worker.initlambda),
                      'mini{id}_l{lam}.gro'.format(id=worker.id, lam=worker.initlambda))
        # First there's equilibration:
        worker.prep_runs('eq', stride, temperature, init_nst, nst, lincs, dt, extra_args)
        worker.log('running equilibration...')
        worker.grompp_me('eq')
        worker.mdrun_me('eq')
        # Then the actual alchemistry:
        worker.prep_runs('prod', stride, temperature, init_nst, nst, lincs, dt, extra_args)
        worker.log('running grompp for all systems...')
        worker.grompp_me('prod')
        worker.mdrun_me('prod')
        worker.log('running mdrun for all systems...')

    def analyze(self):
        """
        The analysis routine: performed after the run,
        reads all data files, processes them and
        prints out the results and/or convergence data
        :return:
        """
        print('analyzing results...')
        all_workers = self.analysis_only + self.workers
        for worker in all_workers:
            worker.analyze_me()
        works_0 = np.array([crk.work for crk in all_workers if crk.initlambda == 0]).reshape(-1, 1)
        # TODO sign
        works_1 = np.array([-crk.work for crk in all_workers if crk.initlambda == 1]).reshape(-1, 1)
        np.savetxt('work_0.dat', works_0, fmt='%10.5f')
        np.savetxt('work_1.dat', works_1, fmt='%10.5f')
        min_val = np.min(np.vstack([works_0, works_1]))
        max_val = np.max(np.vstack([works_0, works_1]))
        diff = max_val - min_val
        min_range = min_val - 0.1 * diff
        max_range = max_val + 0.1 * diff
        grid = np.linspace([min_range], [max_range], 200).reshape(-1, 1)
        np.savetxt('prob_0.dat', np.hstack([grid, dens0]), fmt='%10.5f')
        np.savetxt('prob_1.dat', np.hstack([grid, dens1]), fmt='%10.5f')
        result_bar_pmx = self.solve_bar_pmx(works_0, works_1, *self.weights)
        np.savetxt('normalized_overlap_{}.dat'.format(self.name), self.bhattacharyya(dens0, dens1, grid),
                   fmt='%10.5f')
        if self.convergence:
            self.analyze_convergence(works_0, works_1, grid)

    @staticmethod
    def bhattacharyya(d0, d1, grid):
        """
        Calculates the Bhattacharyya (similarity) coefficient between distributions d0 and d1
        :param d0: np.array, rho_0, density of A->B work values
        :param d1: np.array, rho_1, density of B->A -work values
        :param grid: np.array, the corresponding X values
        :return: np.array of shape [1,1]
        """
        return simps(np.sqrt(d0.reshape(-1) * d1.reshape(-1)), grid.reshape(-1)).reshape(-1, 1)

    def analyze_convergence(self, w0, w1, grid):
        """
        Performs a series of bootstraps for a number of
        subsample sizes, then plots the resulting means/
        standard deviations to a file
        :param w0: np.array, work values for A->B
        :param w1: np.array, -work values for A->B
        :param grid: np.array, the X values, as in the original plot
        :return: None
        """
        samples = np.linspace(4, self.nruns, 10).astype(int)
        conv = {ns: [] for ns in samples}
        if self.debug:
            try:
                os.mkdir('bootstrap')
            except FileExistsError:
                pass
        for nsamp in samples:
            print('bootstrap for {} samples...'.format(nsamp))
            with ProcessPoolExecutor(max_workers=10) as executor:
                results = np.array(list(executor.map(self.boot, [(w0, w1, nsamp, grid, i)
                                                                 for i in range(10)]))).reshape(-1)
            conv[nsamp] = [np.mean(results), np.std(results)]

    def boot(self, params):
        """
        Performs a single bootstrap resampling,
        parallelized through ProcessPoolExecutor
        :param params: a 5-element tuple of input parameters:
        A->B work, B->A -work, number of subsamples to draw, original X values, random seed
        :return: float, the KDE solution for the intersection
        """
        w0, w1, nsamp, grid, n = params
        np.random.seed(n)
        w0_sampled = np.random.choice(w0.reshape(-1), nsamp, replace=True).reshape(-1, 1)
        w1_sampled = np.random.choice(w1.reshape(-1), nsamp, replace=True).reshape(-1, 1)
        dens0_sampled = self.get_opt_kde(w0_sampled, grid, self.bandwidth)
        dens1_sampled = self.get_opt_kde(w1_sampled, grid, self.bandwidth)
        result_kde = self.solve_kde(dens0_sampled, dens1_sampled, grid)
        result_cgi = self.solve_bar_pmx(w0_sampled, w1_sampled)
        if self.debug:
            np.savetxt('bootstrap/conv_{}_{}_{}_{}.dat'.format(nsamp, n, str(result_kde), str(result_cgi)),
                       np.hstack([grid, dens0_sampled, dens1_sampled]), fmt='%10.5f')
        return result_kde

    @staticmethod
    def solve_kde(d0, d1, grid):
        """
        Finds the intersection point based
        on KDE densities
        :param d0: np.array,
        :param d1: np.array,
        :param grid: np.array,
        :return: float, abscissa of intersection
        """
        diff = d0 - d1
        avg = 0.5 * (np.argmax(d0) + np.argmax(d1))
        solution = [x for x in np.where(diff[1:] * diff[:-1] < 0)[0]]
        arg = solution[np.argmin([np.abs(x - avg) for x in solution])]
        result = round(0.5 * (grid.reshape(-1)[arg] + grid.reshape(-1)[arg]), 3)
        return result

    @staticmethod
    def solve_cgi(w0, w1):
        """
        Finds the intersection point based on the CGI
        estimator (intersection of two Gaussians)
        :param w0: np.array, work values for A->B
        :param w1: np.array, -work values for B->A
        :return: float, abscissa of intersection
        """
        w0, w1 = w0.reshape(-1), w1.reshape(-1)
        m0, s0 = np.mean(w0), np.std(w0)
        m1, s1 = np.mean(w1), np.std(w1)
        x1 = m0 / s0 ** 2 - m1 / s1 ** 2
        den = 1 / s0 ** 2 - 1 / s1 ** 2
        sq = (m0 - m1) ** 2 / (s0 ** 2 * s1 ** 2) + 2 * den * np.log(s1 / s0)
        mid = (m0 + m1) / 2
        r1, r2 = (x1 - np.sqrt(sq)) / den, (x1 + np.sqrt(sq)) / den
        return round(r1, 3) if np.abs(r1 - mid) < np.abs(r2 - mid) else round(r2, 3)

    def solve_bar(self, w0, w1, guess):
        kbt = 2.494 * self.temperature / 300

        def fun_ssq(f, work0, work1, beta):
            ni = len(work0)
            nj = len(work1)
            logf = np.log(ni / nj)
            term1 = np.sum(1 / (1 + np.exp(logf + beta * work0 - beta * f)))
            term2 = np.sum(1 / (1 + np.exp(-logf - beta * work1 + beta * f)))
            return (term1 - term2) ** 2

        opt = minimize(fun_ssq, guess, args=(w0, w1, 1 / kbt), options={'disp': False}, tol=0.0000001)
        return opt.x[0]

    def solve_bar_pmx(self, w0, w1):
        return BAR(w0, w1, self.temperature).dg

    def get_opt_kde(self, y, x, bandwidth=None):
        """
        Performs KDE (kernel density estimation) with a Gaussian
        kernel, and optimizes the bandwidth using the LeaveOneOut
        cross-validation if requested
        :param y: np.array, work values
        :param x: np.array, X values for the density
        :param bandwidth: if None, optimize; if a float, use this value
        :return: np.array, the estimated density
        """
        if bandwidth is None:
            bandwidths = 10 ** np.linspace([-1], [1], 20).reshape(-1)
            grid = GridSearchCV(estimator=KernelDensity(kernel='gaussian'), param_grid={'bandwidth': bandwidths},
                                cv=LeaveOneOut())
            grid.fit(y.reshape(-1, 1))
            kde = KernelDensity(kernel='gaussian', **grid.best_params_).fit(y.reshape(-1, 1))
            self.bandwidth = grid.best_params_['bandwidth']
        else:
            kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(y.reshape(-1, 1))
        return np.exp(kde.score_samples(x)).reshape(-1, 1)

    def plot_results(self, w0, w1, p0, p1, result_kde, result_cgi, result_bar, result_bar_pmx):
        """
        Plots the two densities, their intersection
        and all individual work values at the bottom;
        also writes the precise solution to a file
        (this is admittedly a bit confusing)
        :param w0: np.array, work values for A->B
        :param w1: np.array, -work values for B->A
        :param p0: np.array, X:Y two columns for the first density
        :param p1: np.array, X:Y two columns for the second density
        :param result_kde: float, KDE result
        :param result_cgi: float, CGI result
        :param result_bar: float, BAR result
        :return: None
        """
        fig, axes = plt.subplots(2, 1)
        axes[0].plot(*p0.T, c='C0')
        axes[0].plot(*p1.T, c='C1')
        with open('dG_{}.dat'.format(self.name), 'w') as outfile:
            axes[0].axvline(result_kde, ls='--', c='k', lw=3)
            outfile.write("KDE: " + str(result_kde) + '\n')
            outfile.write("CGI: " + str(result_cgi) + '\n')
            outfile.write("BAR: " + str(result_bar) + '\n')
            outfile.write("BAR-PMX: " + str(result_bar_pmx) + '\n')
        axes[0].fill_between(*p0.T, color='C0', alpha=0.25)
        axes[0].fill_between(*p1.T, color='C1', alpha=0.25)
        axes[0].plot(*p1.T, c='C1')
        axes[1].scatter(w0, np.zeros(len(w0)), c='C0')
        axes[1].scatter(w1, np.zeros(len(w1)), c='C1')
        axes[0].set_xlim(np.min(p0[:, 0]), np.max(p0[:, 0]))
        axes[1].set_xlim(np.min(p0[:, 0]), np.max(p0[:, 0]))
        datarange = np.max(p0[:, 0]) - np.min(p0[:, 0])
        axes[1].set_ylim(-datarange / 2, datarange / 2)
        axes[1].get_yaxis().set_visible(False)
        axes[0].get_xaxis().set_visible(False)
        axes[1].set_aspect(0.05)
        axes[1].spines['top'].set_visible(False)
        plt.subplots_adjust(hspace=-0.47)
        plt.savefig('{}.svg'.format(self.name))
        plt.close()

    def drop_frames(self):
        for i in [0, 1]:
            for worker in [w for w in self.workers if w.initlambda == i]:
                worker.mkdir()
                plu = self.plumed if i == 0 else self.plumed2
                if plu:
                    if plu.startswith('/'):
                        plumed = plu
                    else:
                        plumed = os.getcwd() + os.sep + plu
                    worker.cp(plumed, 'plumed.dat')
                    worker.plumed = True
            if not all([os.path.isfile('run{}_l{}/frame{}_l{}.gro'.format(num, i, num, i))
                        for num in range(self.offset, self.offset + self.nruns)]):
                struct = self.struct2 if i == 1 else self.struct
                traj = md.load(self.xtc[i], top=struct)  # TODO optional trjconv if mdtraj is missing
                nfr = len(traj)
                if not self.random:
                    frames = np.linspace([0], [nfr - 0.00001], self.nmax).astype(int)[
                             self.offset:self.offset + self.nruns]
                else:
                    frames = np.random.choice(nfr, size=self.nmax, replace=False)
                for num, fr in enumerate(frames, self.offset):
                    curr = traj[fr]
                    curr.save_gro('run{}_l{}/frame{}_l{}.gro'.format(num, i, num, i))
            else:
                print("frames for lambda {} already present, skipping generating new ones".format(i))

    @staticmethod
    def apply_pmx(data):
        mutate, wt_only, num, lam, pmxff = data
        resp = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS',
                'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN',
                'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}
        resnums, targets = [], []
        for mut in mutate.split('/'):
            if mut[1] in '1234567890':
                if lam == 0:
                    target = resp[mut[-1]]
                elif lam == 1 and not wt_only:
                    target = resp[mut[0]]
                else:
                    target = None
                resnum = int(mut[1:-1])
            else:
                if lam == 0:
                    target = mut[-2:]
                elif lam == 1 and not wt_only:
                    target = mut[:2]
                else:
                    target = None
                resnum = int(mut[2:-2])
            resnums.append(resnum)
            targets.append(target)
        import pmx
        model = pmx.Model('run{}_l{}/frame{}_l{}.gro'.format(num, lam, num, lam))
        for resnum, target in zip(resnums, targets):  # TODO OP1, OP2 to O1P, O2P
            if target is not None:
                model = pmx.mutate(model, resnum, target, pmxff)
        model.writeGRO('run{}_l{}/frame{}_l{}.gro'.format(num, lam, num, lam))

    @staticmethod
    def grompp(data):
        crook, runtype = data
        crook.grompp_me(runtype)

    @staticmethod
    def mdrun(data):
        crook, runtype = data
        crook.mdrun_me(runtype)

    @staticmethod
    def read_results(crook):
        crook.analyze_me()


class Crooks:
    def __init__(self, master, num, initlambda):
        """
        A single worker instance
        :param master: a CrooksPool instance
        :param num: ID of the worker
        :param initlambda: initial alchemical lambda (0 or 1)
        """
        self.id = num
        self.master = master
        self.initlambda = initlambda
        self.plumed = False
        self.work = 0

    def mkdir(self):
        """
        Attempts to create a directory if not there yet
        :return: None
        """
        try:
            os.mkdir('run{}_l{}'.format(self.id, self.initlambda))
        except FileExistsError:
            pass

    def cp(self, filename, target_filename):
        copy(filename, 'run{}_l{}/{}'.format(self.id, self.initlambda, target_filename))

    def log(self, command):
        with open('run{}_l{}/crooks.log'.format(self.id, self.initlambda), 'a') as out:
            out.write(command.strip() + '\n')

    def prep_runs(self, runtype, stride, temperature, init_nst, nst, lincs, dt, extra_args):
        """
        Prepares the run by writing customized .mdp files
        :return: None
        """
        print('preparing runs in mode "{}"...'.format(runtype))
        mdp_defaults = {"nstcomm": 100, "nstenergy": 5000, "nstlog": 5000, "nstcalcenergy": 100,
                        "nstxout-compressed": int(1 / dt) * stride, "compressed-x-grps": "System",
                        "compressed-x-precision": 2500, "free-energy": "yes", "sc-alpha": 0.3, "sc-coul": "yes",
                        "sc-sigma": 0.25, "dt": dt, "constraints": lincs, "coulombtype": "PME",
                        "tcoupl": "v-rescale", "ref-t": temperature, "tau-t": 0, "ref-p": 1.0,
                        "rlist": 1.2, "rcoulomb": 1.2, "vdw-type": "Cut-off", "rvdw_switch": 0.8, "rvdw": 1.2,
                        "compressibility": "4.5e-5", "tau-p": 1.0, "tc-grps": "System"}
        mini_defaults = {"integrator": "steep", "free-energy": "yes", "sc-alpha": 0.3, "sc-coul": "yes",
                         "sc-sigma": 0.25, "nsteps": 1000, "emtol": 200, "emstep": 0.001, "nstlist": 10,
                         "pbc": "xyz", "coulombtype": "PME", "vdw-type": "Cut-off"}
        if runtype == 'eq':
            extra_defaults = {"gen-vel": "yes", "gen-temp": temperature, "pcoupl": "Berendsen",
                              "separate-dhdl-file": "no", "nsteps": init_nst, "nstxout": 10000, "nstvout": 10000}
        else:
            extra_defaults = {"gen-vel": "no", "pcoupl": "Parrinello-Rahman", "separate-dhdl-file": "yes",
                              "nsteps": nst, "nstxout": 100000, "nstvout": 100000}
        mdp_defaults.update(extra_defaults)
        mdp_defaults.update(extra_args)
        default = mini_defaults if runtype == 'mini' else mdp_defaults
        pref = 1 if self.initlambda == 0 else -1
        dl = 0 if runtype in ['eq', 'mini'] else pref / nst
        default.update({"init-lambda": self.initlambda, "delta-lambda": dl})
        mdp = '\n'.join(["{} = {}".format(param, value) for param, value in default.items()])
        name = "mini" if runtype == 'mini' else 'md'
        with open('run{}_l{}/{}.mdp'.format(self.id, self.initlambda, name), 'w') as outfile:
            outfile.write(mdp.format(nst=nst, init=self.initlambda, growth=pref / nst))

    def grompp_me(self, runtype):
        """
        Checks for the .tpr file, runs grompp if not present
        :return: None
        """
        tpr = 'eq' if runtype == 'eq' else 'mini' if runtype == 'mini' else 'dyn'
        mdp = 'mini' if runtype == 'mini' else 'md'
        trr = '' if runtype in ['mini', 'eq'] else ' -t eq{}_l{}.trr -time {} '.format(self.id, self.initlambda,
                                                                                       self.master.init_length)
        os.chdir('run{}_l{}'.format(self.id, self.initlambda))
        master_top = self.master.top if int(self.initlambda) == 0 else self.master.top2
        frame = 'mini{}_l{}.gro'.format(self.id, self.initlambda) if runtype == 'eq' \
            else 'frame{}_l{}.gro'.format(self.id, self.initlambda)
        if '{}{}_l{}.tpr'.format(tpr, self.id, self.initlambda) not in os.listdir('.'):
            top = master_top if master_top.startswith('/') else '../' + master_top
            call('{gmx} grompp -f {mdp}.mdp -p {top} -c {frame} -o {tpr}{n}_l{l}.tpr {trr}'
                 '-maxwarn {mw} >> gmp.log 2>&1'.format(gmx=self.master.gmx, l=self.initlambda, n=self.id, trr=trr,
                                                        top=top, mw=self.master.maxwarn, tpr=tpr, mdp=mdp,
                                                        frame=frame), shell=True)
            os.remove('mdout.mdp')
        os.chdir('..')

    def mdrun_me(self, runtype):
        """
        Checks for the final .gro file, runs mdrun if not present
        :return: None
        """
        tpr = 'eq' if runtype == 'eq' else 'mini' if runtype == 'mini' else 'dyn'
        if self.plumed and runtype == 'eq':
            plu = ' -plumed plumed.dat '
        else:
            plu = ''
        if self.master.tmpi:
            tmpi = ' -ntmpi 1 '
        else:
            tmpi = ''
        os.chdir('run{}_l{}'.format(self.id, self.initlambda))
        if '{}{}_l{}.gro'.format(tpr, self.id, self.initlambda) not in os.listdir('.'):
            call('{gmx} mdrun -deffnm {tpr}{n}_l{l} -v -ntomp 1 {tmpi} {plu} -cpi >> mdr.log '
                 '2>&1'.format(gmx=self.master.gmx, l=self.initlambda, n=self.id, tpr=tpr, plu=plu, tmpi=tmpi),
                 shell=True)
        os.chdir('..')

    def analyze_me(self):
        """
        Reads the .xvg file, removes potential duplicates (if the
        job was restarted), and integrates dH/dl to calculate work
        :return: None
        """
        dhdl = np.loadtxt('run{n}_l{l}/dyn{n}_l{l}.xvg'.format(n=self.id, l=self.initlambda), comments=['#', '@'])
        if not int(dhdl[-1, 0]) == int(self.master.sim_length):
            raise RuntimeError("In file run{}_l{}/dyn{}_l{}.xvg last line reads {}, less than the requested total "
                               "fime of {} ps".format(self.id, self.initlambda, self.id, self.initlambda, dhdl[-1],
                                                      self.master.sim_length))
        dhdl_dict = {i: j for i, j in dhdl}
        dhdl = np.array([dhdl_dict[i] for i in sorted(list(dhdl_dict.keys()))])  # avoid duplicates in case of restarts
        endlambda = 1 if self.initlambda == 0 else 0
        self.work = simps(dhdl, np.linspace([self.initlambda], [endlambda], len(dhdl)).reshape(-1))


class BAR:
    def __init__(self, wf, wr, T):
        self.kb = 0.00831447215  # kJ/(K*mol)
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nf = len(wf)
        self.nr = len(wr)
        self.beta = 1. / (self.kb * self.T)
        self.M = self.kb * self.T * np.log(float(self.nf) / float(self.nr))
        self.dg = self.calc_dg(self.wf, self.wr, self.T)

    def calc_dg(self, wf, wr, T):
        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1. / (self.kb * T)
        M = self.kb * T * np.log(nf / nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1 / (1 + np.exp(beta * (M + v - x)))
            sr = 0
            for v in wr:
                sr += 1 / (1 + np.exp(-beta * (M + v - x)))
            r = sf - sr
            return r ** 2

        avA = np.average(wf)
        avB = np.average(wr)
        x0 = (avA + avB) / 2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)
        return float(dg)



    args, extra = parse_args()
    pool = CrooksPool(struct=args.gro, top=args.top, xtc0=args.xtc, xtc1=args.xtc2, nruns=args.njobs, length=args.time,
                      gmx=args.gmx, maxwarn=args.maxwarn, convergence=args.conv, alias=args.alias, debug=args.debug,
                      offset=args.offset, mpi=not args.nompi, nmax=args.nmax, init_eq=args.inittime,
                      temperature=args.temperature, stride=args.stride, plumed=args.plumed, top2=args.top2,
                      tmpi=not args.notmpi, mutate=args.mutate, pmxff=args.pmxff, wt_only=args.wt_only, mini=args.mini,
                      struct2=args.gro2, plumed2=args.plumed2, hbond=args.hbond, dt=args.timestep, weights=args.weights,
                      random=args.random, **extra)
    if not args.norun:
        import mdtraj as md

        pool.run()
    if not args.noan:
        pool.analyze()
