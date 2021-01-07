# Gromologist

Gromologist is a package designed to facilitate handling, editing and manipulating GROMACS topology files 
(.top and .itp), as well as compatible structures (.pdb and .gro).

## Installation

The latest "official" release can be obtained directly through `pip` by typing `pip install gromologist`.

To get the latest development version, first locally clone the git repo (`git clone https://gitlab.com/KomBioMol/gromologist.git`),
then install the package into Python by typing `pip install .` in the main Gromologist directory.
If you're using Anaconda, the same will work with `/path/to/anaconda/bin/pip`.

## Usage

##### Reading and writing files

`Top` and `Pdb` are the core classes of the library, and are supposed to provide representation
for topology and structure objects, respectively. To initialize them, a path to the file
should be passed to the constructor:

```
>>> from gromologist import Top, Pdb
>>> t = Top('md/topol.top')
>>> p = Pdb('md/conf.pdb')
```

Since all .itp files are by default included into the `Top` object, sometimes it is
necessary to specify a custom path to Gromacs database:

```
>>> t = Top('md/topol.top', gmx_dir='/home/user/gromacs/share/gromacs/top')
```

Alternatively, `Top` can be initialized with both paths, or `Pdb` can be supplied later.
Note that one `Top` instance can only hold one `Pdb` instance at a time.

```
>>> t = Top('md/topol.top', pdb='md/conf.pdb')
>>> t.pdb
PDB file md/conf.pdb with 100 atoms
>>> t.add_pdb('md/other.pdb')
>>> t.pdb
PDB file md/other.pdb with 105 atoms
```

After changes have been made, modified files can be saved:

```
>>> t.save_top('md/new_topol.top')
>>> t.pdb.save_pdb('md/new_conf.pdb')
```

##### File inspection, checks and printing

If `Pdb` is bound to `Top`, a number of diagnostic and fixing options are available,
including name consistency checks:

```
>>> t.check_pdb()
Check passed, all names match
```

as well as renumbering, automatic addition of chains and fixing ordering issues in PDB
using topology ordering as a template:

```
>>> t.pdb.renumber_all()
>>> t.pdb.add_chains()
>>> t.pdb.match_order_by_top_names()
```

Missing atoms in (standard) protein residues can also be identified rapidly:

```
>>> t.pdb.find_missing()
```

With `Pdb.select_atoms()`, selections can be made in a VMD-like manner:

```
>>> t.pdb.select_atoms('name CA and (resname ASP or chain B)')
[5, 60, 72, 88]
```

Several 'convenience' functions exist to list relevant properties of the topology:

```
>>> t.list_molecules()
Protein                      1
Other                        1
>>> protein = t.molecules[0]
>>> protein.print_molecule()
# prints all atoms in the molecule
>>> protein.list_bonds()
# lists bonds, labeling bonded atoms by atom name
>>> protein.list_bonds(by_types=True)
# lists bonds, labeling bonded atoms by atom type
>>> protein.list_bonds(by_params=True)
# lists bonds, adding FF parameter values alongside
```

By analogy, the `.list_bonds()` method can be used to `list_angles`, `list_dihedrals`
and `list_impropers`.


##### Producing lightweight files

If the topology contains too many parts irrelevant to the system at hand,
a leaner version can be produced that lacks unused molecule definitions:

```
>>> t.clear_sections()
```

Similarly, to remove FF parameters that are not used in the molecule definition,
another 'clearing' method can be used:

```
>>> t.clear_ff_params()
```

To save a 'lightweight' .top file with all contents split into separate .itp files, 
use the `split` parameter of `Top.save_top`:

```
>>> t.save_top('md/new_topol.top', split=True) 
```

##### Dealing with unspecified 'define' keywords in topologies

If some FF terms are assumed to be defined elsewhere, e.g. in .mdp files, their values
can be explicitly specified at construction:

```
>>> t = Top('topol.top', define={'POSRES_FC_BB':400.0})
```

On the other hand, some `#define` keywords are included in topology files, and are correctly
read/processed by Gromologist - cf. the case of ILDN dihedral values:
```
#define torsion_ILE_N_CA_CB_CG2_mult1 0.0 0.8158800 1  ; Proteins 78, 1950 (2010)
```

To convert the keywords/variable names (like `torsion_ILE_N_CA_CB_CG2_mult1`) into their corresponding
values (like `0.0 0.8158800 1`) in the topology at hand, use:

```
>>> t.explicit_defines()
```

### Editing topologies

Let's start with a generic topology file:

```
>>> t = Top('md/topol.top')
```

##### Adding bonds within or between molecules

One useful application of Gromologist is adding bonds (and, automatically, other bonded terms)
either within a molecule or between them:

```
>>> protein = t.get_molecule("Protein")
>>> ligand = t.get_molecule("Other")
>>> t.list_molecules()
Protein                      1
Other                        1
>>> protein.merge_two(ligand, anchor_own=5, anchor_other=1)
>>> t.list_molecules()
Protein                      1
```

The above script merges Protein and Other into a single Protein molecule, adding a bond
between atom 5 of Protein and atom 1 of Other (here, indices are 1-based, corresponding
to numbering in .itp files).

To add a bond within a single e.g. Protein molecule, one can use `protein.merge_two(protein, anchor_own=2, anchor_other=3)`
or, more simply, `protein.add_bond(5,3)`.

##### Adding and removing atoms while maintaining ordered numbering

When an atom is removed, other atom numbers are modified accordingly, something that has to be
considered when removing multiple atoms. For instance, one can remove the first three atoms
in the following manner:

```
>>> protein.del_atom(1)
>>> protein.del_atom(1)
>>> protein.del_atom(1)
```

Note that all bonds, pairs, angles and dihedrals involving this atom are automatically removed as well.

To add an atom, one should specify its desired placement within the molecule, and at least 
a minimal set of parameters:

```
>>> protein.add_atom(atom_number=20, atom_name="CA", atom_type="CT")
By default, atom will be assigned to residue MET1. Proceed? [y/n]
y
>>> protein.add_bond(20,8)
```

If residue data is not specified, Gromologist will attempt to guess the residue based on
neighboring atoms.

##### Adding alchemical B-states

To generate alchemical states for a subset of atoms, one can use `gen_state_b`:

```
>>> protein.gen_state_b(atomtype='CT',new_type="CE")
```

The arguments for `gen_state_b` are divided into two subgroups:

 + `atomname`, `resname`, `resid` and `atomtype` behave as selectors, allowing to specify
 one or many atoms that should have its B-type specified;
 + `new_type`, `new_charge` and `new_mass` act as setters, allowing to specify the values
 in the B-state.
 
If the selectors point to multiple atoms (e.g. `atomtype=CT` selects all atoms with type CT),
all will be modified as specified. In turn, if a given setter is not specified, the respective 
value will be identical to that for state A.

##### Removing or swapping alchemical states

To make an alchemical topology non-alchemical again, one has two options:

+ To preserve state A, use:

```
>>> protein.drop_state_b()
```

+ To preserve state B as the only non-alchemical state, use:

```
>>> protein.drop_state_a()
```

If you want to invert the direction of the alchemical change by swapping states A and B, use:

```
>>> protein.swap_states()
```

##### Duplicating types

Often it's useful to duplicate an atomtype exactly, i.e., assign it a different name while
retaining all bonded and nonbonded parameters of the original. This can be done easily with:

```
>>> params = t.parameters
>>> params.clone_type("CT", prefix="Y")
```

This will create a type "YCT" that shares all properties with "CT" but can be modified independently.

##### Adding NBFIX terms

To generate an NBFIX (custom combination rule) entry, use the following snippet:

```
>>> t.parameters.add_nbfix(type1=CT, type2=HA, mod_sigma=0.01, mod_epsilon=-0.1)
```

This will introduce a term modifying the CT-HA Lennard-Jones interaction, increasing the default 
(Lorenz-Berthelot) sigma by 0.01 nm, and decreasing the default epsilon by 0.1 kJ/mol.

##### Explicitly listing parameters in topology 

To explicitly include all parameters in sections `[ bonds ]`, `[ angles ]` and `[ dihedrals ]`,
one can use:

```
>>> t.add_ff_params()
>>> t.save_top('with_params.top')
```

##### Finding missing parameters in topology

To find FF parameters that are missing (e.g. to include them by analogy, or optimize), run:

```
>>> t.find_missing_ff_params()
```

Note that both `add_ff_params()` and `find_missing_ff_params()` have an optional `section` parameter
that can specify you only want to look at `bonds`, `angles`, `dihedrals` or `impropers`.

### Dihedral optimization

With a completed Gaussian dihedral scan results at hand (.log file), we can use Gromologist
to run dihedral fitting. To select dihedral terms for refinement, add the `DIHOPT` keyword
anywhere in the dihedral term's comment (as many as you like, within reason) in the `.top` file
(or in `ffbonded.itp` in the FF directory you're using):

```
    CT    CT     N     C    9    0.000000   2.217520   1  ; phi,psi,parm94 DIHOPT
```

To run the optimization, simply use:

```
>>> from gromologist import DihOpt
>>> d = DihOpt(top='topol.top', qm_ref='gaussian_scan.log')
>>> d.optimize()
```

Upon termination, you will see a brief summary, and the resulting `opt1_topol.top` will 
contain optimized parameters. You can run `d.plot_fit()` to visualize the results,
and `d.restart()` to run refinement again starting from the optimized values. To control how
exhaustive the optimization is, both `.optimize()` and `.restart()` methods accept a `maxiter=N`
parameter determining the maximum number of iterations.

To perform multiple optimizations in parallel, add `processes=N` as a parameter to `DihOpt()`;
in this case, `N` runs will be initialized with different random seeds, and the best result
will be kept.

With Molywood installed, it is possible to use `d.make_movie()` to produce a movie illustrating
the structural aspects of the optimization (actively optimized dihedrals are highlighted in green)
along with a plot of the energy values.

### Editing structures

Let's start by reading a PDB file:

```
>>> p = Pdb('md/other.pdb')
```

##### Adding atoms along a vector specified by other atoms

To add e.g. a hydrogen atom to an existing atom CB, with a bond length of 1 A in the direction
specified by a vector from atom C to atom CA, one can use:

```
>>> p.insert_atom(len(p.atoms), base_atom=p.atoms[11], atomname='HX')
>>> p.reposition_atom_from_hook('name HX', 'name CB', 1.0, p1_sel='name C', p2_sel='name CA')
```

All the selections should be unique (corresponding to a single atom), and combinations can be
used like in VMD, e.g. `name CB and resid 10 and chain A`. This way you can e.g. automate
the addition of dummy atoms, DNA/RNA conversions etc.

##### Interpolating between two pre-aligned structures

To generate intermediate structures emulating a continuous conformational transition,
try the following snippet:

```
>>> p1 = Pdb('conf1_aligned.pdb')
>>> p2 = Pdb('conf2_aligned.pdb')
>>> p1.interpolate_struct(p2, num_inter=50, write=True)
```

This will create a total of 52 structures (1 starting + 50 intermediate + 1 final) named 
`interpolated_structure_{0..51}.pdb` that sample the transition through linear interpolation.

##### Filling beta-values with custom data (for visualization)

To use the PDB's beta column for color-mapping of observables e.g. in VMD, use the following:

```
>>> ca_atom_indices = p.select_atoms('name CA')
>>> p.fill_beta(per_residue_data, serials=[x+1 for x in ca_atom_indices])
>>> p.save_pdb('with_betas.pdb')
```

By adding the `smooth=...` parameter to `Pdb.fill_beta`, data can be spatially smoothed
using a Gaussian kernel with a specified standard deviation (in A).

##### Creating new PDB as a subset of existing one

To choose and save e.g. only the DNA atoms from a protein-DNA complex, use:

```
>>> dna =  Pdb.from_selection(p, 'dna')
>>> dna.save_pdb('dna.pdb')
```
