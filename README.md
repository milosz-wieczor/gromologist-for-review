# Gromologist

Gromologist is a package designed to facilitate handling, editing and manipulating GROMACS topology files 
(.top and .itp), as well as compatible structures (.pdb and .gro).

- [Gromologist](#gromologist)
  * [Installation](#installation)
  * [How to cite](#how-to-cite)
  * [Usage](#usage)
    + [Basic features](#basic-features)
        * [Reading and writing files](#reading-and-writing-files)
        * [File inspection, checks and printing](#file-inspection--checks-and-printing)
        * [Producing lightweight files](#producing-lightweight-files)
        * [Creating subsets of existing systems](#creating-subsets-of-existing-systems)
        * [Dealing with unspecified 'define' keywords in topologies](#dealing-with-unspecified--define--keywords-in-topologies)
    + [Editing topologies](#editing-topologies)
        * [Adding mutations to proteins](#adding-mutations-to-proteins)
        * [Adding bonds within or between molecules](#adding-bonds-within-or-between-molecules)
        * [Adding disulfide bonds and coordination bonds with transition metals](#adding-disulfide-bonds-and-coordination-bonds-with-transition-metals)
        * [Adding and removing atoms while maintaining ordered numbering](#adding-and-removing-atoms-while-maintaining-ordered-numbering)
        * [Adding alchemical B-states (including altered protonation)](#adding-alchemical-b-states)
        * [Adding position restraints](#adding-position-restraints)
        * [Removing or swapping alchemical states](#removing-or-swapping-alchemical-states)
        * [Duplicating and reassigning types](#duplicating-and-reassigning-types)
        * [Adding parameters or molecules from other topology files](#adding-external-params)
        * [Editing the contents of the system](#editing-contents-system)
        * [Modifying Lennard-Jones parameters](#modifying-lj-terms)
        * [Adding NBFIX terms](#adding-nbfix-terms)
        * [Adding H-mass repartitioning](#adding-hmass-repartitioning)
        * [Explicitly listing parameters in topology & finding missing parameters](#explicitly-listing-parameters-in-topology---finding-missing-parameters)
        * [Making .rtp entries](#making-rtp-entries)
        * [Preparing REST2 topologies](#preparing-rest2-topologies)
    + [Dihedral optimization](#dihedral-optimization)
    + [Editing structures](#editing-structures)
        * [Adding atoms along a vector specified by other atoms, and deleting them](#adding-atoms-along-a-vector-specified-by-other-atoms--and-deleting-them)
        * [Interpolating between two pre-aligned structures](#interpolating-between-two-pre-aligned-structures)
        * [Filling beta-values with custom data (for visualization)](#filling-beta-values-with-custom-data--for-visualization-)
        * [Renumbering atoms or residues in a structure](#renumbering-atoms-or-residues-in-a-structure)
        * [Adding chain, CONECT or element information](#adding-chain--conect-or-element-information)
        * [Converting a 3-point water model to a 4-point one](#converting-a-3-point-water-model-to-a-4-point-one)
    + [Selection language syntax](#selection-language-syntax)
    + [Access to Gromacs utilities](#access-to-gromacs-utilities)
        * [Energy decomposition for a structure or trajectory](#energy-decomposition-for-a-structure-or-trajectory)
        * [Automated system preparation](#automated-system-preparation)
        * [Miscellanous](#miscellanous)
    * [Sensitivity analysis](#sensitivity-analysis)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Installation
<a name="installation"/>

The latest "official" release can be obtained directly through `pip` by typing `pip install gromologist`.

To get the latest development version, first locally clone the git repo (`git clone https://gitlab.com/KomBioMol/gromologist.git`),
then install the package into Python by typing `pip install .` in the main Gromologist directory.
If you're using Anaconda, the same will work with `/path/to/anaconda/bin/pip`.

Note that Gromologist has no explicit dependencies, although a few specific functionalities will require 
basic scientific libraries such as `numpy` or `sklearn`.

## How to cite
<a name="how-to-cite"/>

If you found Gromologist useful, cite our preprint! [Gromologist: a Gromacs-Oriented Utility Library for Structure and Topology Manipulation
](https://chemrxiv.org/engage/chemrxiv/article-details/628bf40843d1f07e2e36c6e1)

## Usage
<a name="usage"/>

##### Reading and writing files
<a name="reading-and-writing-files"/>

`Top` and `Pdb` are the core classes of the library, and are supposed to provide representation
for topology and structure objects, respectively. To initialize them, a path to the file
should be passed to the constructor:

```
>>> from gromologist import Top, Pdb
>>> t = Top('md/topol.top')
>>> p = Pdb('md/conf.pdb')
>>> p
PDB file md/conf.pdb with 100 atoms
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
>>> t.pdb.save_gro('md/new_conf.gro')
```

##### File inspection, checks and printing
<a name="file-inspection--checks-and-printing"/>

If `Pdb` is bound to `Top`, a number of diagnostic and fixing options are available,
including name consistency checks:

```
>>> t.check_pdb()
Check passed, all names match
### or, alternatively:
>>> t.pdb.check_top()
Check passed, all names match
```

By default, up to 20 mismatches are shown, but more can be enabled with `maxwarn=...`.
If there are mismatches, in order to keep one naming convention and discard the other
use `fix_pdb=True` or `fix_top=True`.

If atom and residue names match but atom order does not, one can enforce correct ordering in Pdb
using topology ordering as a template:

```
>>> t.pdb.match_order_by_top_names()
```

Missing atoms in (standard) protein residues can also be identified rapidly:

```
>>> t.pdb.find_missing()
```

To check if all protein residues assume correct chirality, run:

```
>>> t.pdb.check_chiral_aa()
```

(By analogy to `Pdb.check_chiral_aa`, one can devise their own tests with `Pdb.check_chiral`.)

With `Pdb.get_atom_indices()`, selections can be made in a VMD-like manner (see full syntax
description at the bottom) to retrieve 0-based atom indices:

```
>>> t.pdb.get_atom_indices('name CA and (resname ASP or chain B)')
[5, 60, 72, 88]

```

while `Pdb.get_atoms()` returns a list of `Atom` instances instead of indices. Correspondingly,
`Pdb.get_atom_index()` and `Pdb.get_atom()` return a single object (int or Atom), raising an error
if the selection is compatible with more than one atom.

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
>>> protein.residues
# returns a list of residue_name-residue_id strings
```

By analogy, the `.list_bonds()` method can be used to `list_angles`, `list_dihedrals`
and `list_impropers`.

For structure files, several functions facilitate e.g. the identification of missing
atoms (useful when parsing files freshly downloaded from PDB):

```
>>> p = Pdb('1BCD.pdb')
>>> p.find_missing()
```

To print the sequence of the macromolecules in the file, use:

```
>>> p.print_nucleic_sequence()
>>> p.print_protein_sequence()
```

With `Pdb.print_protein_sequence()`, you can add `gaps=True` to fill in missing residues
with dashes (-), e.g. for alignment in Modeller.

To list all molecules in the Pdb instance (assigned from chains), use:

```
>>> p.print_mols()
```

and to simply print atoms in the structure (optionally: atoms corresponding to a selection), use:

```
>>> p.list_atoms()
>>> p.list_atoms(selection='within 5 of resid 8')
```

##### Producing lightweight files
<a name="producing-lightweight-files"/>

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

##### Creating subsets of existing systems
<a name="creating-subsets-of-existing-systems"/>

Gromologist features dedicated functions to create/save both topologies and structures 
of subsets of existing ones, where the subset can be defined with a simple text-based selection.

For example, to produce a structure containing DNA atoms from a protein-DNA complex, use:

```
>>> dna_structure = p.from_selection('dna')
# alternatively, if one simply wants to save the subset as PDB, there's another convenience function:
>>> p.save_from_selection('dna', 'dna_subsystem.pdb', renum=True)
# the above will save atoms corresponding to the 'dna' selection 
# in a file called 'dna_subsystem.pdb', renumbering the atoms from 1
```

Analogous functions exist for topologies. Importantly, the selections can easily cut through
molecules, and all the bonded terms will be processed correctly, yielding a simulation-ready topology:

```
>>> trinucleotide = t.from_selection('dna and resid 4 5 6')
# same convenience function can be used here:
>>> t.save_from_selection('dna and resid 4 5 6')
```

Molecules and parameters not relevant to the new subsystem will be removed.

##### Dealing with unspecified 'define' keywords in topologies
<a name="dealing-with-unspecified--define--keywords-in-topologies"/>

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
<a name="editing-topologies"/>

Let's start with a generic topology file:

```
>>> t = Top('md/topol.top')
```

##### Adding mutations to proteins
<a name="adding-mutations-to-proteins"/>

For certain complex systems, having to pass through pdb2gmx or CHARMM-GUI 
for every mutant is a major drawback. To avoid this, Gromologist allows to insert
amino acid mutations into existing topologies, preserving all their existing features.
This is as easy as the following snippet shows:

```
>>> protein = t.get_molecule("Protein")
>>> protein.mutate_protein_residue(2, "Y")
>>> t.save_top("x2y_mutant.top")
```

If Gromacs files are not found by Gromologist, the `Top.Section.mutate_protein_residue` 
function can accept an optional `rtp=/path/to/rtp/file` argument to specify residue parameters.

If the `Top` object has a `Pdb` object bound to it, by default the mutation will be introduced
to both the topology and structure: 

```
>>> t = Top('md/topol.top', pdb='md/conf.pdb')
>>> protein = t.get_molecule("Protein")
>>> protein.mutate_protein_residue(2, "Y")
>>> t.save_top("x2y_mutant.top")
>>> t.pdb.save_pdb("x2y_mutant.pdb")
```

Structures alone can be mutated using the associated
`Pdb.mutate_protein_residue(resid, target, chain='')` function.

Note that mutations in the structure are not guaranteed to be free from clashes, so always
make sure your final structure is acceptable and run an energy minimization prior to running
dynamics (note the double precision version of Gromacs, `gmx_d`, works better for atoms 
that almost overlap in space).

##### Adding bonds within or between molecules
<a name="adding-bonds-within-or-between-molecules"/>

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

##### Adding disulfide bonds and coordination bonds with transition metals
<a name="adding-disulfide-bonds-and-coordination-bonds-with-transition-metals"/>

To fix issues with Gromacs' specbonds directives, Gromologist can automatically add disulfide
bonds between two cysteine residues. An .rtp file is required (either selected interactively
or passed as an argument) to make sure that charges and types are assigned properly:

```
>>> protein_a = t.get_molecule("Protein_chain_A")
>>> protein_b = t.get_molecule("Protein_chain_B")
# the following line adds an intramolecular disulfide between residues 15 and 30 of chain A:
>>> protein_a.add_disulfide(15, 30)
# the following line adds an intermolecular disulfide between residues 15 of chain A and 30 of chain B:
>>> protein_a.add_disulfide(15, 30, other=protein_b)
```

The same can be done with coordination bonds involving Cys/His and Zn/Fe:

```
>>> protein = t.get_molecule("Protein_chain_A")
>>> ion = t.get_molecule("Zn")
>>> protein.add_coordinated_ion(15, 1, other=ion)
```

When the amino acid is Cys, an .rtp file is required (either path specified as `rtp=...` or selected 
interactively), as Cys has to deprotonate and types/charges change. In case of His, the proper protonation
state is recognized and the bond is added between the metal ion and the nitrogen with a lone pair, but 
no changes are made to types/charges.

Note that even though excess atoms (e.g. HG of cysteine) are automatically removed from the accompanying structure file,
the resulting bonds can be very long and generate high energies/forces when simulated or minimized. Make sure you know
what you are doing!

##### Adding and removing atoms while maintaining ordered numbering
<a name="adding-and-removing-atoms-while-maintaining-ordered-numbering"/>

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

##### Adding alchemical B-states (including altered protonation)
<a name="adding-alchemical-b-states"/>

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

To make an alchemical residue that switches between a protonated and deprotonated version,
use the following for any Asp or Glu residue (be it protonated or deprotonated in the first place):

```
>>> protein.alchemical_proton(31)
```

You will be asked for an .rtp file that contains the respective protonated/deprotonated
residue pair, and Gromologist will try to make residue number 31 alchemical.

##### Adding position restraints
<a name="adding-position-restraints"/>

If you want to specify new position restrains for your molecule (or you deleted them previously),
use:

```
>>> t.add_posres(keyword = 'POSRES', value = 1000)
```

will add a position restraint with a force constant of 1000 to all heavy atoms of molecules 
larger than 5 atoms (to exclude water models), conditional on the `-DPOSRES` directive in .mdp.

The same add_posres() method of a molecule object can be used, and selections can be specified
to restrict the number of heavy atoms to which the restraint is applied. When `keyword` is omitted,
the section becomes unconditional (will always be turned on during a simulation). Moreover, a list
of force constants can be supplied to specify restraints differently in X, Y and Z dimensions:

```
>>> protein.add_posres(value = [200, 200, 0], selection='resid 47 to 83')
```

will apply an unconditional restraint in the XY plane to heavy atoms of residues 47-83 in `protein`. 

##### Removing or swapping alchemical states
<a name="removing-or-swapping-alchemical-states"/>

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

##### Duplicating and reassigning types
<a name="duplicating-and-reassigning-types"/>

Often it's useful to duplicate an atomtype exactly, i.e., assign it a different name while
retaining all bonded and nonbonded parameters of the original. This can be done easily with:

```
>>> params = t.parameters
>>> params.clone_type("CT", prefix="Y")
>>> params.clone_type("CT", new_type="CY")
```

This will create a type "YCT" (2nd line) and a type "CY" (3rd line) that both share 
all properties with "CT" but can be modified independently.

To then set e.g. all CA atoms in the 1st molecule to the new type, run the following:

```
>>> t.molecules[0].set_type("CY", atomname="CA")  # sets all CAs's type as CY
>>> t.molecules[0].set_type("CY", atomname="CA", resname="ALA")  # sets CAs in ALA as CY
>>> t.molecules[0].set_type("CY", atomname="CA", resname=["ALA", "LYS"])  # sets CAs in ALA and LYS as CY
```

##### Adding parameters or molecules from other topology files
<a name="adding-external-params"/>

To incorporate and merge parameters from another topology (e.g. a ligand parametrization
generated elsewhere), use the following:

```
>>> t.add_parameters_from_file('ligand.top')
```

To "import" molecule definitions in the same way, go for:

```
>>> t.add_molecule_from_file('ligand.top')
>>> t.add_molecule_from_file('ligand.top', ['LIG']) # to import only selected molecules
```

##### Editing the contents of the system
<a name="editing-contents-system"/>

To add a specified number of molecules to the system (note, this is different than 
just adding molecule definitions!), use:
 
```
>>> t.add_molecules_to_system(molname="LIG", nmol=3)
```
In turn, to completely remove the molecule (both definition & in system contents),
you can use:

```
>>> t.remove_molecule(molname="LIG")
```

##### Modifying Lennard-Jones parameters
<a name="modifying-lj-terms"/>

To change the values of sigma or epsilon for a given type, use:

```
>>> t.parameters.edit_atomtype('CT', mod_sigma=0.01, mod_epsilon=-0.1)
```

##### Adding NBFIX terms
<a name="adding-nbfix-terms"/>

To generate an NBFIX (custom combination rule) entry, use the following snippet:

```
>>> t.parameters.add_nbfix(type1='CT', type2='HA', mod_sigma=0.01, mod_epsilon=-0.1)
```

This will introduce a term modifying the CT-HA Lennard-Jones interaction, increasing the default 
(Lorenz-Berthelot) sigma by 0.01 nm, and decreasing the default epsilon by 0.1 kJ/mol.

##### Adding H-mass repartitioning
<a name="adding-hmass-repartitioning"/>

Hydrogen mass repartitioning is a popular method of extending the timestep in atomistic molecular
simulations to 4 fs by shifting the mass from heavy atoms to hydrogens, thereby slowing down the 
most rapid oscillations in the system (ones involving hydrogen atoms). 

One can apply hydrogen mass repartitioning in `gmx pdb2gmx` with the `-heavyh` keyword, but occasionally
it is more convenient to apply it to an existing topology:

```
>>> t.hydrogen_mass_repartitioning()
```

By default, the masses of hydrogens are set to 4.032 Da; this can be controlled with the `hmass` keyword.
The method skips all molecules with up to 5 atoms (intended not to affect water molecules). 


##### Explicitly listing parameters in topology & finding missing parameters
<a name="explicitly-listing-parameters-in-topology---finding-missing-parameters"/>

To explicitly include all parameters in sections `[ bonds ]`, `[ angles ]` and `[ dihedrals ]`,
one can use:

```
>>> t.add_ff_params()
>>> t.save_top('with_params.top')
```

To find FF parameters that are missing (e.g. to include them by analogy, or optimize), run:

```
>>> t.find_missing_ff_params()
```

Note that both `add_ff_params()` and `find_missing_ff_params()` have an optional `section` parameter
that can specify you only want to look at `bonds`, `angles`, `dihedrals` or `impropers`.

Additionally, one can recalculate the cumulative charge ("qtot" comment in molecule definitions)
by using:

```
>>> t.recalculate_qtot()
```

This does not affect the energy of the system in any way, and only serves for informative purposes.

##### Making .rtp entries
<a name="making-rtp-entries"/>

To convert a residue (possibly defined in an .itp file coming from a parametrization tool) 
into an .rtp entry (that can be used by `gmx pdb2gmx`), you can load the topology and use:

```
itp = Top('new_residue.itp')
itp.molecules[0].tp_rtp('new_entry.rtp')
```

Then copy the entry to the .rtp database of your chosen force field. Do not forget to update your
`residuetypes.dat`! 

##### Preparing REST2 topologies
<a name="preparing-rest2-topologies"/>

To prepare a molecule for replica exchange/solute tempering simulations, one can use the 
top-level Top object:

```
>>> t.solute_tempering(temperatures=[300, 320, 350], molecules=[0,1])
```

If the "hot" region is localized within one molecule, one can use the lower-level interface 
instead, and specify a selection for the "hot" subsystem:

```
>>> t.molecules[0].solute_tempering(temperatures=[300, 320, 350], selection='resid 1 to 50')
```

### Dihedral optimization
<a name="dihedral-optimization"/>

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
<a name="editing-structures"/>

Let's start by reading a PDB file:

```
>>> p = Pdb('md/other.pdb')
```

##### Adding atoms along a vector specified by other atoms, and deleting them
<a name="adding-atoms-along-a-vector-specified-by-other-atoms--and-deleting-them"/>

To add e.g. a hydrogen atom "hooked" to an existing atom CB in residue 2, with a bond length of 1 A in the direction
specified by a vector from atom C to atom CA, one can use:

```
>>> p.insert_atom(serial=11, name='HX', hooksel='name CB and resid 2', 1.0, p1_sel='name C and resid 2', p2_sel='name CA and resid 2')
```

The new atom will have its residue name, residue number and chain copied from the "hook" atom.

All the selections should be unique (corresponding to a single atom), and combinations can be
used like in VMD when necessary, e.g. `name CB and resid 2 and chain A`. This way you can e.g. automate
the addition of dummy atoms, DNA/RNA conversions etc.

For advanced users, it is possible to directly specify the `vector=...` parameter instead of
`p1_sel` and `p2_sel`, while the vector can be conveniently calculated using `Pdb._vector()`
(in fact, this is how the mutation module is implemented).

Atoms can be easily deleted with `Pdb.delete_atom()` using the serial number, with the optional
`renumber=True` parameter to renumber the remaining atoms from 1. The `Pdb.remove_hydrogens()` function
automatically removes all hydrogen atoms from the structure.

##### Interpolating between two pre-aligned structures
<a name="interpolating-between-two-pre-aligned-structures"/>

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
<a name="filling-beta-values-with-custom-data--for-visualization-"/>

To use the PDB's beta column for color-mapping of observables e.g. in VMD, use the following:

```
>>> p.set_beta(per_residue_data, selection='name CA')
>>> p.save_pdb('with_betas.pdb')
```

By adding the `smooth=...` parameter to `Pdb.set_beta`, data can be spatially smoothed
using a Gaussian kernel with a specified standard deviation (in A).


##### Renumbering atoms, residues or chains in a structure
<a name="renumbering-atoms-or-residues-in-a-structure"/>

`Pdb.renumber_atoms()` and `Pdb.renumber_residues()` serve to easily reset the numbering
of atoms or residues, respectively. The renumbering can be modified with the `offset` 
and `selection` parameters:

```
>>> p.renumber_residues(offset=20)  # starts numbering residues from 20
>>> p.renumber_atoms(selection='chain B')  # only renumbers atoms in chain B
```

For when it is desirable to change the order of chains in a system, use `Pdb.permute_chains()`:

```
>>> p.permute_chains([2, 1, 0], rename=True)
```

This will put chain nr 3 first, then nr 2, and then nr 1 (C, B, A if they were originally A, B, C).

##### Adding chain, CONECT or element information
<a name="adding-chain--conect-or-element-information"/>

When chain information goes missing (common issue with conversion between .pdb and .gro),
this information can be easily recovered with `Pdb.add_chains()`. Note: if the `Pdb` instance
does not have a bound `Top` object, it will try to guess chain end/start based on the `cutoff`
parameter (default is 10 A); otherwise, it will first check Pdb/Top for consistency and use 
molecule information from `Top` to assign chains.

```
>>> p.add_chains()
```

In order to only assign chains to proteins, use `protein_only=True`. To start numbering chains
from a letter different than A, use `offset=X`, where X is a 0-based integer.

If the distance calculation fails due to PBC issues, use `nopbc=True` and **make sure** that
your molecule is whole.

To guess the information about elements (e.g. to add element-based coloring in VMD),
use:

```
>>> p.add_elements()
```

Finally, if CONECT entries are needed in the PDB (as required by some programs), one can run:

```
>>> p.add_conect()
```

with a default distance cut-off of 1.55 A to define a chemical bond. Note that this last
feature depends on `numpy` for speed!

##### Converting a 3-point water model to a 4-point one
<a name="converting-a-3-point-water-model-to-a-4-point-one"/>

If your system was prepared e.g. with TIP3P water and you want to change it to OPC (or other 4-point)
without rerunning `gmx pdb2gmx`, the following function will do the job:

```
>>> p.tip3_to_opc()
```

The default `offset` parameter is set to 0.147722363, typical for OPC; for TIP4P, use
`offset=0.128012065`, and in general check the `[ virtual_sites3 ]` section in the solvent's .itp file.

### Selection language syntax
<a name="selection-language-syntax"/>

The custom selection language was meant to be as similar as possible to that 
available in VMD: 
+ keywords such as "name", "resid", "resname", "serial", "chain" (for structures) or "type" (for topologies)
+ phrases "same residue as ..." and "within ... of" 
+ predefined selections include "protein", "dna", "rna" and "solvent"
+ logical operators as "or" & "and" work on sets of atoms; "not" inverses the selection
+ ranges can be specified as e.g. "resid 5 to 25" or by simple enumeration
+ parentheses can be used to customize order of operations

Examples: 
+ "(resid 1 to 100 and name CA) or (same residue as within 5 of resname LIG)"
+ "chain A B D and not solvent"

### Access to Gromacs utilities
<a name="access-to-gromacs-utilities"/>

##### Energy decomposition for a structure or trajectory
<a name="energy-decomposition-for-a-structure-or-trajectory"/>

To perform energy decomposition using the Gromacs rerun module, use the
`calc_gmx_energy()` utility function:

```
>>> import gromologist as gml
>>> gml.calc_gmx_energy(struct='md/conf.pdb', topfile='md/topol.top', traj='traj.xtc', terms='all')
```

The `terms` keyword can be "all" (returns all available energy/pressure/volume terms), 
any specific keyword allowed by `gmx energy` (e.g. "potential"), 
or a list of these (e.g. ["bonds", "angles", "potential"]).

To specifically calculate the interaction energy between two groups, similarly use:

```
>>> gml.calc_gmx_energy(struct='md/conf.pdb', topfile='md/topol.top', traj='traj.xtc', group_a = 'protein',
                    group_b='dna', sum_output=True, savetxt='energies.dat')
```

This will calculate the Coulombic and Lennard-Jones terms for the interaction between groups defined
by selections `protein` and `dna`, add a dataset corresponding to their sum, and save these 3 terms
to a file `energies.dat`.


##### Automated system preparation
<a name="automated-system-preparation"/>

To prepare a standard simulation box starting from a .pdb structure, it is enough to run:

```
>>> gml.prepare_system('starting.pdb')
```

You will be prompted to choose a force field based on the .ff directories found in the Gromacs 
installation and the local directory, and a standard sequence of Gromacs commands will be executed:
+ `gmx pdb2gmx` to generate the topology
+ `gmx editconf` to set boxsize (default is dodecahedron with a 1.5 nm buffer)
+ `gmx solvate` to add the water model compatible with the desired force field
+ `gmx grompp` and `gmx genion` to add selected ions at a specified concentration (by default 150 mM KCl)
+ `gmx grompp` and `gmx mdrun` to perform energy minimization
+ `gmx trjconv` to make the resulting molecule whole, 
eventually yielding `minimized_structure.pdb` and `merged_topology.top`.

To customize or automate the workflow, the following options are available:

+ ff: a string ending with .ff to specify the force field
+ water: a string naming the water model
+ box: a string specifying the box type, e.g. 'cubic
+ cation: a string specifying the cation type, e.g. 'Na' 
+ anion: a string specifying the anion type, e.g. 'Br'
+ ion_conc: a float specifying the ion concentration in molar units
+ maxsol: an int specifying the maximum number of water molecules
+ resize_box: a bool specifying whether or not to set a new size for the box
+ box_margin: a float defining the margin around the solute in nm
+ explicit_box_size: a list containing the a, b, c components of the cubic box vectors in nm 
+ topology: a string with a path to an existing topology of a solute (will skip `gmx pdb2gmx`)

Any additional key:value pairs will be passed to pdb2gmx (e.g. `heavyh='yes'` to turn on
hydrogen mass repartitioning).

##### Miscellanous
<a name="miscellanous"/>

`gml.frames_count('file.xtc')` quickly returns the number of frames in an .xtc file.

`gml.load_frcmod(gml.top, ff.frcmod)` loads all the additional force field entries
from Amber's `frcmod` files into a Gromacs `.top` topology.

`gml.compare_topologies_by_energy('struct.pdb', 'top1.top', 'top2.top')` checks if
two topologies yield identical potential energies when evaluated with a chosen structure
(useful when testing converters or other topology modifiers).


### Sensitivity analysis
<a name="sensitivity-analysis"/>
To perform a full sensitivity analysis in the NBFIX space, start by calculating
the energy derivatives for each frame for each possible NBFIX:

```
>>> import gromologist as gml
>>> td = gml.ThermoDiff()
>>> # this adds all possible NBFIXes to the list of calculated sensitivities:
>>> td.add_all_nbfix_mods(top='md/topol.top', structure='md/conf.pdb')
>>> # this specifies a trajectory on which the sensitivity will be calculated, as well as relevant datasets:
>>> hdata = np.loadtxt('helix_content.dat')[:,1]
>>> td.add_traj(top='md/topol.top', traj='md/traj.xtc', datasets={'helicity': hdata})
>>> td.run() # this part will take some time
>>> # let's find the difference between the binned derivatives for the lower and upper half of the dataset:
>>> hmin, hmax = np.min(hdata), np.max(hdata)
>>> hmid = 0.5 * (hmin + hmax)
>>> td.calc_discrete_derivatives(dataset='helicity', threshold=[hmin, hmid, hmid, hmax])
```