Optative flags
###########################

Optional flags for each package:

    - `All packages <all_packages/index.html>`_
    - `Protein-Protein Inhibitors <ppi/index.html>`_
    - `Pocket Exploration (Allosteric) <pocket_exploration/index.html>`_
    - `HT-Fragment Growing <frag/index.html>`_
    - `AquaPELE (water perturbation) <water/index.html>`_

.. toctree::
   all_packages/index.rst
   :hidden:

.. toctree::
   ppi/index.rst
   :hidden:

.. toctree::
   pocket_exploration/index.rst
   :hidden:

.. toctree::
   frag/index.rst
   :hidden:

.. toctree::
   water/index.rst
   :hidden:

Box parameters
=================

Parameters to set the exploration Box:

- **box_radius**: Radius of the box. Default=[induced_fit (10), local_exploration (30), global_exploration (50)]

- **box_center**: Center of the box. Default=[indeuced_fit&local_exploration (CM of the ligand), global (calculater center)]


..  code-block:: yaml

  box_radius: 30
  box_center: 
    - 20
    - 30
    - 50


Simulation params
====================

- **seed**: Seed of the job for reproducibility. Default=12345

- **log**: Retrieve PELE logfiles during simulation. Default=False

- **verbose**: Set to true to activate verbose mode in PELE. DEfault=False

- **anm_freq**: Every how many steps to perform anm. Default=4

- **anm_displacement**: Angstrom to displace carbon alphas in each ANM movement. Default=0.75

- **anm_modes_change**: Number of steps before we change to a new normal mode movement. Default=4

- **sidechain_freq**: Every how many steps to perform sidechain sampling. Default=2

- **min_freq**: Every how many steps to perform minimization. Default=1

- **water_freq**: Every how many steps to perform water perturbation. Default=1

- **temperature**: Temperature of the simulation. Default=1500

- **solvent**: Solvent of the simulation. (OBC or VDGBNP). Default=VDGBNP

- **sidechain_res**: Receptor sidechain resolution. Default=10

- **overlap_factor**: Vanderwals overlap factor (More in PELE docs). Default=0.65

- **steric_trials**: Number of steric trials (More in PELE docs). Default=250

..  code-block:: yaml

  seed: 312312
  log: true
  verbose: true
  anm_freq: 4
  anm_displacement: 0.5
  anm_modes_change: 3
  sidechain_freq: 2
  min_freq: 1
  water_freq: 1
  temperature: 1500
  solvent: "VDGBNP"
  sidechain_res: 30
  overlap_factor: 0.65
  steric_trials: 250



PELE params
===================

**These flags are exclusive of the PELE modes not fragPELE**

- **iterations**: Adaptive epochs to run. Set to 1 by default if using PELE

- **steps**: Pele steps in each iteration

- **debug**: Use this flag to only create the inputs of the simulation. No simulation is run. (Usefull to transport it to another machine)

- **spawning**: Spawning type ([independent, inverselyProportional or epsilon so far]). Default: inverselyProportional

- **density**: Density type ([null, exitContinuous...]. More in AdaptivePELE docs). Default: null

- **cluster_values**: Clusterization values. More in AdaptivePELE. Default: Depending on simulation type

- **cluster_conditions**: Clusterization condition. More in AdaptivePELE. Default: Depending on simulation type

- **equilibration**: Whether to run initial equilibration or not. Default: false

- **equilibration_steps**: Equilibration steps. Default: 2
  
- **adaptive_restart**: Use adaptive restart with the working folder option to restart the simulation. Default: false

- **report**: Change the name of the report file. Default: report

- **traj**: Change the name of the trajectory file. Default: trajectory.pdb

..  code-block:: yaml

    iterations: 30
    steps: 12
    debug: true
    spawning: "epsilon"
    density: "exitContinuous"
    cluster_values: [2,3,4]
    cluster_conditions: [0.8, 0.6, 0.2]
    equilibration: false
    equilibration_steps: 10
    adaptive_restart: true
    working_folder: "folder_to_restart"
    report: report
    traj: trajectory.xtc

FragPELE params
===================

**These flags are exclusive of the FragPele modes not PELE**

- **growing_steps**: Number of steps to grow the fragment with.

- **steps_in_gs**: Number of pele steps within each growing step

- **sampling_steps**: Number of pele steps in the final sampling simulation

- **protocol**: Type of protocol. options = [HT, ES]. For more info please refere here.


..  code-block:: yaml

    growing_steps: 6
    steps_in_gs: 6
    sampling_steps: 20
    protocol: HT
    cpus: 24

PPI params
===============

**These flags are exclusive of the ppi: true mode**

- n_components: Number of clusters after global exploration. In other words, number of inputs for the refinment exploration after the global simulation. Default: 10


..  code-block:: yaml

    n_components: 10


Constraints
==================

- **water_constr**: Water constraints. Default=5

- **constrain_smiles**: SMILES string to indicate what part of the molecule to constrain. Default=None

- **smiles_constr**: Numeric value of the SMILES constraints. Default=10

- **external_constraints**: You can specify 2 types of constraints. Positional constraints or atom-atom constraint. (Example below)

  - The positional constraints are given either by: 
        - springConstant-atomnumber. i.e. "10-17"
        - springConstant-chain:resnum:atomname. i.e. "5-A:1:H"

  - The atom-atom constraints are specified either by: 
        - springConstant-equilibriumDistance-atomnumber1-atomnumber2. i.e. "50-2.34-17-4159"
        - springConstant-equilibriumDistance-chain1:resnum1:atomname1-chain2:resnum2:atomname2. i.e. "50-2.34-A:1:H-L:1:C21"

- **remove_constraints**: Do not place constraints on the carbon-alpha of the protein. Default: False


..  code-block:: yaml

    water_constr: 5
    constrain_smiles: "C2CCC1CCCCC1C2"
    smiles_constr: 5
    external_constraints:
    - "10-17" #constrain of 10kcal/mol at atomnumber 17
    - "5-A:1:H" ##constrain of 10kcal/mol at atom with chain A residuenumber 1 and atomname H
    - "50-2.34-17-4159" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atomnumbers 17 & 4159
    - "50-2.34-A:1:H-L:1:C21" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    remove_constraints: true

Carbon-alpha constraints
+++++++++++++++++++++++++

Each package in the platform has its own predefined constraint parameters which are likely to be the best choice in each
type of study. However, the platform provides the users with several different levels of constraining the alpha carbons
of the protein backbone with varying spring constants and intervals:

- **level 0** - no constraints

- **level 1** - terminal CAs constrained with a spring constant of 5 kcal/mol, the rest of the CAs in the backbone with 0.5 kcal/mol at an interval of 10, i.e. every 10 residues (default)

- **level 2** - terminal CAs constrained at 5 kcal/mol, the rest of the CAs with 2.5 kcal/mol at the interval of 8 (default for the ``rescoring`` package)

- **level 3** - the whole backbone is constrained every 5 atoms with 5 kcal/mol (default for the ``gpcr_orth`` package)

We strongly suggest relying on the default settings for each package. However, in case of studying a system where the
defaults are not optimal (more flexibility or rigidity required), the users can change the level, for example:

..  code-block:: yaml

    constraint_level: 3

Alternatively, advanced users can manipulate the constraint parameters individually at their own risk, using the following flags:

- **terminal_constr** - sets the spring constant for the terminal C-alpha constraints, default = 5 kcal/mol

- **ca_constr** - sets the spring constant for the remaining C-alphas in the backbone, default = 0.5 kcal/mol

- **ca_interval** - interval at which the backbone C-alphas should be constrained, default = 10 (i.e. every 10 residues).

Take into account that specific modifiers of constraint parameters will prevail over the settings coming from the
constraints levels and those predefined in each package.

..  code-block:: yaml

    terminal_constr: 10.5
    ca_constr: 6.0
    ca_interval: 3

Metal constraints
+++++++++++++++++++++

Algorithm to automatically set metal constraints around the ligand.

- **no_metal_constraints**: Ignore all metals in the PDB file, no constraints will be set automatically. Default=False

- **permissive_metal_constr**: Expand the search for coordinated atoms by allowing 35% deviation from “ideal” angles. If the algorithm finds a valid geometry it will include the metal constraint into the simulation. Default=False

- **constrain_all_metals**: Constrain all atoms around the metal, regardless of the angles or coordination number. Default=False

- **external_constraints**: Set a manual constraint containing a metal atom to disable search for this particular metal. Default=[]


..  code-block:: yaml

    no_metal_constraints: true
    permissive_metal_constr: true
    constrain_all_metals: true
    external_constraints:
        - "50-2.34-A:1:H-L:1:MG" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    constrain_core: "CN(C)C(=O)c1ccc(F)cc1"  # SMILES or SMARTS pattern
    constrain_core_spring: 30  # optional, default 50.0


Core constraints
+++++++++++++++++++++

You can constrain the core of your ligand by specifying either SMILES or SMARTS pattern using ``constrain_core`` flag.
The default spring constant is 50 but you can choose your own.

..  code-block:: yaml

    constrain_core: "CN(C)C(=O)c1ccc(F)cc1"  # SMILES or SMARTS pattern
    constrain_core_spring: 30  # optional, default 50.0

WaterPerturbation
======================

- **n_waters**: Number of waters to randomly add into your simulation. Compulsory when running MonteCarlo with water perturbation. Default=0

- **box_water**: Center of the box for the waters. Default: Centroid of the center of masses of all water molecules.

- **water_radius**: Radius of the water box. Default=7

- **water_trials**: Numerical trials on water perturbation. Default=10000

- **water_constr**: COM constrain applied to th water molecule after perturbation. Default=0

- **water_temp**: Temperature of the water perturbation step. Default=5000

- **water_overlap**: Overlap factor of water. Default=0.78


..  code-block:: yaml

    n_waters: 3 #Compulsory
    box_water:
    - 20
    - 30
    - 20
    water_radius: 8
    water_trials: 500
    water_constr: 0.5
    water_temp: 2000
    water_overlap: 0.5

Interaction restrictions
=========================

Interaction restrictions allow for biased exploration, where the simulation results are limited to those that fit the specified conditions.

Users can define two types of conditions using the atom strings (format "chain:resnum:atomname", e.g. A:2:CA) to select the atoms:

- **distance**: Distance between two atoms, which can be limited to a user-defined maximum, minimum or both.

- **angle**: Angle between three atoms with a user-defined maximum, minimum or both.


..  code-block:: yaml

    interaction_restrictions:
    - distance:  # distance between the two atoms will not exceed 3 A
        max: 3
      atoms:
        - "A:318:OG1"   # chain A, residue number 318, atom OG1
        - "Z:201:O3"
    - angle:  # angle between those three atoms will remain betwenn 90 and 180 degrees
        min: 90
        max: 180
      atoms:
        - "A:318:OG1"
        - "A:318:HG1"
        - "Z:201:O3"


Metrics
=============

Metrics to track along the simulation

- **atom_dist**: Calculate distance between two atomnumbers. To calculate more than one append them in column as the example below. Default=None

    - The atomdist can be specified via chain:resnum:atomname i.e. A:2:CA

- **rmsd_pdb**: Calculate rmsd of the ligand to a native pdb structure


..  code-block:: yaml

    atom_dist:
        # Distance between the A:2:CA and B:3:CG also between A:5:N and B:3:CG. Append more if desired.
        - "A:2:CA"
        - "B:3:CG"
        - "A:5:N"
        - "B:3:CG"
    rmsd_pdb: "/home/dsoler/native.pdb"


Analysis
=============

Run a post simulation analysis to extract plots, top poses and clusters.

- **only_analysis**: Analyse PELE simulation without running it.

- **analysis_nclust**: Numbers of clusters out of the simulation, if using the standard clustering method. Default: 10

- **be_column**: Column of the binding energy in the reports starting by 1. Default: 5

- **te_column**: Column of the total energy in the reports starting by 1. Default: 4

- **limit_column**: Specify the column where your external metrics start. Default: 6

- **mae**: To extract the best energy and cluster poses as .mae files with the metrics as properties (schrodinger need it). Default: false

- **analysis**: Whether to run or not the analysis at the end of the simulation. Default: true

- **clustering_method**: If you want to override the default clustering method (Gaussian mixture model), you can set this flag to ``MeanShift`` or ``HDBSCAN``.

- **bandwidth**: Value for the Mean Shift bandwidth (when using the Mean Shift algorithm) or epsilon (when using the HDBSCAN clustering); default = 5.0

..  code-block:: yaml

    only_analysis: true
    be_column: 5
    te_column: 4
    limit_column: 6
    mae: true
    clustering_method: "meanshift"
    bandwidth: 7.0

The bandwidth parameter hugely influences the clustering results, therefore, it might be worth trying out different values depending on your system.
In case of the mean shift algorithm, the bandwidth refers to the maximum RMSD allowed within the cluster, whereas in HDBSCAN to distances between your data points.

Output
==========

Configure the output

- **working_folder**: Name of the main working folder where to store the processed input, control files and the simulation folder. Default="resname_Pele_X" where X is a number.

- **output**: Output folder of the simulation. Default=output

..  code-block:: yaml

    working_folder: "NOR_solvent_OBC"
    output: "output_sim"
