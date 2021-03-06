import RNAfoldz.utils
import textwrap as tw
import compbio.config as cfg
import os



def run_mcmc(run_id):
    run_id = '{0}_mcmc'.format(run_id)

def make_mcmc(run_id, struct, seqs, seqnames):
    run_id = '{0}_mcmc'.format(run_id)
    _clear_folder(run_id)
    _write_rna(run_id, struct, seqs, seqnames)
    outgroup_name = seqnames[0]
    _write_mcmc_ctl(run_id, outgroup_name)
    
def make_ml(run_id, struct, seqs, seqnames):
    run_id = '{0}_ml'.format(run_id)
    _clear_folder(run_id)
    _write_rna(run_id, struct, seqs, seqnames)
    outgroup_name = seqnames[0]
    _write_ml_ctl(run_id,outgroup_name)

def _clear_folder(run_id):
    datadir = cfg.dataPath('phase/{0}/'.format(run_id))
    for f in os.listdir(datadir):
        os.remove(os.path.join(datadir, f))
    

def _write_ml_ctl(run_id, outgroup_name):
    datafile = cfg.dataPath('phase/{0}/datafile.rna'.format(run_id))
    outfile = cfg.dataPath('phase/{0}/outfile.phylip'.format(run_id))
    cfgfile = cfg.dataPath('phase/{0}/control.ml'.format(run_id))

    data =  '''
#Phylogenetic tree reconstruction in the ML framework with mlphase
#The dataset in this example is small and mlphase can be used.

{DATAFILE}
Data file = %(datafile)s
Interleaved data file = no
#Use the "automatic method" to analyse this dataset:
#unpaired nucleotides ('.' in the secondary structure) are
#handled by the MODEL1 of the MIXED model (see below).
#pairs (corresponding parenthesis in the secondary structure)
#are handled by the MODEL2 of the MIXED model (see balow)
Heterogeneous data models = auto
{\DATAFILE}

''' % {'datafile':datafile}
    model = '''#Set up a MIXED model with REV for loops and 7D for stems
{MODEL}
Model = MIXED
Number of models = 2
  {MODEL1}
  Model = REV
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
  {\MODEL1}
  {MODEL2}
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
{\MODEL2}
{\MODEL}'''
    tree = '''
#A TREE block
{TREE}
#You must specify an outgroup although it is used for representation
#purpose only and it does not affect the results.
#This outgroup must be the name of a species in your datafile or the name
#of a clade in your clusters file (see below).
Outgroup = %(outgroup)s
#See manual for the available heuristic/exhaustive search method.
Search algorithm = Stepwise addition
#Optional: we specify a file that contains monophyletic clades. Tree topologies
#that do not match these constraints are not evaluated.
#Clusters file = sequence-data/hiv6.cls
{\TREE}
''' % {'outgroup': outgroup_name}
 
    run_cfg = '''
Random seed=9

Output file   = %(outfile)s

'''  % {'outfile':outfile}
    
    all_text = '\n'.join([data, model, tree, run_cfg])
    fopen = open(cfgfile, 'w')
    fopen.write(all_text)

    

def _write_mcmc_ctl(run_id, outgroup_name):
    datafile = cfg.dataPath('phase/{0}/datafile.rna'.format(run_id))
    outfile = cfg.dataPath('phase/{0}/outfile.phylip'.format(run_id))
    cfgfile = cfg.dataPath('phase/{0}/control.mcmc'.format(run_id))
    '''Write a control file for the mcmc tree builder in phase

'''
    data = '''
#A standard DATAFILE block for RNA sequences having a secondary structure.
#see also the sequence file sequence-data/mammals69.rna 

{DATAFILE}
Data file = %(datafile)s
Interleaved data file = no
#Use the "automatic method" to analyse this dataset:
#unpaired nucleotides ('.' in the secondary structure) are
#handled by the MODEL1 of the MIXED model (see below).
#pairs (corresponding parenthesis in the secondary structure)
#are handled by the MODEL2 of the MIXED model (see balow)
Heterogeneous data models = auto
{\DATAFILE}

''' % {'datafile':datafile}
    model = '''

#Set up a MIXED model with REV for loops and 7D for stems
{MODEL}
Model = MIXED
Number of models = 2
  {MODEL1}
  Model = REV
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
  {\MODEL1}
  {MODEL2}
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
{\MODEL2}
{\MODEL}

'''

    tree = '''

#Use a standard unrooted tree. The outgroup is compulsory but do not affect the results.
{TREE}
Tree = Unrooted MCMC tree
Outgroup = %(outgroup)s
{\TREE}

''' % {'outgroup': outgroup_name}
    perturbation = '''
    
#Tuning parameters for the MCMC runs. 
{PERTURBATION}

#relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 8
Model, proposal priority = 1

{PERTURBATION_TREE}
#We use 10/40 for topology change vs branch length changes.
#It is not exactly equivalent to 1/4 because this is also given relative to the
#proposal priority for hyperparameters that are introduced with the
#the prior on branch lengths (Hyperpriors, proposal priority)
Topology changes, proposal priority = 10
Branch lengths, proposal priority = 40
Hyperpriors, proposal priority = 1

#We use a vague prior exp(lambda) on branch lengths rather than the default exp(10)
Branch lengths, prior = exponential(uniform(0,100))
#A lambda hyperparameter has been introduced. It needs a "proposal priority"
#but this is not used because it is the only hyperparameter
Branch lengths exponential hyperparameter, proposal priority = 1
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
#relative probabilities for the proposals on the two models and the average substitution rate of MODEL2 
Model 1, proposal priority = 10
Model 2, proposal priority = 10
Average rates, proposal priority = 1
{PERTURBATION_MODEL1}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL1}
{PERTURBATION_MODEL2}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL2}
{\PERTURBATION_MODEL}

{\PERTURBATION}
'''
    run_cfg = '''

Random seed = 11


Output file   = %(outfile)s
Output format = phylip
'''  % {'outfile':outfile}


    #OLD VALUES:
    #Burnin iterations = 750000
    #Sampling iterations = 1500000
    #Sampling period = 150

    all_text = '\n'.join([data, model, tree, perturbation, run_cfg])
    fopen = open(cfgfile, 'w')
    fopen.write(all_text)

    
    return
    
def _write_mcmc_ctl(run_id, outgroup_name):
    datafile = cfg.dataPath('phase/{0}/datafile.rna'.format(run_id))
    outfile = cfg.dataPath('phase/{0}/outfile.phylip'.format(run_id))
    cfgfile = cfg.dataPath('phase/{0}/control.mcmc'.format(run_id))
    '''Write a control file for the mcmc tree builder in phase

'''
    data = '''
#A standard DATAFILE block for RNA sequences having a secondary structure.
#see also the sequence file sequence-data/mammals69.rna 

{DATAFILE}
Data file = %(datafile)s
Interleaved data file = no
#Use the "automatic method" to analyse this dataset:
#unpaired nucleotides ('.' in the secondary structure) are
#handled by the MODEL1 of the MIXED model (see below).
#pairs (corresponding parenthesis in the secondary structure)
#are handled by the MODEL2 of the MIXED model (see balow)
Heterogeneous data models = auto
{\DATAFILE}

''' % {'datafile':datafile}
    model = '''

#Set up a MIXED model with REV for loops and 7D for stems
{MODEL}
Model = MIXED
Number of models = 2
  {MODEL1}
  Model = REV
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
  {\MODEL1}
  {MODEL2}
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
{\MODEL2}
{\MODEL}

'''

    tree = '''

#Use a standard unrooted tree. The outgroup is compulsory but do not affect the results.
{TREE}
Tree = Unrooted MCMC tree
Outgroup = %(outgroup)s
{\TREE}

''' % {'outgroup': outgroup_name}
    perturbation = '''
    
#Tuning parameters for the MCMC runs. 
{PERTURBATION}

#relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 8
Model, proposal priority = 1

{PERTURBATION_TREE}
#We use 10/40 for topology change vs branch length changes.
#It is not exactly equivalent to 1/4 because this is also given relative to the
#proposal priority for hyperparameters that are introduced with the
#the prior on branch lengths (Hyperpriors, proposal priority)
Topology changes, proposal priority = 10
Branch lengths, proposal priority = 40
Hyperpriors, proposal priority = 1

#We use a vague prior exp(lambda) on branch lengths rather than the default exp(10)
Branch lengths, prior = exponential(uniform(0,100))
#A lambda hyperparameter has been introduced. It needs a "proposal priority"
#but this is not used because it is the only hyperparameter
Branch lengths exponential hyperparameter, proposal priority = 1
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
#relative probabilities for the proposals on the two models and the average substitution rate of MODEL2 
Model 1, proposal priority = 10
Model 2, proposal priority = 10
Average rates, proposal priority = 1
{PERTURBATION_MODEL1}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL1}
{PERTURBATION_MODEL2}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL2}
{\PERTURBATION_MODEL}

{\PERTURBATION}
'''
    run_cfg = '''

Random seed = 11

Burnin iterations = 750
Sampling iterations = 150
Sampling period = 150


Output file   = %(outfile)s
Output format = phylip
'''  % {'outfile':outfile}


    #OLD VALUES:
    #Burnin iterations = 750000
    #Sampling iterations = 1500000
    #Sampling period = 150

    all_text = '\n'.join([data, model, tree, perturbation, run_cfg])
    fopen = open(cfgfile, 'w')
    fopen.write(all_text)

    
    return

def _write_rna(run_id, struct, seqs, seqnames):
    '''
Write a datafile for the mcmc tree builder in phase.

Seqs should be specified simply as an (ascii) list of 
strings having values AUGC.

The file itself should have a first line giving:

nseqs, lenseqs, seqtype 
  eg:   '16 3571 STRUCT'

then the structure should be spec'd with '(.)'

and then the seqs in format: 

'NAME'   AUGC...
GCUGUGUGUGCUU... 

'NAME2'  AUGU...
AUAUAUUAUAAUA...

...

INPUTS:
 struct   (specified as pairs)
 seqs     (specified as strlist)
 seqnames (specified as s.sttrlist)

'''
    rutils = RNAfoldz.utils

    l = len(seqs[0])
    n = len(seqs)
    dtype = 'STRUCT'
    
    lines = '{0} {1} {2}\n'.format(n, l, dtype)
    lines += '\n'

    stk = rutils.pairs_stk(struct, l)    
    lines += '\n'.join(tw.wrap(stk)) + '\n\n'
    
    for seq, name in zip(seqs, seqnames):
        lines += name + '\n'
        lines += '\n'.join(tw.wrap(seq)) + '\n'
        lines += '\n'
    datafile = cfg.dataPath('phase/{0}/datafile.rna'.format(run_id))
    fopen = open(datafile, 'w')
    fopen.write(lines)

    return
