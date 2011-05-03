from modeller import *
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='query-3BIX.ali',
              knowns='3bix_A', sequence='query',
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
