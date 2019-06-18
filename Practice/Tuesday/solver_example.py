from java.util import ArrayList
from xal.extension.solver import Scorer, Trial, Variable, Stopper, Solver

from xal.extension.solver import SolveStopperFactory, ProblemFactory, Problem

from xal.extension.solver.algorithm import SimplexSearchAlgorithm 
from xal.extension.solver.algorithm import RandomShrinkSearch, RandomSearch
from xal.extension.solver.hint import Hint, InitialDelta


import sys
import math
import types
import time
import random

from jarray import *
from java.lang import *
from java.util import *
from java.io import *

from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

class MyScorer(Scorer):
    def __init__(self,variables,func):
        self.variables = variables
        self.func = func
		#--------------------------------
        self.diff2_min = Double.MAX_VALUE
        self.count = 0
    def score(self,trial,variables):
        self.count += 1
        var = self.variables.get(0)
        x =  trial.getTrialPoint().getValue(var)
        return self.func(x)
    
def fun_to_opt(x):
    return (x-3)*(x-3)+1

variables = ArrayList()
xvar = Variable( "x", 0.0, -25.0, 25.0 )
variables.add(xvar)

scorer = MyScorer(variables,fun_to_opt)
n_iterations = 50
maxSolutionStopper = SolveStopperFactory.maxEvaluationsStopper(n_iterations)
# solver = Solver(SimplexSearchAlgorithm(),maxSolutionStopper)
solver = Solver(RandomSearch(),maxSolutionStopper)
problem = ProblemFactory.getInverseSquareMinimizerProblem(variables,scorer,0.00001)
delta_hint = InitialDelta()
problem.addHint(delta_hint)
solver.solve(problem)

trial = solver.getScoreBoard().getBestSolution()
result = trial.getTrialPoint().getValue(xvar)

print 'Result is', result


