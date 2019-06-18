#
# This script will v correct the orbit without any knowledge 
# of the lattice. It will use BPM signals to flatten the orbit
# and the dipole correctors as independent parameters.
#

import sys
import math
import types
import time

from jarray import *
from java.lang import *
from java.util import *
from java.io import *
from java.util import ArrayList

from xal.smf.data import XMLDataManager
from xal.smf import AcceleratorSeqCombo

from xal.smf.impl import BPM, VDipoleCorr

#----------------------------------------
# Calsses of OpenXAL Solver Package 
#----------------------------------------
from xal.extension.solver import Scorer
from xal.extension.solver import Trial
from xal.extension.solver import Variable
from xal.extension.solver import Stopper
from xal.extension.solver import SolveStopperFactory
from xal.extension.solver import ProblemFactory
from xal.extension.solver import Solver
from xal.extension.solver import Problem
from xal.extension.solver.algorithm import SimplexSearchAlgorithm
from xal.extension.solver.algorithm import RandomShrinkSearch
from xal.extension.solver.hint import Hint
from xal.extension.solver.hint import InitialDelta

#===============================================================
#              MAIN PROGRAM
#===============================================================
# read the accelerator & make the sequence
accl = XMLDataManager.loadDefaultAccelerator()
"""
#====== Let's construct accelerator ===========
ccl1 = accl.getSequence("CCL1")
ccl2 = accl.getSequence("CCL2")
ccl3 = accl.getSequence("CCL3")
ccl4 = accl.getSequence("CCL4")


#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [ccl1,ccl2,ccl3,ccl4])
"""
DTL1 = accl.getSequence("DTL1")
DTL2 = accl.getSequence("DTL2")
DTL3 = accl.getSequence("DTL3")
DTL4 = accl.getSequence("DTL4")
DTL5 = accl.getSequence("DTL5")
DTL6 = accl.getSequence("DTL6")

#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [DTL1,DTL2,DTL3,DTL4,DTL5,DTL6])

	
bpms = accSeq.getAllNodesOfType("BPM")
for bpm in bpms:
	print "debug bpm=",bpm.getId()," pos[m]=",accSeq.getPosition(bpm)
	
print "=============================================="

dchs = accSeq.getAllNodesOfType("DCH")
for dch in dchs:
	print "debug dch=",dch.getId(),"  pos[m]=",accSeq.getPosition(dch)
	#---- initial set all dch to zero as a test
	#---- dch.setField(0.)
	
print "=============================================="

def getBPM_X_Arr(bpms):
	bpm_x_arr = []
	for bpm in bpms:
		x = bpm.getXAvg()
		bpm_x_arr.append(x)
	return bpm_x_arr

bpm_x_init_arr = getBPM_X_Arr(bpms)
st = " Initial X[mm] = "
for x_val in bpm_x_init_arr:
	st += " %+6.2f "%x_val
print st
print "=========="

#-------------------------------------------------------------------
# Start of orbit optimization ======================================
#-------------------------------------------------------------------

#-------------------------------------
# Scorer Interface implementation
#-------------------------------------
class OrbitScorer(Scorer):
	""" 
	Calculate the avreage deviation of the orbit (at BPM positions)
	from the center by using data from trial_point variables
	variables is Java's ArrayList()
	"""
	def __init__(self,bpms,dchs,variables):
		self.bpms = bpms
		self.dchs = dchs
		self.variables = variables
		#--------------------------------
		self.diff2_min = Double.MAX_VALUE
		self.count = 0
		
	def score(self,trial,variables):
		self.count += 1
		#---------------------------------
		diff2_dch = 0.
		limit = 0.98
		penalty_rate = 1000000.0
		for dch_ind in range(self.variables.size()):
			var = self.variables.get(dch_ind)
			field =  trial.getTrialPoint().getValue(var)
			self.dchs[dch_ind].setField(field)			
			min_field = self.variables.get(dch_ind).getLowerLimit()*limit
			max_field = self.variables.get(dch_ind).getUpperLimit()*limit
			if(field < min_field):
				diff2_dch = penalty_rate*(field - min_field)**2
			if(field > max_field):
				diff2_dch = penalty_rate*(field - max_field)**2
		diff2_dch /= self.variables.size()
		#-----------------------------------------------
		time.sleep(1.5)
		#-----------------------------------------------
		bpm_x_arr = getBPM_X_Arr(self.bpms)
		diff2 = 0.
		for x in bpm_x_arr:
			diff2 += x**2
		diff2 /= len(bpm_x_arr)
		#--------------------------------------
		diff2 += diff2_dch
		if(self.diff2_min > diff2):
			self.diff2_min = diff2
			print "debug solver count = ",self.count,"  sqrt(diff2)= %12.5g "%math.sqrt(diff2)
		return diff2

	def getDCH_Field_Arr_for_Trial(self,trial):
		#------ return dch field array for the trial point 
		field_arr = []
		for dch_ind in range(self.variables.size()):
			var = self.variables.get(dch_ind)
			field =  trial.getTrialPoint().getValue(var)
			field_arr.append(field)
		return field_arr

#---- Initial step in parameters. During optimization
#---- these steps will be reduced inside the optimizer. 
delta_hint = InitialDelta()

#---- optimizing variabes
variables = ArrayList()

field_max =  0.012
field_min = -0.012

field_step = (field_max - field_min)/30

for dch_ind in range(len(dchs)):
	dch = dchs[dch_ind]
	field = dch.getField()
	var = Variable(dch.getId(),field, field_min, field_max)
	variables.add(var)
	delta_hint.addInitialDelta(var,field_step)
	
scorer = OrbitScorer(bpms,dchs,variables)

n_iterations = 200
maxSolutionStopper = SolveStopperFactory.maxEvaluationsStopper(n_iterations) 
#solver = Solver(SimplexSearchAlgorithm(),maxSolutionStopper)
solver = Solver(RandomShrinkSearch(),maxSolutionStopper)
problem = ProblemFactory.getInverseSquareMinimizerProblem(variables,scorer,0.00000001)
problem.addHint(delta_hint)
solver.solve(problem)

#------- get optimization results
trial = solver.getScoreBoard().getBestSolution()
dch_field_arr = scorer.getDCH_Field_Arr_for_Trial(trial)

print "========== put new dch fields into VA ======"
#---- send all results to EPICS
for dch_ind in range(len(dchs)):
	dch = dchs[dch_ind]
	field = dch_field_arr [dch_ind]
	dch.setField(field)
	print "dch=",dch.getId()," field= %+8.6f "%field
	
print "Done."

sys.exit()
