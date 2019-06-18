#==================================================================
# This script will create a 3 kicks bump in the vertical trajectory.
# The kick value from the 1st dipole corrector is user defined 
# value. The last two kicks are defined by closed bump condition.
# It means the initial orbit should be restored after the 3-kickers
# area. The bump closing procedure is done by the OpenXAL optimizer.
# This script is working with VA.
#==================================================================

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

from xal.smf.data import XMLDataManager

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

accl = XMLDataManager.loadDefaultAccelerator()
accSeq = accl.findSequence("MEBT")

print "========================================="

dcvs = accSeq.getAllNodesOfType("DCV")

#---- we will use only 3 first correctors
dcvs = dcvs[:3]
for dcv in dcvs:
	print "debug dcv=",dcv.getId()," pos[m]=",accSeq.getPosition(dcv)," effLength[m]=",dcv.getEffLength()
	#dcv.setField(0.)

time.sleep(1.5)

pos_last_dcv =  accSeq.getPosition(dcvs[len(dcvs)-1])

print "========================================="
bpms_all = accSeq.getAllNodesOfType("BPM")

#--- we will use only BPMs after 3 kicks dipoles
bpms = []
for bpm in bpms_all:
	pos_bpm = accSeq.getPosition(bpm)
	if(pos_bpm > pos_last_dcv):
		bpms.append(bpm)

for bpm in bpms:
	print "debug bpm=",bpm.getId()," pos[m]=",accSeq.getPosition(bpm)
	
print "========================================="	
	
field_start = 0.01
dcvs[0].setField(0.01)

time.sleep(1.5)


def getBPM_Y_Arr(bpms):
	bpm_y_arr = []
	for bpm in bpms:
		y = bpm.getYAvg()
		bpm_y_arr.append(y)
	return bpm_y_arr
	
#-------------------------------------
# Scorer Interface implementation
#-------------------------------------
class OrbitScorer(Scorer):
	""" 
	Calculate the avreage deviation of the orbit (at BPM positions)
	from the center by using data from trial_point variables
	variables is Java's ArrayList()
	"""
	def __init__(self,bpms,dcvs,variables):
		self.bpms = bpms
		self.dcvs = dcvs
		self.variables = variables
		#--------------------------------
		self.diff2_min = Double.MAX_VALUE
		self.count = 0
		
	def score(self,trial,variables):
		self.count += 1
		#---------------------------------
		diff2_dcv = 0.
		limit = 0.98
		penalty_rate = 1000000.0
		for dcv_ind in range(self.variables.size()):
			var = self.variables.get(dcv_ind)
			field =  trial.getTrialPoint().getValue(var)
			self.dcvs[dcv_ind].setField(field)			
			min_field = self.variables.get(dcv_ind).getLowerLimit()*limit
			max_field = self.variables.get(dcv_ind).getUpperLimit()*limit
			if(field < min_field):
				diff2_dcv = penalty_rate*(field - min_field)**2
			if(field > max_field):
				diff2_dcv = penalty_rate*(field - max_field)**2
		diff2_dcv /= self.variables.size()
		#-----------------------------------------------
		time.sleep(1.5)
		#-----------------------------------------------
		bpm_y_arr = getBPM_Y_Arr(self.bpms)
		diff2 = 0.
		for y in bpm_y_arr:
			diff2 += y**2
		diff2 /= len(bpm_y_arr)
		#--------------------------------------
		diff2 += diff2_dcv
		if(self.diff2_min > diff2):
			self.diff2_min = diff2
			print "debug solver count = ",self.count,"  sqrt(diff2)= %12.5g "%math.sqrt(diff2)
		return diff2

	def getDCV_Field_Arr_for_Trial(self,trial):
		#------ return dcv field array for the trial point 
		field_arr = []
		for dcv_ind in range(self.variables.size()):
			var = self.variables.get(dcv_ind)
			field =  trial.getTrialPoint().getValue(var)
			field_arr.append(field)
		return field_arr

#-------------------------------------------------------
#---- in optimization we will use only 2 last correctors
#-------------------------------------------------------
dcvs_opt = dcvs[1:]

#---- Initial step in parameters. During optimization
#---- these steps will be reduced inside the optimizer. 
delta_hint = InitialDelta()

#---- optimizing variabes
variables = ArrayList()

field_max =  0.012
field_min = -0.012

field_step = (field_max - field_min)/30

for dcv_ind in range(len(dcvs_opt)):
	dcv = dcvs_opt[dcv_ind]
	field = dcv.getField()
	var = Variable(dcv.getId(),field, field_min, field_max)
	variables.add(var)
	delta_hint.addInitialDelta(var,field_step)
	
scorer = OrbitScorer(bpms,dcvs_opt,variables)

n_iterations = 50
maxSolutionStopper = SolveStopperFactory.maxEvaluationsStopper(n_iterations) 
solver = Solver(SimplexSearchAlgorithm(),maxSolutionStopper)
#solver = Solver(RandomShrinkSearch(),maxSolutionStopper)
problem = ProblemFactory.getInverseSquareMinimizerProblem(variables,scorer,0.00000001)
problem.addHint(delta_hint)
solver.solve(problem)

#------- get optimization results
trial = solver.getScoreBoard().getBestSolution()
dcv_field_arr = scorer.getDCV_Field_Arr_for_Trial(trial)

print "========== put new dcv fields into VA ======"
print "dcv=",dcvs[0].getId()," field= %+8.6f "%dcvs[0].getField()
#---- send all results to EPICS
for dcv_ind in range(len(dcvs_opt)):
	dcv = dcvs_opt[dcv_ind]
	field = dcv_field_arr[dcv_ind]
	dcv.setField(field)
	print "dcv=",dcv.getId()," field= %+8.6f "%field
	
print "Done."

sys.exit()
	