#
# This script measures the BPMs' response to the DCV in CCL and 
# optimizes the DCV fields to minimize the vertical orbit offsets
# at the BPMs' locations
#
# This script performs measurements of the response matrix and 
# the orbit optimization
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
from xal.extension.solver.hint import Hint
from xal.extension.solver.hint import InitialDelta

#===============================================================
#              MAIN PROGRAM
#===============================================================
# read the accelerator & make the sequence
accl = XMLDataManager.loadDefaultAccelerator()

#====== Let's construct accelerator ===========
ccl1 = accl.getSequence("CCL1")
ccl2 = accl.getSequence("CCL2")
ccl3 = accl.getSequence("CCL3")
ccl4 = accl.getSequence("CCL4")


#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [ccl1,ccl2,ccl3,ccl4])

	
bpms = accSeq.getAllNodesOfType("BPM")
for bpm in bpms:
	print "debug bpm=",bpm.getId()," pos[m]=",accSeq.getPosition(bpm)
	
print "=============================================="

dcvs = accSeq.getAllNodesOfType("DCV")
for dcv in dcvs:
	print "debug dcv=",dcv.getId(),"  pos[m]=",accSeq.getPosition(dcv)
	#---- initial set all dcv to zero as a test
	#dcv.setField(0.)
	
print "=============================================="	
print "Let's measure the response matrix."

def getBPM_Y_Arr(bpms):
	bpm_y_arr = []
	for bpm in bpms:
		y = bpm.getYAvg()
		bpm_y_arr.append(y)
	return bpm_y_arr

bpm_y_init_arr = getBPM_Y_Arr(bpms)
st = " Initial Y[mm] = "
for y_val in bpm_y_init_arr:
	st += " %+6.2f "%y_val
print st
print "=========="

dcv_field_init_arr = []
for dcv in dcvs:
	field = dcv.getField()
	dcv_field_init_arr.append(field)
	
field_max =  0.012
field_min = -0.012

#----- let's create the response matrix
#----- response_matrix = [[bpm_dy/dB[mm/T],...]] with dcv index
response_matrix = []

for dcv_ind in range(len(dcvs)):
	dcv = dcvs[dcv_ind]
	dcv.setField(field_max)
	time.sleep(2.0)
	bpm_y_max_arr = getBPM_Y_Arr(bpms)
	dcv.setField(field_min)
	time.sleep(2.0)
	bpm_y_min_arr = getBPM_Y_Arr(bpms)
	bpm_coeff_arr = []
	st = " = [ "
	for bpm_ind in range(len(bpms)):
		bpm_coeff = (bpm_y_max_arr[bpm_ind] - bpm_y_min_arr[bpm_ind])/(field_max-field_min)
		bpm_coeff_arr.append(bpm_coeff)
		if(bpm_ind != 0): st += ","
		st += " %+6.1f "%bpm_coeff
	st += " ] "
	print "dcv=",dcv.getId(),st
	#---- restore the initial field
	dcv.setField(dcv_field_init_arr[dcv_ind])
	response_matrix.append(bpm_coeff_arr)

"""
response_matrix.append([ +0.0 ,   +5.7 ,   -35.4 ,  +111.9,    -15.5 ,  +114.2 ,  -130.1,    -27.7 ,  +201.4 ,   +78.9])
response_matrix.append([ +0.0 ,   +0.0 ,  -119.6 ,   -55.4,   -144.9 ,   -55.0 ,  -104.4,   -186.3 ,  +128.7 ,  +273.5]) 
response_matrix.append([ +0.0 ,   +0.0 ,  +106.2 ,   -17.5,   +114.4 ,   -19.0 ,  +144.2,   +151.1 ,  -199.2 ,  -241.8]) 
response_matrix.append([ +0.0 ,   +0.0 ,   +80.3 ,  +130.9,   +117.3 ,  +132.3 ,    -2.3,   +145.2 ,   +32.8 ,  -185.3]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,   -160.1 ,   -61.9 ,  -114.2,   -205.8 ,  +140.4 ,  +301.8]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,    -82.7 ,  -133.2 ,   +41.1,    -99.9 ,   -85.4 ,  +114.2]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,   +148.3 ,   +87.6 ,   +76.0,   +188.7 ,   -83.0 ,  -267.2]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,    +76.8 ,  +147.2 ,   -61.5,    +91.2 ,  +116.1 ,   -96.4]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,   -97.2,   -184.8 ,  +117.6 ,  +269.2]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,   +35.6,   -109.6 ,   -78.7 ,  +128.5]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,  +160.4,   +155.0 ,  -224.2 ,  -252.1]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,   +73.3,   +175.1 ,   -81.5 ,  -248.8]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,    +72.8 ,  +110.8 ,   -73.1]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,     +0.0 ,  +252.6 ,  +132.8]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,     +0.0 ,  +281.2 ,  +278.1]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,     +0.0 ,  +184.0 ,  +298.0]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +4.0 ,  +181.4]) 
response_matrix.append([ +0.0 ,   +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +0.0,     +0.0 ,    +0.0 ,    +4.0]) 
"""

#---------------------------------------------
# The response matrix is ready
#---------------------------------------------

#-------------------------------------------------------------------
# Start of orbit optimization ======================================
#-------------------------------------------------------------------

#-------------------------------------
# Scorer Interface implementation
#-------------------------------------
class OrbitScorer(Scorer):
	""" 
	Calculate the difference between predicted inverse change in the orbit
	and the initial BPMs vertical values
	variables is Java's ArrayList()
	"""
	def __init__(self, bpm_y_init_arr, dcv_field_init_arr, response_matrix, variables):
		self.bpm_y_init_arr = bpm_y_init_arr
		self.dcv_field_init_arr = dcv_field_init_arr
		self.response_matrix = response_matrix
		self.variables = variables
		self.diff2_min = Double.MAX_VALUE
		self.count = 0
		
	def score(self,trial,variables):
		self.count += 1
		#-----------------------------------------------
		bpm_y_arr = self.bpm_y_init_arr[:]
		n_bpm = len(bpm_y_arr)
		field_arr = []
		for dcv_ind in range(self.variables.size()):
			var = self.variables.get(dcv_ind)
			field =  trial.getTrialPoint().getValue(var)
			delta_field = field - self.dcv_field_init_arr[dcv_ind]
			bpm_coeff_arr =  response_matrix[dcv_ind]
			for bpm_ind in range(n_bpm):
				bpm_y_arr[bpm_ind] += delta_field*bpm_coeff_arr[bpm_ind]
			field_arr.append(field)
		#---- make score value
		diff2 = 0.
		for bpm_ind in range(n_bpm):
			diff2 += bpm_y_arr[bpm_ind]**2
		diff2 /= n_bpm
		diff2_dcv = 0.
		limit = 0.98
		penalty_rate = 1000000.0
		for dcv_ind in range(self.variables.size()):
			field =  field_arr[dcv_ind]
			min_field = self.variables.get(dcv_ind).getLowerLimit()*limit
			max_field = self.variables.get(dcv_ind).getUpperLimit()*limit
			if(field < min_field):
				diff2_dcv = penalty_rate*(field - min_field)**2
			if(field > max_field):
				diff2_dcv = penalty_rate*(field - max_field)**2
		diff2_dcv /= self.variables.size()
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

#---- Initial step in parameters. During optimization
#---- these steps will be reduced inside the optimizer. 
delta_hint = InitialDelta()

#---- optimizing variabes
variables = ArrayList()

field_max =  0.012
field_min = -0.012

field_step = (field_max - field_min)/30

for dcv_ind in range(len(dcvs)):
	dcv = dcvs[dcv_ind]
	field = dcv_field_init_arr[dcv_ind]
	var = Variable(dcv.getId(),field, field_min, field_max)
	variables.add(var)
	delta_hint.addInitialDelta(var,field_step)
	
scorer = OrbitScorer(bpm_y_init_arr, dcv_field_init_arr, response_matrix, variables)

n_iterations = 1000
maxSolutionStopper = SolveStopperFactory.maxEvaluationsStopper(n_iterations) 
solver = Solver(SimplexSearchAlgorithm(),maxSolutionStopper)
problem = ProblemFactory.getInverseSquareMinimizerProblem(variables,scorer,0.00000001)
problem.addHint(delta_hint)
solver.solve(problem)

#------- get optimization results
trial = solver.getScoreBoard().getBestSolution()
dcv_field_arr = scorer.getDCV_Field_Arr_for_Trial(trial)

print "========== debug new dcv fields ======"
#---- send all results to EPICS
for dcv_ind in range(len(dcvs)):
	dcv = dcvs[dcv_ind]
	field = dcv_field_arr [dcv_ind]
	dcv.setField(field)
	print "dcv=",dcv.getId()," field= %+8.6f "%field
	
print "Done."

sys.exit()
