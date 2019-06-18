#==================================================================
# This script will calculate the 3 correctors bump by using 
# the OpenXAL Online Model transport matrices.
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
	dcv.setField(0.)
	
#-------------------------------------------------
# Transport matrix calculations
#-------------------------------------------------

scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_LIVE)

particle_tracker = AlgorithmFactory.createParticleTracker(accSeq)
particle_tracker.setRfGapPhaseCalculation(True)
probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

mtrx_arr = []
for state_ind in range(traj.numStates()):
	state = traj.stateWithIndex(state_ind)
	for dcv in dcvs:
		if(state.getElementId().find(dcv.getId()) >= 0):
			mtrx = state.getResponseMatrix()
			mtrx_arr.append(mtrx)

mtrx_0_1 = mtrx_arr[1].times(mtrx_arr[0].inverse())
mtrx_1_2 = mtrx_arr[2].times(mtrx_arr[1].inverse())

print "========================================="

#--------------------------------
# Dipole corrector kick
# dp_p [radian] = Leff*B/(Brho) , Leff in [m] and B in [T]
# Brho = 3.33564*(momentum [GeV/c]) 
# Brho = 3.33564*sqrt(Tkin*(2m+Tkin)) , Tkin and m in GeV
#--------------------------------

state_init = traj.initialState()
momentum = state_init.getMomentum()/1.0e+9 # GeV/c
Leff_0 = dcvs[0].getEffLength()

#---- kick from corrector #0
field_dcv_0 = 0.01
kick_0 = Leff_0*field_dcv_0/(3.33564*momentum)
print "dcv=",dcvs[0].getId()," B[T]= %+7.6f "%field_dcv_0,"  kick[mrad]=",kick_0*1.0e+3

sys.exit(0)
