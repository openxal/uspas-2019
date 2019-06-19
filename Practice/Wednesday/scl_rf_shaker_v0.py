#
# This script changes the start time of the probe and transports the beam
# through the SCL recording the BPM phase differences.
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
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

from xal.smf.impl import BPM, Quadrupole, RfCavity


def makePhaseNear(phase, phase0):
	""" It will add or subtruct any number of 360 deg from phase to get close to phase0 """
	n = int(phase0/360.)
	phase = phase%360.
	min_x = 1.0e+38
	n_min = 0
	for i0 in range(5):
		i = i0 - 3
		d = math.fabs(phase + 360.*(i+n) - phase0)
		if(d < min_x):
			n_min = i
			min_x = d
	return (phase + 360.*(n_min+n))	


#===============================================================
#              MAIN PROGRAM
#===============================================================
# read the accelerator & make the sequence
accl = XMLDataManager.loadDefaultAccelerator()

#====== Let's construct accelerator ===========
sclMed = accl.getSequence("SCLMed")
sclHigh = accl.getSequence("SCLHigh")
hebt1 = accl.getSequence("HEBT1")

#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [sclMed,sclHigh,hebt1])

#----- all RF cavities
cavs = accSeq.getAllNodesOfType(RfCavity.s_strType)
for cav in cavs:
	print "debug cav=",cav.getId()," phase=",cav.getDfltCavPhase(),"  1st gap=", cav.getGaps()[0]," phase=",cav.getGaps()[0].getGapDfltPhase()
	

bpms = accSeq.getAllNodesOfType("BPM")
bpm_names_dict = {}
for bpm in bpms:
	bpm_names_dict[bpm.getId()] = bpm
	#print "debug bpm=",bpm.getId()," pos[m]=",accSeq.getPosition(bpm)," freq.=",bpm.getBPMBucket().getFrequency()

#---- array of all 1st RF gaps
rf_gaps = []
for cav in cavs:
	rf_gaps.append(cav.getGaps()[0])

#----- New Online Model for the acc. sequence
scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)

particle_tracker = AlgorithmFactory.createParticleTracker(accSeq)
particle_tracker.setRfGapPhaseCalculation(True)
probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

"""
#----- bpm_phases_arr = [bpm, position, phase]
bpm_phases_init_arr = []
for bpm in bpms:
	print "bpm=",bpm.getId()
	#state = traj.stateForElement(bpm.getId())
	state = traj.statesForElement(bpm.getId())[0]
	tm = state.getTime()
"""	

#----- bpm_phases_init_arr = [[bpm,position,phase],...]
bpm_phases_init_arr = []
for state_ind in range(traj.numStates()):
	state = traj.stateWithIndex(state_ind)
	if(state.getElementId().find("BPM") >= 0):
		bpm = bpm_names_dict[state.getElementId()]
		tm = state.getTime()
		bpm_phase = 360.0*bpm.getBPMBucket().getFrequency()*1.0e+6*tm
		bpm_phase = makePhaseNear(bpm_phase,0.)
		bpm_phases_init_arr.append([bpm,state.getPosition(),bpm_phase])
		#print "elem Id=",state.getElementId()," pos =",state.getPosition()," time=",state.getTime()

#----------------------------------------------------------
#   Now we add 1 deg (RF frequency) to the particle probe
#----------------------------------------------------------
probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)
tm = probe.getTime() + 1./(360.*805.0e+6)
probe.setTime(tm)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

bpm_phases_new_arr = []
for state_ind in range(traj.numStates()):
	state = traj.stateWithIndex(state_ind)
	if(state.getElementId().find("BPM") >= 0):
		bpm = bpm_names_dict[state.getElementId()]
		tm = state.getTime()
		bpm_phase = 360.0*bpm.getBPMBucket().getFrequency()*1.0e+6*tm
		bpm_phase = makePhaseNear(bpm_phase,0.)
		bpm_phases_new_arr.append([bpm,state.getPosition(),bpm_phase])
		#print "elem Id=",state.getElementId()," pos =",state.getPosition()," time=",state.getTime()
		

for ind in range(len(bpm_phases_init_arr)):
	phase_diff = bpm_phases_new_arr[ind][2] - bpm_phases_init_arr[ind][2]
	phase_diff = makePhaseNear(phase_diff,0.)
	bpm = bpm_phases_new_arr[ind][0]
	pos = bpm_phases_new_arr[ind][1]
	print "bpm=",bpm.getId()," pos[m]=",pos," phase diff. [deg] =",phase_diff

print "Done."