#
# This script changes the phases of all RF cavities and transports the beam
# through the SCL recording the BPM phase differencies.
# It does it by calculating cavities phase shift as a change in 
# the arrival time of the probe particle.
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

#----- all cavities phase shift
cav_pahse_shift = 1.0

#---- let's create a dictionary with arrival time at the 1st gap at each cavity
#---- for design optics

#-------- rfgap_arrival_tm_dict[gapId] = arrival_time
rfgap_arrival_tm_dict = {}
for state_ind in range(traj.numStates()):
	state = traj.stateWithIndex(state_ind)
	if(state.getElementId().find("Rg01") > 0):
		tm = state.getTime()
		rfgap_arrival_tm_dict[state.getElementId()] = tm

#----------------------------------------------------------
#   Now we add 1 deg (RF frequency) to all RF cavities
#----------------------------------------------------------
for cav in cavs:
	cav_phase =  cav.getDfltCavPhase()
	cav.setDfltCavPhase(cav_phase + 1.0)

#=================================================
# Let's set up all new phases for SCL cavities

for rfgap in rf_gaps:
	probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)
	
	scenario.setStopElementId(rfgap.getId())
	scenario.setIncludeStopElement(False)
	
	scenario.resync()
	scenario.setProbe(probe)
	scenario.run()
	
	traj = scenario.getProbe().getTrajectory()
	tm = traj.finalState().getTime()
	
	tm_ini = rfgap_arrival_tm_dict[rfgap.getId()]
	
	delta_t = tm - tm_ini
	
	delta_phase = 360.*805.e+6*delta_t
	
	cav = rfgap.getParent()
	
	
	cav.setDfltCavPhase(cav.getDfltCavPhase() + delta_phase)
	
	#print "debug cav=",cav.getId()," elem=",state.getElementId()," delta_phase=",delta_phase
	

scenario.unsetStopNode()

#=======================================================================
# Now all RF phases are set, and we repeat the bpm phases calculations 
#=======================================================================

probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

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
	print "bpm=",bpm.getId()," pos[m]=",pos," bpm phase diff.[deg] =",phase_diff

print "Done."