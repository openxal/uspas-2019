#
# This script will find synchronous phases of SCL RF cavities
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

from xal.smf.impl import BPM, RfCavity


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

#+++++++++++++++++++++++++++++++++++++++++++++
accSeq = AcceleratorSeqCombo("SEQUENCE", [sclMed,sclHigh])

#----- all RF cavities
cavs = accSeq.getAllNodesOfType(RfCavity.s_strType)
for cav in cavs:
	print "debug cav=",cav.getId()," phase=",cav.getDfltCavPhase(),"  1st gap=", cav.getGaps()[0]," phase=",cav.getGaps()[0].getGapDfltPhase()
	
#----- Online Model for the acc. sequence
scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)

particle_tracker = AlgorithmFactory.createParticleTracker(accSeq)
particle_tracker.setRfGapPhaseCalculation(True)
probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

cav_init_ampl_phase_dict = {}
for cav in cavs:
	cav_phase = cav.getDfltCavPhase()
	cav_amp = cav.getDfltCavAmp()
	cav_init_ampl_phase_dict[cav] = [cav_amp,cav_phase]
	
#switch off all cavities
for cav in cavs:
	cav.setDfltCavAmp(0.)

#=======================================================

def findMaxAccelerationPhase(accSeq,scenario,particle_tracker,cav,min_cav_phase,max_cav_phase,cav_phase_step):
	"""
	this function will find the cavity's phase for maximal acceleration 
	"""
	n_steps = int((max_cav_phase-min_cav_phase)/cav_phase_step) + 1
	max_accel_phase = min_cav_phase
	cav.setDfltCavPhase(min_cav_phase)
	probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)
	scenario.resync()	
	scenario.setProbe(probe)
	scenario.run()
	traj = scenario.getProbe().getTrajectory()	
	eKin_max = traj.finalState().getKineticEnergy()
	for ind in range(n_steps):
		cav_phase = min_cav_phase + ind*cav_phase_step
		cav.setDfltCavPhase(cav_phase)
		probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)
		scenario.resync()	
		scenario.setProbe(probe)
		scenario.run()
		traj = scenario.getProbe().getTrajectory()	
		eKin = traj.finalState().getKineticEnergy()
		if(eKin_max < eKin):
			max_accel_phase = cav_phase
			eKin_max = eKin
	#print "debug cav=",cav.getId()," eKin_max=",eKin_max/1.0e+6
	return max_accel_phase

#=============================================================


#---- start to go over all cavities calculating synch. phase
for cav in cavs:
	cav_amp = cav_init_ampl_phase_dict[cav][0]
	cav_phase = cav_init_ampl_phase_dict[cav][1]
	cav.setDfltCavAmp(cav_amp)
	max_accel_phase = findMaxAccelerationPhase(accSeq,scenario,particle_tracker,cav,-180.,+180.,3.0)
	max_accel_phase = findMaxAccelerationPhase(accSeq,scenario,particle_tracker,cav,max_accel_phase-10.,max_accel_phase+10.,1.0)
	max_accel_phase = findMaxAccelerationPhase(accSeq,scenario,particle_tracker,cav,max_accel_phase-3.,max_accel_phase+3.0,0.1)
	synch_phase = - makePhaseNear(max_accel_phase - cav_phase,0.)
	cav.setDfltCavPhase(cav_phase)
	print "cav=",cav.getId()," synch_phase [deg]= %+6.1f "%synch_phase
	


print "Done."