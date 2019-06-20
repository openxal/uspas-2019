#==================================================================
# This script calculates dispersion function along HEBT
#==================================================================

import sys
import math
import types
import time

from jarray import *
from java.lang import *
from java.util import *
from java.io import *

from Jama import Matrix

from xal.extension.widgets.plot import BasicGraphData, FunctionGraphsJPanel

from xal.smf.impl import Marker, Quadrupole, RfGap, BPM, RfCavity, ProfileMonitor
from xal.model.probe import ParticleProbe
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

from xal.tools.beam.calc import CalculationsOnBeams

from xal.smf.data import XMLDataManager

#===============================================================
#              MAIN PROGRAM
#===============================================================
accl = XMLDataManager.loadDefaultAccelerator()
accSeq = accl.findSequence("HEBT")

bpms = accSeq.getAllNodesOfType("BPM")

print "==========================="

scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)


particle_tracker = AlgorithmFactory.createParticleTracker(accSeq)
particle_tracker.setRfGapPhaseCalculation(True)
probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

#----- bpm_phases_arr = [bpm, position, phase]
bpm_phases_init_arr = []
for bpm in bpms:
	#print "bpm=",bpm.getId()
	state = traj.stateForElement(bpm.getId())
	gamma = state.getGamma()
	phase_matr = state.getResponseMatrix()
	disp_func = phase_matr.getElem(0,5)/gamma**2
	print "bpm = ",bpm.getId()," x[mm] = %5.2f "%(state.getPhaseCoordinates().getx()*1000.)," D[m] = %6.3f "%disp_func 
	
	
print "=========================================="

probe = ProbeFactory.createParticleProbe(accSeq,particle_tracker)

delta_eKin = 1.0e+6
eKin_init = probe.getKineticEnergy()
gamma = probe.getGamma()
beta = probe.getBeta()
mass = probe.getSpeciesRestEnergy()
dp_p =  (delta_eKin/(eKin_init+mass))/(beta**2)

probe.setKineticEnergy(eKin_init+delta_eKin)

scenario.resync()	
scenario.setProbe(probe)
scenario.run()

traj = scenario.getProbe().getTrajectory()

#----- bpm_phases_arr = [bpm, position, phase]
bpm_phases_init_arr = []
for bpm in bpms:
	#print "bpm=",bpm.getId()
	state = traj.stateForElement(bpm.getId())
	disp_func = state.getPhaseCoordinates().getx()/dp_p
	delta_x = state.getPhaseCoordinates().getx()*1000.
	print "bpm = ",bpm.getId()," x[mm] = %5.2f "%(delta_x)," D[m] = %6.3f "%(disp_func)
	
print "Done."