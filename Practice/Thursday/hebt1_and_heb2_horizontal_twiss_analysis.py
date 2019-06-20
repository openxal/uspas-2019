#==================================================================
# This script will calculate X Twiss parameters and energy spread
# at the entrance of HEBT by using WS rms sizes
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
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

from xal.tools.beam.calc import CalculationsOnBeams

from xal.smf.data import XMLDataManager

#===============================================================
#              MAIN PROGRAM
#===============================================================

accl = XMLDataManager.loadDefaultAccelerator()
accSeq = accl.findSequence("HEBT")

wss = accSeq.getAllNodesOfType(ProfileMonitor.s_strType)
for ws in wss:
	print "ws=",ws.getId()," pos=",accSeq.getPosition(ws)

print "==========================="

quads = accSeq.getAllNodesOfType(Quadrupole.s_strType)
for quad in quads:
	print "quad=",quad.getId()," B=",quad.getDfltField()

print "==========================="

scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)

tracker = AlgorithmFactory.createEnvTrackerAdapt(accSeq)
tracker.setProbeUpdatePolicy(tracker.UPDATE_ALWAYS)
		
probe = ProbeFactory.getEnvelopeProbe(accSeq,tracker)
#---- peak current in [A]
peak_current = 0.000
probe.setBeamCurrent(peak_current)

#---- initial kinetic energy in [MeV]
eKin_init = probe.getKineticEnergy()/1.0e+6

scenario.setProbe(probe)
scenario.resync()
scenario.run()

traj = scenario.getProbe().getTrajectory()

beam_calculator = CalculationsOnBeams(traj)

for ws in wss:
	state = traj.stateForElement(ws.getId())
	pos = state.getPosition()
	phase_arr = beam_calculator.computeBetatronPhase(state)
	phaseX = phase_arr.getx()*180./math.pi
	phaseY = phase_arr.gety()*180./math.pi
	print "ws=",ws.getId()," pos[m]= %8.3f "%pos," (PhaseX, PhaseY) =  (%6.1f,%6.1f)"%(phaseX,phaseY)

#--------- sizes for each WS and transport matrices
matrx_all_X = []
sizes_X_arr = []
for ws in wss:
	state = traj.stateForElement(ws.getId())
	pos = state.getPosition()
	xRMS_Size = state.twissParameters()[0].getEnvelopeRadius()
	mtrx = state.getResponseMatrix()
	#----------------------
	sizes_X_arr.append(xRMS_Size)
	#--------elements of the transport matrix
	a11 = mtrx.getElem(0,0)
	a12 = mtrx.getElem(0,1)
	a16 = mtrx.getElem(0,5)
	matrx_all_X.append([a11**2, 2*a11*a12,a12**2, a16**2])
	
#=========================================
#  Only x-axis analysis
#=========================================
sigma_rel_err = 0.01

n_ws = len(matrx_all_X)

mMatrX = Matrix(matrx_all_X,n_ws,4)

sigma2Vector = Matrix(n_ws,1)
weightM = Matrix.identity(n_ws,n_ws)

for ind in range(n_ws):
	sigma = sizes_X_arr[ind]
	sigma2Vector.set(ind,0,sigma**2)
	err2 = (2*sigma*sigma*sigma_rel_err)**2
	weightM.set(ind,ind,1.0/err2)
	
#=== mwmMatr = (M^T*W*M) =======
mwmMatr = ((mMatrX.transpose()).times(weightM)).times(mMatrX)

#=== corr2ValMatr = [(M^T*W*M)^-1] * M^T * W * Vector(sigma**2)
corr2ValMatr = (((mwmMatr.inverse()).times(mMatrX.transpose())).times(weightM)).times(sigma2Vector)

corr2ErrMatr = mwmMatr.inverse()

print "========================================"
print "<x^2>            = ","%12.5e"%corr2ValMatr.get(0,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(0,0)))
print "<x*x'>           = ","%12.5e"%corr2ValMatr.get(1,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(1,1)))
print "<x'^2>           = ","%12.5e"%corr2ValMatr.get(2,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(2,2)))
print "<(dbeta/beta)^2> = ","%12.5e"%corr2ValMatr.get(3,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(3,3)))
print "========================================"

x2 = corr2ValMatr.get(0,0)
xxp = corr2ValMatr.get(1,0)
xp2 = corr2ValMatr.get(2,0)
dbb2 = corr2ValMatr.get(3,0)

gam = 1. + eKin_init/938.
dpp = math.sqrt(dbb2)*gam*gam
print "<dp/p> = ","%12.5e" %dpp

print "========================================"
print "x2, xxp, xp2, dbb2, gam, dpp = ",   x2, xxp, xp2, dbb2, gam, dpp
print "========================================"

# Calculate some other parameters at the start of HEBT2
emitX = math.sqrt(x2*xp2-xxp*xxp)
alphX = -xxp/emitX
betaX = x2/emitX
gammaX = xp2/emitX
print "emittance =", emitX, ", apha = ", alphX, ", beta = ", betaX, ", gamma = ", gammaX, ", dp/p =", dpp
print "========================================"

#------------------------------------------
state_init = traj.initialState()
x_rms_ini = state_init.twissParameters()[0].getEnvelopeRadius()

twissX = state_init.twissParameters()[0]
print "Initial Twiss X=",twissX

twissZ = state_init.twissParameters()[2]
print "Initial Twiss Z=",twissZ

print "Done."
