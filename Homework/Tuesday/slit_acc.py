#
# This ia a Jython script
#
# This slit test Virtual Accelerator
# PV names
# Slit:x_set     - set slit position [mm]
# Slit:x         - read back slit position [mm]
# Slit:fc_signal - Faraday Cup (FC) signal [C]

import sys
import math
import types
import time
import random

from xal.ca import ChannelFactory
from xal.ca import Monitor
from xal.ca import Channel
from xal.ca import IEventSinkValue

channelFactory = ChannelFactory.newServerFactory()

class GaussDistrubution:
	def __init__(self,N,center,sigma,rel_base):
		"""
		g(x) = (N/(sigma*sqrt(2*pi)))*exp(-(x-center)**2/(2*sigma**2))
		"""
		self.N = N
		self.center = center
		self.sigma = sigma
		self.base = rel_base*N/(sigma*math.sqrt(2*math.pi))
		self.distance_to_slit_edge = 32.0 # mm
		self.beam_pipe_radius = 15 # mm
		self.step = sigma/100
		
	def getValue(self,x):
		z = (x-self.center)**2/(2*self.sigma**2)
		if(z > 36.): return 0.
		g = self.N/(self.sigma*math.sqrt(2*math.pi))*math.exp(-z)
		return g

	def getIntegral(self,x0,x1):
		"""
		Calculates integral from x0 to x1
		"""
		if(x0 >= x1): return 0.
		n_step = int(abs(x1 - x0)/self.step) + 1
		if(n_step < 3): n_step = 3
		step = (x1 - x0)/(n_step - 1)
		sum_val = (self.getValue(x1) + self.getValue(x0))/2
		for ind in range(1,n_step - 1):
			x = x0 + ind*step
			sum_val += self.getValue(x)
		sum_val *= step
		return sum_val
		
	def getFC_Signal(self,x,width,noise_level,base_noise_level):
		"""
		Returns Faraday Cup signal
		"""
		if(x > - self.beam_pipe_radius):
			signal = self.getIntegral(x-width/2,x+width/2)*(1.0+random.gauss(0.,noise_level))
		else:
			dist = x + self.distance_to_slit_edge
			signal = self.getIntegral(dist,self.beam_pipe_radius)*(1.0+random.gauss(0.,noise_level))
		base = self.base*(1.0+random.gauss(0.,base_noise_level))
		return (signal + base)	

#---- width of the slit [mm]
width = 0.5

#----- total change in the macro-pulse
Q_total = 1.5e+11*1.602177e-19   # C
#print "Q total [C] =",Q_total

#--- these are beam parameters: position of the center and the RMS size [mm]
center = 2.5
sigma = 2.36

#----average noise base level relative to peak of the Gaussian
relative_base = 2/100.

#----- noise levels for base and beam (gaussian) signals [%/100]
noise_level = 2.5/100
base_noise_level = 40./100

gauss_distr = GaussDistrubution(Q_total,center,sigma,relative_base)

pv_x = channelFactory.getChannel("Slit:x")
pv_x.connectAndWait(3.0)
pv_x.putVal(0.)

pv_x_set = channelFactory.getChannel("Slit:x_set")
pv_x_set.connectAndWait(3.0)
pv_x_set.putVal(0.)

pv_signal = channelFactory.getChannel("Slit:fc_signal")
pv_signal.connectAndWait(3.0)
pv_signal.putVal(gauss_distr.getFC_Signal(0.,width,noise_level,base_noise_level))

pv_x_set.setSettable(True)
pv_signal.setSettable(True)

def moveSlit_one_step(pv_x_set,pv_x):
	"""
	Moves readback to the set value.
	Parameters are for the actuator mechanics.
	It returns True if readback value in 
	"""
	set_accuracy = 0.02 #mm
	move_step = 0.2 #mm
	time_sleep = 0.2 #mm
	x_set = pv_x_set.getValDbl()
	x_rb = pv_x.getValDbl()
	if(abs(x_set - x_rb) <= move_step):
		if(abs(x_set - x_rb) >= set_accuracy):
			pv_x.putVal(x_set + random.gauss(0.,set_accuracy))
		return True
	else:
		if(x_rb > x_set):
			pv_x.putVal(x_rb - move_step + random.gauss(0.,set_accuracy))
		else:
			pv_x.putVal(x_rb + move_step + random.gauss(0.,set_accuracy))
	time.sleep(time_sleep)
	return False

def generateSignal(pv_x,pv_signal):
	x = pv_x.getValDbl()
	signal = gauss_distr.getFC_Signal(x,width,noise_level,base_noise_level)
	pv_signal.putVal(signal)
	return signal

count = 0
#-----------start VA loop
while(1 < 2):
	time.sleep(0.2)
	while(not moveSlit_one_step(pv_x_set,pv_x)):
		signal = generateSignal(pv_x,pv_signal)
	signal = generateSignal(pv_x,pv_signal)
	x = pv_x.getValDbl()
	#print "count=",count," x=",x," signal=",signal
	count += 1

