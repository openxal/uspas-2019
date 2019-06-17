import sys
import math
import types
import time

from jarray import *
from java.lang import *
from java.util import *
from java.io import *

from xal.ca import BatchGetValueRequest, ChannelFactory

from xal.ca import Monitor
from xal.ca import Channel
from xal.ca import IEventSinkValue



pv = ChannelFactory.defaultFactory().getChannel("MEBT_Mag:PS_QH01:B_Set")
pv.connectAndWait(1.5)
print "ca=",pv.channelName()," val=",pv.getValDbl()
pv.putVal(0.005)

class ParameterListener(IEventSinkValue):
        def __init__(self):
                self.a = 1

        def eventValue(self,record,pv):
                val = record.doubleValue()
                print "debug pv=",pv.channelName()," val=",val


monitor = pv.addMonitorValue(ParameterListener(), Monitor.VALUE)
