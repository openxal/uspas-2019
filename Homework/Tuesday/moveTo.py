import sys
import math
import types
import time

from jarray import *
from java.lang import *
from java.util import *
from java.io import *

from xal.ca import ChannelFactory

from xal.ca import Monitor
from xal.ca import Channel



pv_x_set = ChannelFactory.defaultFactory().getChannel("Slit:x_set")
pv_x_set.connectAndWait(1.5)
pv_x = ChannelFactory.defaultFactory().getChannel("Slit:x")
pv_x.connectAndWait(1.5)

if len(sys.argv) >1 :
    dest = float(sys.argv[1])
    pv_x_set.putVal(dest)
else:
   rb=pv_x.getValDbl()
   print 'Current position',rb
   exit(0)

while 2 > 1:
    rb=pv_x.getValDbl()
    print 'Current position',rb
    if abs(rb-dest) < 0.1:
        print 'Destination ', dest, ' reached'
        break
    time.sleep(0.5)

exit(0)

