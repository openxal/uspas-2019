
from xal.smf.data import XMLDataManager
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

accl = XMLDataManager.acceleratorWithPath('/home/student/lib/openxal/site/optics/production/main.xal')
seq = accl.findSequence("MEBT")
scenario = Scenario.newScenarioFor(seq)
tracker = AlgorithmFactory.createEnvTrackerAdapt(seq)
probe = ProbeFactory.getEnvelopeProbe(seq,tracker)
scenario.setProbe(probe)
scenario.run()
traj = scenario.getProbe().getTrajectory()

xsz = map(lambda s: (s.getPosition(),s.twissParameters()[0].getEnvelopeRadius()),traj.iterator())
print(list(xsz))

