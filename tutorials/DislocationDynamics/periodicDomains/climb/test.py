import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('./')
import pyMoDELib

simulationDir = "./"

ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)

# Microstructure Generation
microstructureGenerator=pyMoDELib.MicrostructureGenerator(ddBase)

spec1=pyMoDELib.ShearLoopIndividualSpecification()
spec1.slipSystemIDs=[0,-1]
spec1.loopRadii=[27.0e-9,27.0e-9]
spec1.loopCenters=np.array([[20.0,0.0,0.0],[0.0,0.0,0.0]])
spec1.loopSides=[10,10]
print(f"calling microstructureGenerator.addShearLoopIndividual(spec1)")
microstructureGenerator.addShearLoopIndividual(spec1)
print(f"finished call microstructureGenerator.addShearLoopIndividual(spec1)")

print(f"calling microstructureGenerator.writeConfigFiles(0)")
microstructureGenerator.writeConfigFiles(0) # write evl_0.txt (optional)
print(f"finished call to microstructureGenerator.writeConfigFiles(0)")

# instantiate/initialize Defective Crystal
print(f"calling defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)")
defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
print(f"finished call to defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)")
print(f"calling defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)")
defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)
print(f"finished call to defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)")

# TODO: run DDD, or read pre-existing configuration from evl_*.txt files
# read using DDconfigIO.read( const size_t& runID )

# retrieve DislocationNetwork, compute displacements
DN=defectiveCrystal.dislocationNetwork()

points = np.array([[0,1,2],[3,4,5.]])

meshSize=100.
lowerMeshSize = 10.
for loopID in DN.loops():
    loop=DN.loops().getRef(loopID)
    meshedLoopVector=loop.meshed(meshSize, lowerMeshSize)
    for meshedLoop in meshedLoopVector:
        disp=meshedLoop.plasticDisplacement(points)
        print(disp)
