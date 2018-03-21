# Main Project 2 File
# Created by Brian Pickens, Nathan Chrisman, and Kelsie Crawford
# Runs the functions defined in simulation.py

import numpy as np
import parameters as pm
import matplotlib.pyplot as plt
import simulation as sim

# Plotting Preferences
plotOrbits = True
plotMomentum = False # <-- Must be the same value (plotted together)
plotEnergy = False   # <-- Must be the same value (plotted together)
plotEnergyTot = False
plotEnergyLog = False
plotProbsHist = False
plotProbsP4 = False
plotWobbleTop = False
plotWobbleCoM = False
plotWobbleRad = False
plotWobbleObs = False

# Name of saved figure
figName = "Default.pdf"

# Calculate new time steps - these override the simulation.py values.
sim.timeDaySim = 24 # Days Simulated
sim.timeStep = 0.01 # Step Size
sim.timeMax = int(sim.timeDaySim / sim.timeStep)

# --- Main two calls to simulate system ---

# Create the system
systemTRAPPIST_1 = sim.createSystem("TRAPPIST-1")

# Calculate orbit data, which are stored within the GBody objects
sim.simulateOrbits(systemTRAPPIST_1)

# --- Optional calls to analyze orbital data in some way ---

# Gets array of system momentum calculations
if plotMomentum == True:
    systemMomentum = sim.getSystemMomentum(systemTRAPPIST_1)
else:
    systemMomentum = 0

# Gets array of system energy calculations    
if plotEnergy == True or plotEnergyTot == True or plotEnergyLog == True:
    systemEnergy = sim.getSystemEnergy(systemTRAPPIST_1)
else:
    systemEnergy = 0
    
# Gets array of planet detection data
if plotProbsP4 == True or plotProbsHist == True:
    viewProbs = sim.getViewProbs(systemTRAPPIST_1)
else:
    viewProbs = 0

# Gets array for system star wobble data
if plotWobbleTop == True or plotWobbleCoM == True or plotWobbleRad == True or plotWobbleObs == True:
    starWobble = sim.getStarWobble(systemTRAPPIST_1)
else:
    starWobble = 0

# --- Plotting Function ---
sim.plotEverything(figName,
                   systemTRAPPIST_1, plotOrbits,
                   systemMomentum, plotMomentum,
                   systemEnergy, plotEnergy, plotEnergyTot, plotEnergyLog,
                   viewProbs, plotProbsP4, plotProbsHist,
                   starWobble, plotWobbleTop, plotWobbleCoM, plotWobbleRad, plotWobbleObs)
