# Simulation.py definition file
# Created by Brian Pickens, Nathan Chrisman, and Kelsie Crawford
# Contains all definitions for simulating TRAPPIST-1 system

import numpy as np
import parameters as pm
import matplotlib.pyplot as plt

# Calculate default time steps (for two weeks, just enough for quick 'py.test' testing)
# Overwritten by project02.py
timeDaySim = 14 # Days Simulated
timeStep = 0.01 # Step Size
timeMax = int(timeDaySim / timeStep)

# Begin GBody class definition ----------------------------------------------------------
class GBody():
    # Gravitational Body Object Class
    # To be used for instantiation of TRAPPIST-1 star and planets
    # Designed to be reusable for whatever reason

    # Global variable for counting number of instances
    GBodyCount = 0

    def __init__(self, name):
        # Called when object is created from class, performs initial calculatons

        # Increment number of instanced classes
        GBody.GBodyCount += 1

        # Assign its given name
        self.Name = name

        #Is it the star or one of the planets? Assume star until proven otherwise
        self.isPlanet = False

        # Determine what parameters to draw from and implement
        if "TRAPPIST-1" in name:

            if name == "TRAPPIST-1":
                self.Mass = pm.star['mass']
                self.IsPlanet = False

                # All of the object's movement data is stored here.
                # Dimensions (timeMax, 6) which contains x, y, vx, vy
                self.MovementData = np.zeros([timeMax, 4], dtype=np.float64)

            elif name[10] == "b":
                self.isPlanet = True
                self.pmIndex = 0

            elif name[10] == "c":
                self.isPlanet = True
                self.pmIndex = 1

            elif name[10] == "d":
                self.isPlanet = True
                self.pmIndex = 2

            elif name[10] == "e":
                self.isPlanet = True
                self.pmIndex = 3

            elif name[10] == "f":
                self.isPlanet = True
                self.pmIndex = 4

            elif name[10] == "g":
                self.isPlanet = True
                self.pmIndex = 5

            else:
                print("ERROR: Invalid 'name' For TRAPPIST-1 System detected in GBody.__init__")
                self.isPlanet = False   # Get upset
                self.Mass = 99999999999 # Throw tantrum

            # The planet is known, now calculate and assign initial conditions.
            # Check so we don't apply planet calculations to the star
            if self.isPlanet == True:
                self.Mass = pm.planets['mass'][self.pmIndex]

                self.MovementData = np.zeros([timeMax, 4], dtype=np.float64)

                # Assign initial conditions
                # First, calculate Aphelion and set position to x (all planets begin on x-axis)
                self.Aphelion = pm.planets['semi-major'][self.pmIndex] * (1 + pm.planets['eccentricity'][self.pmIndex])
                self.MovementData[0, 0] = self.Aphelion

                # Then, get velocity at the Aphelion and set to positive y direction (all planets going counter-clockwise)
                self.MovementData[0, 3] = np.sqrt((pm.G_local * pm.star['mass'] / pm.planets['semi-major'][self.pmIndex]) * ((1 - pm.planets['eccentricity'][self.pmIndex]) / (1 + pm.planets['eccentricity'][self.pmIndex])))

                #print(self.MovementData[0, :])
        else:
            print("ERROR: No recognized system name passed to GBody()")

    def totalGravForce(self, systemList, time):
        # Calculates total gravitational force on self object from the list passed into calcNewPosition
        # Called from calcNewPosition() and nowhere else

        totalF = np.zeros([2])

        for i in range(self.GBodyCount):
            if systemList[i].Name != self.Name: # Don't factor self on self in total gravity magnitude

                # print("Calculating Gravity of " + systemList[i].Name + " on " + self.Name)
                # Get distance between the two by calculating x2 - x1 and y2 - y1

                r = systemList[i].MovementData[time, :2] - self.MovementData[time, :2]
                rr = np.sum(r*r)
                totalF += (1 * (pm.G_local * systemList[i].Mass * self.Mass) / rr) * (r / np.sqrt(rr))
                # print(r, rr, totalF)

        return totalF

    def calcNewPosition(self, systemList, time, VVPhase):
        # Takes list of GBody objects ('systemList') with current time and calculates orbits based on their positions
        # This is only done for one time step. The main time loop is in simulateOrbits()

        if time == 0: # Start off Velocity Verlet by calculating FG, ignore afterward
            self.FG = self.totalGravForce(systemList, time)

        # One step of Velocity Verlet
        # Warning: Terrible hacked solution below.
        # VV Algorithm cut in half so that self.FG doesn't calc gravity without position data

        if VVPhase == 1:
            self.vHalf = self.MovementData[time, 2:4] + 0.5 * timeStep * self.FG / self.Mass
            self.MovementData[time + 1, :2] = self.MovementData[time, :2] + timeStep * self.vHalf

        elif VVPhase == 2:
            self.FG = self.totalGravForce(systemList, time + 1)
            self.MovementData[time + 1, 2:4] = self.vHalf + 0.5 * timeStep * self.FG / self.Mass

        else:
            print("ERROR: Invalid VVPhase value")

        # print("Did a VV step for " + self.Name)

        # Doesn't return anything
        return

# End GBody class definition ------------------------------------------------------------




# Begin createSystem() ------------------------------------------------------------------
def createSystem(starName):
    # Input the name of the star, generate system as list of objects
    # System must have parameters specified in parameters.py

    if starName == "TRAPPIST-1":
        system = [GBody("TRAPPIST-1"),
                  GBody("TRAPPIST-1b"), GBody("TRAPPIST-1c"),
                  GBody("TRAPPIST-1d"), GBody("TRAPPIST-1e"),
                  GBody("TRAPPIST-1f"), GBody("TRAPPIST-1g")]

        # Get star's initial (y) velocity
        system[0].MovementData[0, 3] = -1 * (system[1].Mass * system[1].MovementData[0, 3] +
                                             system[2].Mass * system[2].MovementData[0, 3] +
                                             system[3].Mass * system[3].MovementData[0, 3] +
                                             system[4].Mass * system[4].MovementData[0, 3] +
                                             system[5].Mass * system[5].MovementData[0, 3] +
                                             system[6].Mass * system[6].MovementData[0, 3])

        # print(system[0].MovementData[0, :])

    else:
        print("ERROR: No Valid System Specified in createSystem(). Returning 0..")
        return 0

    return system
# End createSystem() --------------------------------------------------------------------




# Begin simulateOrbits() ----------------------------------------------------------------
def simulateOrbits(systemList):
    # Takes GBody object list 'systemList', runs main integration loop on it for total global time steps 'timeMax'

    for i in range(0, timeMax - 1):
    # print("------- ONE TIME STEP -------")

        VVPhase = 1

        for j in range(systemList[0].GBodyCount):
            # Run first part of VV on all bodies and get new position data

            systemList[j].calcNewPosition(systemList, i, VVPhase)

        VVPhase = 2

        for j in range(systemList[0].GBodyCount):
            # Run second part of VV on all bodies and use new position data for new gravity values

            systemList[j].calcNewPosition(systemList, i, VVPhase)

        # Progress Bar for those long waits
        if i % (timeMax / 50) == 0:
            print("Orbit calculations are {}% complete...\r" .format(int(i / timeMax * 100)))

    # Returns nothing
    return
# End simulateOrbits() ------------------------------------------------------------------




# Begin getSystemMomentum() -------------------------------------------------------------
def getSystemMomentum(systemList):
    # Combs through systemList[].MovementData[] calculations and gets the system momentum at all time steps
    # In previous version of simulation.py, this was integrated into the simulateOrbits() process
    # It was separated to allow quicker 1000 days simulations when momentum isn't required

    # Array for tracking system momentum (2nd dimension is [total momentum, change from previous momentum])
    systemMomentum = np.zeros([timeMax, 2])

    print("Calculating system momentum...")

    # Main time loop I
    for i in range(0, timeMax - 1):

        # Body Loop
        for j in range(systemList[0].GBodyCount):

            # Add body's momentum to system total at this time step
            systemMomentum[i, 0] += systemList[j].Mass * systemList[j].MovementData[i, 2] + systemList[j].Mass * systemList[j].MovementData[i, 3]

        # At end of J loop, calculate the change from previous momentum (only if i isnt 0)
        if i != 0:
            systemMomentum[i - 1, 1] = systemMomentum[i, 0] - systemMomentum[i - 1, 0]

    # Calculate final timestep system change in momentum (loop wont do it)
    systemMomentum[timeMax - 1, 1] += systemMomentum[timeMax - 1, 0] - systemMomentum[timeMax - 2, 0]

    # Return systemMomentum array for plotEverything() to process
    return systemMomentum
# End getSystemMomentum() ---------------------------------------------------------------




# Begin getSystemEnergy() ---------------------------------------------------------------
def getSystemEnergy(systemList):
    #Does same exact thing as getSystemMomentum except it calcs [total system energy, change in energy from previous entry, relative error in energy]
    systemEnergy = np.zeros([timeMax, 3])

    print("Calculating system energy...")

    # Main time loop I
    for i in range(0, timeMax):

        # Planet Loop
        for j in range(systemList[0].GBodyCount - 1):
            
            kineticEnergy = 0.5 * systemList[j + 1].Mass * (systemList[j + 1].MovementData[i, 2] + systemList[j + 1].MovementData[i, 3])**2
                
            potentialEnergy = -0.5 * pm.G_local * systemList[j + 1].Mass * systemList[0].Mass / np.sqrt(np.sum(systemList[j + 1].MovementData[i, 0]**2 + systemList[j + 1].MovementData[i, 0]**2))


            # Add J's energy to system total at this time step
            systemEnergy[i, 0] += kineticEnergy + potentialEnergy

        #print(systemEnergy[i, :])
        # At end of J loop, calculate the change from previous energy (only if i isnt 0)
        if i != 0:
            systemEnergy[i - 1, 1] = systemEnergy[i, 0] - systemEnergy[i - 1, 0]


    # Calculate final timestep system change in energy (loop wont do it)
    systemEnergy[timeMax - 1, 1] += systemEnergy[timeMax - 1, 0] + systemEnergy[timeMax - 2, 0]
    
    # print(systemEnergy[0, :])
    
    # Generate the log10(relative error) values for system energy
    for i in range (0, timeMax):
        # Get energy and asbolute value it, log can't accept negative values.
        systemEnergy[i, 2] = np.absolute((systemEnergy[i, 0] - systemEnergy[0, 0]) / systemEnergy[0, 0])
        
        # Don't accept zero values either, just leave it zero.
        if systemEnergy[i, 2] != 0:
            systemEnergy[i, 2] = np.log10(systemEnergy[i, 2])
    
    # Return systemEnergy array for plotEverything() to process
    return systemEnergy

# End getSystemEnergy() -----------------------------------------------------------------




# Begin getViewProbs() ------------------------------------------------------------------
def getViewProbs(systemList, viewerResolution = 1, viewObjectSize = 1000):
    # Scours the systemList[].MovementData[] calculations for the n planet viewing probabilities for all planets.
    # Naturally, the orbits need to be calculated first with simulateOrbits()
    # 'viewerResolution' is viewer's arc minutes of resolution (default 1), 'viewObjectSize' is target size in km (default 1000)

    # Determine maximum distance to resolve object of size 'objectResolution' (km)

    #Convert viewerResolution to degrees
    viewerResolution = viewerResolution / 60 # 1 ArcMinute = 1/60 of a degree

    #Calculate with derived formula Max distance = Object Size * ~57.29 / Viewer Resolution and convert to 1/1000 AU
    maxViewDistance = viewObjectSize * 57.29 / viewerResolution # Get max distance in km
    maxViewDistance = maxViewDistance / pm.au # Convert to AU
    maxViewDistance = maxViewDistance * 1000 # Convert to AU^-3

    print("Calculating planet sighting probabilities at max target distance {} AU^-3..." .format(maxViewDistance))

    # Viewing Probabilty Array.
    # First axis is time steps
    # Second axis contains how many planets any planet can see at once
    # i.e. [[3, 4, 4, 5, 4, 3], ...] --> Planet 4 can see 5 other planets at this time step
    # Doesn't include the star, so -1 count
    viewProbs = np.zeros([timeMax, systemList[0].GBodyCount - 1])

    # Analyze Orbit data for probabilties
    # First, check that the orbit data was calculated at all (Is outer planet x value = 0?)
    if systemList[systemList[0].GBodyCount - 1].MovementData[1, 0] == 0:
        print("ERROR: Orbit Data not calculated before passing system into getViewProbs(). Returning 0..")
        return 0

    # Main time loop I
    for i in range(0, timeMax):

        # Viewing Planet J loop (inefficient but can't think of anything better right now)
        for j in range(systemList[0].GBodyCount - 1):
            planetsSeen = 0 # Haven't seen anything yet

            # Target Planet K Loop
            for k in range(systemList[0].GBodyCount - 1):
                # Get Planet J's and Planet K's distance from each other (+1 so the star isn't factored in)
                # Oh, and don't check to see if you can see yourself
                if j != k:
                    jkDistance = np.sqrt((systemList[j + 1].MovementData[i, 0] - systemList[k + 1].MovementData[i, 0])**2 +
                                     (systemList[j + 1].MovementData[i, 1] - systemList[k + 1].MovementData[i, 1])**2)

                    # Is K in viewing range of J?
                    if jkDistance <= maxViewDistance:
                        planetsSeen += 1

            # At end of K loop,
            viewProbs[i, j] = int(planetsSeen)
    #print(viewProbs)

    # Return the array for processing by plotEverything()
    return viewProbs
# End getViewProbs() --------------------------------------------------------------------




# Begin getStarWobble() -----------------------------------------------------------------
def getStarWobble(systemList):
    # Analyzes star orbital data and returns a 2D array for time steps and
    # [x-positions, y-positions, system center mass x, system center mass y, observer vx, observer vy]
    starWobble = np.zeros([timeMax, 6])
    starWobble[:, :2] = systemList[0].MovementData[:, :2]

    totalPlanetMass = 0
    
    print("Calculating star details...")

    # Get total system planet mass for later
    for i in range(systemList[0].GBodyCount - 1):
        totalPlanetMass += systemList[i + 1].Mass

    # Get this too
    totalSystemMass = totalPlanetMass + systemList[0].Mass

    # Get system center mass at all time steps
    for i in range(0, timeMax):

        # Add up all moments mi*xi and mi*yi
        for j in range(systemList[0].GBodyCount):
            starWobble[i, 2] += systemList[j].Mass * systemList[j].MovementData[i, 0]
            starWobble[i, 3] += systemList[j].Mass * systemList[j].MovementData[i, 1]

        # Divide all total moments by system total mass
        starWobble[i, 2:4] = starWobble[i, 2:4] / totalSystemMass

    # Get vx* and vy* from star's vx and vy
    starWobble[:, 4:6] = (totalPlanetMass) / (totalSystemMass) * systemList[0].MovementData[:, 2:4] 

    # Convert all starWobble values to km
    starWobble = (starWobble / 1000) * pm.au

    # Convert observer vx* and vy* into meters
    starWobble[:, 4:6] *= 1000

    # Return array for plotEverything()
    return starWobble
# End getStarWobble() -------------------------------------------------------------------




# Begin plotEverything() ----------------------------------------------------------------
def plotEverything(figName,
                   systemList = [0], plotOrbits = False,
                   systemMomentum = [0], plotMomentum = False,
                   systemEnergy = [0], plotEnergy = False, plotEnergyTot = False, plotEnergyLog = False,
                   viewProbs = [0], plotProbsP4 = False, plotProbsHist = False,
                   starWobble = [0], plotWobbleTop = False, plotWobbleCoM = False, plotWobbleRad = False, plotWobbleObs = False):
    # The All-Plotter. Takes your arrays and plots. Lots of plotting going on in here.
    # Technically able to plot everything at once. Technically don't recommend you try. Technically do one at a time.

    fig = plt.figure()

    # Figure out how many plots we have to deal with
    totalPlots = 0
    plotOrder = 1

    if plotOrbits == True:
        totalPlots += 1

    if plotMomentum == True:
        totalPlots += 1

    if plotEnergy == True:
        totalPlots += 1
        
    if plotEnergyTot == True:
        totalPlots += 1
        
    if plotEnergyLog == True:
        totalPlots += 1
        
    if plotProbsP4 == True:
        totalPlots += 1

    if plotProbsHist == True:
        totalPlots += 1

    if plotWobbleTop == True:
        totalPlots += 1

    if plotWobbleCoM == True:
        totalPlots += 1

    if plotWobbleRad == True:
        totalPlots += 1

    if plotWobbleObs == True:
        totalPlots += 2

    # Arrange the plots into an NxN pattern if possible
    totalRows = int(np.floor(np.sqrt(totalPlots)))
    totalColumns = totalRows
    print("Arranging plots in {0} row(s) and {1} column(s)..." .format(totalRows, totalColumns))

    # Start plotting

    if plotOrbits == True:
        axOrbits = fig.add_subplot(totalRows, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        for i in range(systemList[0].GBodyCount):
            axOrbits.plot(systemList[i].MovementData[:, 0], systemList[i].MovementData[:, 1], label=systemList[i].Name)

        axOrbits.axis((-80, 50, -50, 50))
        axOrbits.plot(0, 0, marker = 'o', color = 'r')
        axOrbits.set_xlabel('$1 x 10^{-3}$ AU')
        axOrbits.set_ylabel('$1 x 10^{-3}$ AU')
        axOrbits.grid()
        axOrbits.legend(loc = 'best', fontsize = 'small')
        axOrbits.suptitle = ('Calculated Orbits of System')

    if plotMomentum == True:
        axMomentum = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        axMomentum.plot(np.arange(0, timeMax - 1, 1), systemMomentum[:timeMax - 1, 1], label = 'Momentum Delta')

        axMomentum.grid()
        axMomentum.set_ylim([-1e-17, 1e-17])
        axMomentum.plot([timeMax, timeMax], axMomentum.get_ylim(), 'r--', lw = 2, label = 'TimeMax')
        axMomentum.set_xlabel('Time Steps for {} Days'.format(timeDaySim))
        axMomentum.set_ylabel('Change in Momentum')
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axMomentum.legend(loc = 'best', fontsize = 'small')
        axMomentum.set_title = ('Conservation of Momentum')

    if plotEnergy == True:
        # I dunno do something, I guess.
        axEnergy = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1

        axEnergy.plot(np.arange(0, timeMax - 1, 1), systemEnergy[:timeMax - 1, 1], label = 'Energy Delta')

        axEnergy.grid()
        axEnergy.plot([timeMax, timeMax], axEnergy.get_ylim(), 'r--', lw = 2, label = 'TimeMax')
        axEnergy.set_xlabel('Time Steps for {} Days'.format(timeDaySim))
        axEnergy.set_ylabel('Change in Energy')
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axEnergy.legend(loc = 'best', fontsize = 'small')
        axEnergy.set_title = ('Conservation of Energy')
        
    if plotEnergyTot == True:
        axEnergyTot = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1

        axEnergyTot.plot(np.arange(0, timeMax, 1), systemEnergy[:, 0], label = 'System Energy (Kinetic + Potential)')

        axEnergyTot.grid()
        axEnergyTot.plot([timeMax, timeMax], axEnergyTot.get_ylim(), 'r--', lw = 2, label = 'TimeMax')
        axEnergyTot.set_xlabel('Time Steps for {} Days'.format(timeDaySim))
        axEnergyTot.set_ylabel('Total Energy for {} System'.format(systemList[0].Name))
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axEnergyTot.legend(loc = 'best', fontsize = 'small')
        axEnergyTot.set_title = ('Total System Energy Over Time')
        
    if plotEnergyLog == True:
        axEnergyLog = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1

        axEnergyLog.plot(np.arange(0, timeMax, 1), systemEnergy[:, 2], label = "Relative Error")
        
        axEnergyLog.grid()
        axEnergyLog.plot([timeMax, timeMax], axEnergyLog.get_ylim(), 'r--', lw = 2, label = 'TimeMax')
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axEnergyLog.set_xlabel('Time Steps for {} Days'.format(timeDaySim))
        axEnergyLog.set_ylabel('Relative Error in Energy')
        axEnergyLog.legend(loc = 'best', fontsize = 'small')
        axEnergyLog.set_title = ('Relative Log 10 Change in Energy')
        
    if plotProbsP4 == True:
        axProbsP4 = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1
        
        axProbsP4.plot(np.arange(0, timeMax, 1), viewProbs[:, 3], 'g-', label = systemList[4].Name)        
        axProbsP4.plot([timeMax, timeMax], axProbsP4.get_ylim(), 'r--', lw = 2, label = 'TimeMax')

        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axProbsP4.grid()
        axProbsP4.set_xlabel('Time Steps for {} Days'.format(timeDaySim))
        axProbsP4.set_ylabel('Number of Planets Seen')
        axProbsP4.legend(loc = 'best', fontsize = 'small')
        axProbsP4.set_title = ('{} Planets Seen'.format(systemList[4].Name))

    if plotProbsHist == True:
        axProbsHist = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        axProbsHist.hist(viewProbs, bins = np.linspace(0, 5, 6) - 0.5, normed = 1, histtype='bar', alpha = 0.7, label = [systemList[1].Name,
                                                                                                                         systemList[2].Name,
                                                                                                                         systemList[3].Name,
                                                                                                                         systemList[4].Name,
                                                                                                                         systemList[5].Name,
                                                                                                                         systemList[6].Name])

        axProbsHist.axis([-0.5, 5.5, 0, 1])
        axProbsHist.set_xlabel('N Planets Seen')
        axProbsHist.set_ylabel('Probability of Occurring on Any Time Step in {0} Days' .format(timeDaySim))
        axProbsHist.legend()
        axProbsHist.set_title = ('Probabilty of Viewing N Planets')

    if plotWobbleTop == True:
        axWobbleTop = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        # Get star radius in km
        starRadius = (pm.star_radius_localunits / 1000) * pm.au

        axWobbleTop.plot(starWobble[:, 0], starWobble[:, 1], lw = 1, label = "{} Movement".format(systemList[0].Name))
        axWobbleTop.plot(852.5, 0, marker = 'o', color = 'r', label = "System Center of Mass")

        axWobbleTop.axis([-1600, 1600, -1600, 1600])
        axWobbleTop.grid()
        axWobbleTop.legend(loc = 'best', fontsize = 'small')
        axWobbleTop.set_xlabel('x (km)')
        axWobbleTop.set_ylabel('y (km)')
        axWobbleTop.set_title = ('Zoom-In of Top Down Orbits Diagram on {}'.format(systemList[0].Name))

    if plotWobbleCoM == True:
        axWobbleCoM = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        # Get star radius in km
        starRadius = (pm.star_radius_localunits / 1000) * pm.au

        axWobbleCoM.plot(starWobble[:, 2], starWobble[:, 3], 'r--', lw = 3, label = "System Center of Mass")

        axWobbleCoM.axis([852, 853, -5.95e-11, 5.95e-11])
        axWobbleCoM.grid()
        axWobbleCoM.legend(loc = 'best', fontsize = 'small')
        axWobbleCoM.set_xticks([852.43564])
        axWobbleCoM.set_xlabel('x (km)')
        axWobbleCoM.set_ylabel('y (km)')


    if plotWobbleRad == True:
        axWobbleRad = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        # Get star radius in km
        starRadius = (pm.star_radius_localunits / 1000) * pm.au
        starDraw = plt.Circle((0, 0), starRadius, color='r', alpha = 0.2, lw = 5)

        axWobbleRad.plot(starWobble[:, 0], starWobble[:, 1], lw = 3, label = "{} Movement".format(systemList[0].Name))
        axWobbleRad.plot(852.5, 0, marker = 'o', color = 'r', label = "System Center of Mass")
        axWobbleRad.add_artist(starDraw)

        axWobbleRad.text(-0.6 * starRadius, starRadius * 1.2, "{}".format(systemList[0].Name), color = "red", fontsize = 24)

        axWobbleRad.axis([-starRadius * 2, starRadius * 2, -starRadius * 2, starRadius * 2])
        axWobbleRad.grid()
        axWobbleRad.legend(loc = 'best', fontsize = 'small')
        axWobbleRad.set_xticks([-starRadius * 2, -starRadius, 0, starRadius, starRadius * 2])
        axWobbleRad.set_yticks([-starRadius * 2, -starRadius, 0, starRadius, starRadius * 2])
        axWobbleRad.set_xlabel('x (km)')
        axWobbleRad.set_ylabel('y (km)')

    if plotWobbleObs == True:
        axWobbleObsBig = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        axWobbleObsBig.plot(np.arange(0, timeMax, 1), starWobble[:, 5], label = "{} y-velocity for observer".format(systemList[0].Name))

        axWobbleObsBig.grid()
        axWobbleObsBig.legend(loc = 'best', fontsize = 'small')
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axWobbleObsBig.set_xlabel('Time Steps For {} Days'.format(timeDaySim))
        axWobbleObsBig.set_ylabel('y (m/s)')
        axWobbleObsBig.set_title = ('Observer View of Star Y-Velocity')

        # Zoomed-in version of above subplot
        axWobbleObsSmall = fig.add_subplot(totalPlots, totalColumns, plotOrder)
        plotOrder += 1 # Spot taken, move on.

        axWobbleObsSmall.plot(np.arange(0, timeMax, 1), starWobble[:, 5], label = "{} y-velocity for observer".format(systemList[0].Name))
        axWobbleObsSmall.plot(np.arange(0, timeMax, 1), -np.ones(timeMax), 'r--', label = "Doppler Spectrometer Detection Limit")
        axWobbleObsSmall.plot(np.arange(0, timeMax, 1), np.ones(timeMax), 'r--')

        axWobbleObsSmall.axis([0, timeMax, -20, 20])
        axWobbleObsSmall.grid()
        axWobbleObsSmall.legend(loc = 'best', fontsize = 'small')
        plt.xticks([0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)], [0, int(timeMax/3), int(2 * timeMax/3), int(timeMax)])
        axWobbleObsSmall.set_xlabel('Time Steps For {} Days'.format(timeDaySim))
        axWobbleObsSmall.set_ylabel('y (m/s)')
        axWobbleObsSmall.set_title = ('Zoom-in of Observer View of Star Y-Velocity')

    if totalPlots != 0:
        plt.tight_layout() # My Best Friend
        plt.savefig(figName, bbox_inches='tight')
        plt.show()

    # Returns nothing
    return
# End plotEverything() ------------------------------------------------------------------
