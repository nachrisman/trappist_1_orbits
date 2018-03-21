import numpy as np
import pytest
import parameters as pm
from numpy.testing import assert_array_almost_equal

def test_import_module():
    try:
        import simulation
    except ImportError:
        raise AssertionError("Cannot import simulation.py module")
    try:
        import parameters
    except ImportError:
        raise AssertionError("Cannot import parameters.py module")
        
def test_import_class():
    try:
        from simulation import GBody
    except ImportError:
        raise AssertionError("Cannot import the simulation.GBody() class")

def test_import_function():
    try:
        from simulation import simulateOrbits
    except ImportError:
        raise AssertionError("Cannot import the simulation.simulateOrbits() function")
     
    try:
        from simulation import createSystem
    except ImportError:
        raise AssertionError("Cannot import the simulation.createSystem() function")   
        
def test_import_variable():
    try:
        from simulation import timeStep
    except ImportError:
        raise AssertionError("Cannot import the timeStep variable")
    try:
        from simulation import timeMax
    except ImportError:
        raise AssertionError("Cannot import the timeMax variable")
        
def test_timeMax_minimum():
    from simulation import timeMax
    
    if timeMax < 1300:
        raise AssertionError("Variable 'timeMax' too small to test orbital periods. (Try 14+ days)")
        
def runOrbitalSimulation():
    from simulation import (GBody, simulateOrbits, createSystem)
    systemTRAPPIST_1 = createSystem("TRAPPIST-1")   
    simulateOrbits(systemTRAPPIST_1)
    return systemTRAPPIST_1

def test_orbital_accuracy():
    from simulation import (timeStep, timeMax, getSystemMomentum)
    system = runOrbitalSimulation()
    
    # Test Orbital Periods
    periodArray = np.array([[12.00991,  0.], [16.47243,  0.], 
                           [22.9408,  0.], [30.56445,  0.], 
                           [39.4373,  0.], [47.8511,  0.]])
    for i in range(system[0].GBodyCount - 1):
        assert_array_almost_equal(system[i + 1].MovementData[int(pm.planets['period'][i] / timeStep), :2], periodArray[i], decimal = 0, 
                                  err_msg="Period for {} is incorrect".format(system[i + 1].Name)) 
    
    # Test Orbital Momentum
    momentum = getSystemMomentum(system)    
    for i in range(timeMax - 1):
        if momentum[i, 0] > 1e-16:
            raise AssertionError("System momentum not zero at timestep {}".format(i))
    
