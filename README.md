# project-2-2017-bogus-students

Last Updated 3/30/2017 by Brian Pickens

To Run:

Run 'Project2/Submission/project02.py' in the terminal, it will produce a graph showing the TRAPPIST-1 system at a default 24 days. This includes the movement of the star.
Inside project02.py are some plotting preference booleans which can be switched on or off for different plots to be shown. The numbers of days calculated can be modified in project02.py as well. The time values in simulation.py are default for quick testing and are overwritten by the project02.py values.
Additionally, the number of days simulated can be edited in project02.py. Run 'py.test' to test orbital periods and system momentum consistency.

Design:

|parameters.py| contains essential measured parameters for determining intial conditions. Used by all other files. Obtained from Prof. Oliver Beckstein

|simulation.py| contains the meat of the project, entirely consisting of class and function definitions.

|test_orbits.py| contains unit tests for determining the consistency of the orbits created by simulation.py. Runs with 'py.test' command. Tests orbital periods and system momentum consistency.

|project02.py| actually runs the simulation, using simulation.py as a module

The system TRAPPIST-1 is kept as a list of objects defined by the 'GBody' class. The class has a function for calculating one step of the Velocity Verlet algorithm (with its own function for calculating total gravity on itself). Each object simulates a step one at a time, through the total time steps in a loop found in simulateOrbits().

From this we obtain a list of objects containing their own orbital data. Several functions exist in simulations.py for interpreting, analyzing, and plotting this data.

An algorithm inefficiency is present in that redundant gravity calculations for velocity verlet are present in the current design. This could be resolved by storing the calculated values in a global GBody() array for use each time step, but there was not enough time implement this.
