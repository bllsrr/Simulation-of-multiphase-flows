# Simulation of multiphase flows

These codes are the result of a month long internship at Institut d'Alembert, under the supervision of Daniel Fuster. The goal was to simulate multiphase flows and investigate the thermal and dynamic interactions between a droplet and a surrounding fluid. Modeling multiphase flows is central in modern engineering, as it is useful to solve many problems: cooling systems, nuclear reactor safety, etc. Accurately modeling heat transfer across such interfaces remains challenging, motivating comparisons between different simulation methods.  

### Approach

After running two simple experiments — putting a bubble of air and oil in water — we modified the values of the Reynolds and Péclet numbers to study their influence on the convergence of heat equation solutions. Three methods to reconstruct temperature were used, therefore each file has three variations: one where temperature is reconstructed from enthalpy (method 1), one where we consider the system as one fluid but with spatially-varying temperature (method 2) and one where we consider the system as two fluids, each with their own temperature (method 3). Finally, we used a simplified calculation method to run much quicker simulations.  

### Results

Our work showed that:  
- Advection plays a critical role in the thermal homogenization process, and significantly dictates the time required to reach thermal equilibrium.   
- Methods 2 and 3 are rigorously equivalent.  
- Methods 2 and 3 seem to produce overall fewer errors than method 1.  

### Prerequisites

- The CFD software *Basilisk* to run the code.  
- The module *SciencePlots* to plot results (can be manually removed in the code). 

### Repository structure

- [`air_in_water/`](air_in_water/) and [`oil_in_water/`](oil_in_water/) are the first two experiments made to get to grips with the software *Basilisk*.   
- [`peclet_influence/`](peclet_influence/) features the code used to study the Péclet number influence on the simulations.  
- [`reynolds_influence/`](reynolds_influence/) applies the same analysis to the Reynolds number.  
- [`simplified_calculations/`](simplified_calculations/) calculates a constant smoothed volume fraction to make calculations much quicker.  
- [`report.pdf`](report.pdf) provides more details about our simulations and our results.

### Disclaimers

- The provided code only features the code necessary to produce data, as well as the code necessary to produce comparative plots.  
- For the code producing comparative plots, the directories need to be changed.  

### Credits

- Daniel Fuster for providing the base code, upon which every code relies.  
- John D. Garrett for the module *SciencePlots*, used to produce high quality plots.
