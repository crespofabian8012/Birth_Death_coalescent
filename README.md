# Birth_Death_coalescent
Coalescent approximations models to Birth Dearth processes

The simulations can run in parallel or sequentially.

The function to simulate from Limit Birt-Death Coalescent and approximations(M_K. M* and Hybrid)
can be found in the file:
/src/CoalSimulationBirthDeath.R

You can find examples of how to call the simulation functions in:
/src/simulator_BD_Coalescent.R

To can modify  the values of 
n: sample size at present time 
Delta: 
Gamma: dimensionless scaled growth rate 
in 
simulator_BD_Coalescent.R

Alternatively you can souce the file /src/simulator_BD_Coalescent.R
and call the functions on your own.

Theres is a Shiny app: app.R in 
https://github.com/crespofabian8012/Birth_Death_coalescent/tree/master/src/LimitBirthDeathCoalescent
where you can play with the paremeter values.



