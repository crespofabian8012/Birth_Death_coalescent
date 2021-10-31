How to cite:
Fausto F. Crespo, David Posada, Carsten Wiuf,
Coalescent models derived from birth–death processes,
Theoretical Population Biology,
Volume 142,
2021,
Pages 1-11,
ISSN 0040-5809,
https://doi.org/10.1016/j.tpb.2021.09.003.
(https://www.sciencedirect.com/science/article/pii/S0040580921000654)
Abstract: A coalescent model of a sample of size n is derived from a birth–death process that originates at a random time in the past from a single founder individual. Over time, the descendants of the founder evolve into a population of large (infinite) size from which a sample of size n is taken. The parameters and time of the birth–death process are scaled in N0, the size of the present-day population, while letting N0→∞, similarly to how the standard Kingman coalescent process arises from the Wright–Fisher model. The model is named the Limit Birth–Death (LBD) coalescent model. Simulations from the LBD coalescent model with sample size n are computationally slow compared to standard coalescent models. Therefore, we suggest different approximations to the LBD coalescent model assuming the population size is a deterministic function of time rather than a stochastic process. Furthermore, we introduce a hybrid LBD coalescent model, that combines the exactness of the LBD coalescent model model with the speed of the approximations.
Keywords: Bernoulli sampling; Conditioned reconstructed process; Founder population; Variable population size coalescence; Kingman coalescence


The function to simulate from Limit Birt-Death Coalescent and approximations(M_K. M* and Hybrid)
can be found in the file:
/src/CoalSimulationBirthDeath.R

The simulation code  can be executed  in parallel or sequentially.
Every method has a parallel version(ended in  .parallel)

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

![shiny](https://github.com/crespofabian8012/Birth_Death_coalescent/blob/master/src/LimitBirthDeathCoalescent/shiny_app.png?raw=true)



