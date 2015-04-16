Host diversity increases symbiont diversity and reduces transmission
=======================================================================

*Maxwell B. Joseph ([maxwell.b.joseph@colorado.edu](mailto:maxwell.b.joseph@colorado.edu)), Dept. of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA, 80309*

*Joseph R. Mihaljevic ([joseph.mihaljevic@colorado.edu](mailto:joseph.mihaljevic@colorado.edu)), Dept. of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA, 80309*

*Pieter T. J. Johnson ([pieter.johnson@colorado.edu](mailto:pieter.johnson@colorado.edu)), Dept. of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA, 80309*

Abstract
--------


Introduction
------------

Debate around effect of host diversity on disease. 
Increased availability of microbial data. 
Need for theoretical framework to make sense of data and predict future observations.

General strategy: 
- formalize the notion that host infection dynamics potentially have similarities to free living metacommunity dynamics
- develop a model that can capture the complexity of multi-host, multi-symbiont systems
- advantages (better realism, adequate complexty)
- disadvantages (less tractable than a simple analytical mathematical model)

Static host communities
--------------------------

### Model structure

We first develop a model for symbiont transmission that assumes static host communities that are neither reproducing nor dying. 
The local community is made up of multiple species, each of which has some abundance within the community. 
Host species vary in traits relevant to symbiont infection probability. 
We assume for simplicity that this variation is unidimensional, as would be the case if there were one primary niche axis for the symbionts (e.g. pH of the host gut). 
Symbiont niches are represented as gaussian functions with niche optima and breadths (means and variances, respectively) along this axis, such that the probability of successful infection for an individual host is a function of the host's value along the niche axis (Figure 1). 
Each symbiont has the same total probability of colonizing in the regional pool as a result of adjusting the niche functions of symbiont species with niche optima near the boundaries of  the niche space of the regional pool (e.g. symbiont 2 in Figure 1).

![Symbiont niches are represented by gaussian functions, such that each symbiont species has a niche optimum and a niche breadth. For any given host community, each host species has a single value of within-host condition represented by vertical tick marks above the x-axis. The range of host conditions represented in a community is the functional diversity of that community. This will always be contained within the total niche space of the regional symbiont pool.](/home/max/Documents/manuscripts/abm/paper/img/niche.png)

Every individual host has some constant probability of a symbiont from the regional pool attempting to infect, with an infection success probability being drawn from the symbiont niche function, given the host species condition. 
If a symbiont successfully colonizes a host, it becomes part of the local community and can either be extirpated through host recovery or persist by transmitting to other hosts (Figure 2). 
Every infected host has some recovery rate that is assumed constant across all species. 
Transmission occurs according to either density-dependent or frequency-dependent dynamics, where contact rates between infected and susceptible hosts are represented by $\beta S I$ or $\dfrac{\beta S I}{N}$, respectively, where $S$ is the number of susceptible hosts, $I$ is the number of infected hosts, $N$ is the total number of hosts, and $\beta$ is the per capita contact rate. 
For the case with static host communities, $N$ is constant, so that the dynamics of density and frequency dependent transmission are not appreciably different, other than a potential reduction in contact rates with frequency-dependent transmission (though this decrease can be offset by increasing $\beta$, leading to identical contact rates). 
Note that these terms are representing contacts between hosts that could result in transmission, but that the success of transmission still depends on the relationship between the condition of the susceptible host relative to the symbiont niche. 
Hosts are infected by symbionts at a particular location that can only be occupied by one type of symbiont (no co-infection). 

![The dynamics of the static host model include symbiont colonization from the regional pool, transmission and recovery of hosts within the local community, and the potential for extirpation of symbionts when all infected hosts are recover.](/home/max/Documents/manuscripts/abm/paper/img/overall.pdf)

We implement the previous model via stochastic simulation, generating continuous-time Markov chains in high-dimensional space, where the dimensions represent the states (number of hosts of each species, and the symbionts that infect each individual). 
Although continuous-time Markov processes are more difficult computationally than discrete-time formulations, they have the advantage of not requiring an arbitrarily specified order of events at each time step. 
Instead, the order of events emerges from the rates of each potential process via the Gillespie algorithm [@Gillespie1976]. 
Specifically, given the state of the system at any time point, every potential event that could happen is represented by some rate $r_e$ for events $e = 1, ..., E$. 
The total rate of events is simply the sum of the component rates $r_{tot} = \sum_{1}^{E} e_r$, and the time until an even happens is exponentially distributed, with overall rate parameter $r_{tot}$. 
This provides a way to generate the timepoints stochastically, and given some time of an event happening, the specific event is be selected with probability $\dfrac{e_r}{r_{tot}}$. 
This algorithm results in an exact solution for the stochastic process of interest. 

### Results

In many ways the structure and results of this static host model mirrors the structure of similar models for free-living species occupying local communities that vary in habitat heterogeneity [@Allouche2012]. 
The aforementioned assumptions lead to a unimodal relationship between symbiont diversity and host functional diversity emerging as a result of the dual impacts of host functional diversity on the potential for symbiont colonization and transmission. 
In low diversity communities, all hosts are very functionally similar, limiting the number of symbionts from the regional pool that can colonize, but facilitating tranmission following colonization. 
Here, host contacts are very likely to result in transmission events, because similar hosts are more able to share symbionts. 
As diversity of the host community increases, more symbiont species from the regional pool are able to colonize.
Because total community density is fixed (we are only changing host functional diversity), the number of hosts of any one condition in the local host community necessarily decreases as diversity increases, analogous to the habitat area-heterogeneity trade-off in free-living species (Figure 3). 

![Results of 1000 simulations of the static host model. Host functional diversity was varied, keeping all other parameters constant. Results of individual simulations encoded as the average symbiont richness across the entire simulated time-series are shown as points, and the line illustrates the fit of a 2nd degree polynomial regression.](/home/max/Documents/manuscripts/abm/paper/img/static.pdf)

These results may be of limited applicability, because in practice host communities are dynamic, and infection often results in changes in fitness due to parasitic or mutualistic relationships. 
However, if infection processes occur on timescales that are extremely rapid relative to host communities, the assumption of static host communities may be fairly reasonable. 
Nonetheless, we next consider the case when host communities are dynamic, with the added benefit of being able to investigate the influence of parasitism and mutualism on host-symbiont diversity relationships. 

Dynamic host communities
-------------------------

### Model structure

Next, we consider a local habitat patch that contains populations of hosts that are infected by symbionts. 
The host populations are open, with deaths, births, emigration and immigration.
Symbionts can still be transmitted among host individuals, or colonize from outside of the local host population from the regional pool, and each infected host as a recovery rate as before. 
There are multiple host species in a regional pool that can colonize the local community, and they vary in traits that are relevant to symbiont infection, but the communities are otherwise neutral, in that all host species have equal colonization, reproduction, and death rates. 
Hosts occur in a homogenous landscape, and can colonize the local habitat from the regional pool, reproduce, and die. 
Offspring attempt to disperse to a random habitat patch, and if it is unnoccupied, they successfully colonize.

Symbionts can be parasitic, commensal, or mutualistic through their effects on host survival.
To ensure that we can investigate a range of conditions, we modify the host survival rate $r$ as a function of infection:

$$r = d(1 + \beta_d * I(infected))$$

where $d$ is the background death rate, and $I(infected)$ is an indicator function that returns $1$ if an individual host is infected. 
The parameter $\beta_d$ controls the effect of infection on survival. 
When $\beta_d = 0$, the symbiont has no effect on survival, positive values impose increases in the death rate, and values in the interval $(-1, 0]$ represent decreases in death rate. 
We would expect that by being virulent and killing hosts, persistence in a local host community becomes more difficult because symbionts are destroying their resource, wheras mutualists can persist more easily because by increasing the longevity of their hosts, infected hosts can produce more offspring over the course of their lifetime.
As a consequence, for a given level of host diversity, parasites may have a lower value of expected value of richness than commensals, which have a lower expected richness compared to mutualists. 

### Results

Host diversity increased symbiont diversity, but the effect was non-linear (Figure 4). 
Furthermore, although we did observe a similar 2nd degree polynomial type of relationship, we did not observe a unimodal relationship between host and symbiont diversity across any range of parameter values. 

![Relationship between host functional diversity and symbiont diversity for commensal symbionts that have no effect on host survival. Each point represents the average symbiont over time richness in the local community (y-axis position), and the range of within host conditions available in the host regional pool (x-axis position).](/home/max/Documents/manuscripts/abm/paper/img/commensals.pdf)

As expected, when symbionts reduced host longevity, persistence within the host community becomes less probable for each symbiont species, resulting in a reduction in mean symbiont richness. 
Symbionts that increased host longevity showed the reverse pattern: persistence becomes more likely when symbionts are mutualists, resulting in higher mean symbiont richness (Figure 5). 

![Contour plot showing the linear predictor from a second degree polynomial linear regression of mean symbiont richness on host functional diversity. Contour lines represent isoclines for mean symbiont richness, and colors indicate the value of the linear predictor.](/home/max/Documents/manuscripts/abm/paper/img/contour.pdf)

Results for specialists vs. generalists...

Discussion
--------------------------------------------------------

Static host communities 
- extension of free living species concept, but not so applicable
- observed unimodal relationship

Dynamic host communities
- more realistic
- difference between potential vs. realized host richness
- not as strong unimodal relationships observed, because it's very unlikely to have tons of very small, very different host populations (small populations face a higher risk of stochastic extinctions)


See Dunn et al. 2010 for non-linear relationship

See Rotstock et al. 2014 for similar resolution of paradox

Conclusion
------------------

Acknowledgments
---------------

References
----------
