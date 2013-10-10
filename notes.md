abm notes
========================================================

To do 
-----
1. Break up the original huge function into sensible, discrete subfunctions     
**Current order of operations**
- Define environmental parameters (E values)
- Define host community traits (w/in host E vals, niche parameters)
- Define parasite community traits (niche parameters)
- Initialize output objects
- For loop begins
- Parasites transmit within community
- Resolve parasite colonization conflicts
- Parasites invade from regional pool
- Hosts recover, parasites die
- Host birth/death process
- Host offspring colonize empty sites
- Resolve conflicts
- Host immigrants colonize empty sites


**Sensible subfunctions**
1) com_init(): define environment, host com traits, parasite com traits  
2) transmit(): parasites transmit  
3) resolve(): resolve simultaneous colonization conflicts  
4) rain_parasites(): parasite colonization from regional pool
5) rain_hosts(): host colonization from regional pool
6) host_lifecycle(): birth/death/recovery process for hosts
7) hosts_colonize(): newborn and immigrating hosts colonize empty sites

2. Remove priority effects for breeders over new immigrating colonists


### Psuedocode for new master function
ibm <- function(args){

	require(packages)
	com_init(pass relevant args)
	
	# initialize output
	...
	
	# begin for-loop
	for (t in 2:timesteps){
		host_lifecycle()
		# within host_lifecycle(), call rain_hosts() then resolve()
		
		transmit()
		# within transmit, call rain_parasites() then resolve()
		
		# calculate relevant metrics 
		richness
		pr.occ
		parasite.richness
		transmission rates
		density, etc.
		}
	# stop if not statements
	# call to return()

}
