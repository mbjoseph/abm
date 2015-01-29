# host_init() initializes the host species traits w.r.t. environmental conditions

# Arguments:
## cells = number of cells in environment (max number of hosts)
## Emin = environmental condition minimum
## Emax = environmental condition maximum
## H = number of hosts
## sig.h = host niche breadth

# Returns:
## hniche.d = data frame of host niche data
## P.col = cell X host species array of colonization probabilities

host_init <- function(cells, H, pcol){
  host.species <- rep(1:H, each=cells)
  environmental.condition <- rep(0, H)
  Pcol <- array(pcol, dim=c(cells, H))  
 return(list(Pcol = Pcol))
}