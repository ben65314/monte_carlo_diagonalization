# Parameters documentaion
In a parameter file, defining a value will always be inside a set of brackets \{\}.

## Always required parameters 
### number_of_sites :
The number of sites ($$N_c$$) of the system.  
Example : `number_of_sites: {2}`
### U:
The Coulomb potential term.
Example : `U: {2}`
### mu:
The chemical potential.
Example : `mu: {2}`
### sites_location:
Called for every site added.
Example : `site_location: {(0,0,0)}`,`sites_location: {(1,0,0)}`
### allowed_jump:
Called for every length of jump. We need to specify the jump energy. Each jump corresponds to a difference of positon between two sites_location.
Example : `allow_jump: {(0,0,0)} {-1}`

## Required parameters for `mcd_solver`
### sampling_space_size:
The number of states in the sampled space ($$N_f$$).
Example : `sampling_space_size: {50}`
### N:
Number of electrons in the bloc computed.
Example : `N: {2}`
### Sz:
Z component of the spin in the bloc computed.
Example : `Sz: {0}`
### beta:
Beta value for the Boltzmann weight.
Example : `beta: {0.2}`
### reticle:
Reduces the reach of each step of the sample. A reticle of 1 will sample like a depth-first-search, an infinite reticle would sample like a breadth-first-search. If the reticle is set to 0, the package will use it like an infinite reticle. 
Example : `reticle: {0}`
### truncated_cutoff:
Weight kept when truncating the computed ground state ($$w_t$$).
Example : `truncated_cutoff: {1}`
### compute_green:
Do you want to compute the Green Q-matrices (YES/NO).
Example : `compute_green: {YES}`
### added_spin:
Spin to use in the Green $$c_{i,\uparrow/\downarrow}$$. $$0$$ is down and $$1$$ is up.
Example : `added_spin: {1}`
### eta
Eta $$(\eta)$$ value in the Green function. 
Example : `eta: {0.05}`


## Optional
### initial_state:
If you want to specified the starting states of the sampling. If not called will use an antiferromagnetic state according to the $$N_e$$ and $$S_z$$.  
Example : `initial_states: {5,6,9}`

