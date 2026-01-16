#include "StatesArr.h"
#include "Structures.h"
#include "basicFunctions.h"
#include "electronManipulationFunctions.h"

#ifndef __StatesR_T_h__
#define __StatesR_T_h__
template<class StateType, class VectorType> class StatesR_T :
public StatesArr<Node,StateType,VectorType>//<Node>
{
private:
	//Allocate MEM
	void allocate_more_nodes(uint64_t MEM_allocated=0) {
        /**************************************************
        * Creates a bigger array for the tree
        *
        * Parameters
        * ----------
        * elec  : (Electrons) Number of up and down electrons
        * sites : (unsigned int) Number of sites of the system
        * size  : (unsigned long) Ignore elec and sites and sets a specific size
        *
        * Returns
        * -------
        * NONE
        ************************************************/
		//Doubles the MEM allocated if <memAllocated> is 0
		uint64_t new_MEM_allocated = (MEM_allocated == 0)
                                    ? this->arr.capacity()*2 : MEM_allocated;

		//If the new MEM allocated is smaller than the current one doesn't change anything
		if (new_MEM_allocated < this->arr.capacity()) return;


		std::vector<StateType> augmented_nodes;
		this->get_num_array(&augmented_nodes);
		this->arr.clear();
		this->arr.reserve(new_MEM_allocated);
		for (unsigned int i = 0; i < augmented_nodes.size(); i++) {
			add(augmented_nodes.at(i));
		}
	}

	void FermiDiracDistributionSampling() {
        /**************************************************
        * Samples states with a breadth-first-search monte carlo method
        *
        * Parameters
        * ----------
        * NONE
        *
        * Returns
        * -------
        * NONE
        ************************************************/
		//Reset the StatesArr
		this->remove_all();

		//Number of states per n
		std::vector<uLong> n_states_for_nU;

		//Get the max number of double occupation possible
		uShort max_nU = (this->electrons.up > this->electrons.down)
                                ? this->electrons.down : this->electrons.up;

		float fermi_mu = max_nU;

		uLong remaining_states = this->sys_sP.sampling_size;

		//Computes the number of states for each double occupation level
		for (uShort i = 0; i <= max_nU; i++) {
			uLong combinations = comb_specified(i, this->sys_hubP.n_sites,
                                    this->electrons.up, this->electrons.down);
			n_states_for_nU.push_back(combinations);

			//Searches for the level that cannot be fully filled according
            //to the samplingSize
			remaining_states -= combinations;
			//Here we use > because remaindingStates is uLong and
            //once it is < 0 it will become very big
			if (remaining_states > this->sys_sP.sampling_size) {
				fermi_mu = i-1;
				break;
			}
		}

		//States stored for the Monte-Carlo sampling
		decltype(this) current_states = new StatesR_T(50);

		//Creating layers of nU for all the full layers who won't
        //need Monte-Carlo sampling
		for (uShort i = 0; i <= fermi_mu; i++) {
			std::vector<StateType> all_combinations_for_this_nU;

			//Calculates all states for a given nUp, nDown, nU
			combination_all(this->electrons.up, this->electrons.down,
                   this->sys_hubP.n_sites, i, &all_combinations_for_this_nU);
			for (uLong j = 0; j < all_combinations_for_this_nU.size(); j++) {
				//Will take the last filled level has the starting point
                //in the Monte-Carlo sampling
				if (i == fermi_mu) {
                    current_states->add(all_combinations_for_this_nU.at(j));
                }

				//Adds the state to the StatesArr
				add(all_combinations_for_this_nU.at(j));
			}
		}

		//Start state if the array is empty
		if (fermi_mu < 0) {
			uLong added_state = create_anti_ferro(this->sys_hubP.n_sites,
                                    this->electrons.up, this->electrons.down);
			current_states->add(added_state);
			add(added_state);
		}

		//Complete sampling with Monte-Carlo
		int min_nU = fermi_mu;
		uLong MH_size = this->get_length();

		decltype(this) next_step_eval = new StatesR_T(50);

		while (MH_size < this->sys_sP.sampling_size) {
			for (uLong i = 0; i < current_states->get_length(); i++) {

				//We dont need to know how many double occupation we have
                //before the jump because we dont set the accepation factor
                //with the differnce of energy but only the energy of the
                //jumped state according to the Fermi-Dirac distribution.
                //But we need it if we use Boltzmann's distribution
				int current_nU = Hu(current_states->get_at(i),
                                    this->sys_hubP.n_sites);

				std::vector<StateType> new_states;
				Ht(current_states->get_at(i), &new_states, &this->sys_hubP);

				for (uLong j = 0; j < new_states.size(); j++) {
					int new_nU = Hu(new_states.at(j),this->sys_hubP.n_sites);
					if(new_nU < min_nU) continue;

					float randomGen = (float)rand() / (float)RAND_MAX;

					float dE = (new_nU - current_nU)*this->sys_hubP.u;
					if(randomGen < boltzmann_distribution_function(dE,
                                                    this->sys_sP.beta_MH)) {
						next_step_eval->add(new_states.at(j));
						if(!this->countains_element(new_states.at(j))){
							add(new_states.at(j));
							MH_size++;
						}
					}

					if (MH_size >= this->sys_sP.sampling_size) break;
				}
				if (MH_size >= this->sys_sP.sampling_size) break;
			}

			decltype(this) temp = current_states;
			current_states = next_step_eval;
			next_step_eval = temp;
		}

		delete current_states;
		delete next_step_eval;
	}

	void MHSamplingOfStates(unsigned long sampling_size, float beta,
                            unsigned long reticle) {
		/****************************************************************
		* Samples states arround the one given and grows a subspace
        * around of desired size
		*
		* Parameters
		* ----------
		* sampling_size : (unsigned long) size of the subspace desired
		* beta : (float) beta value desired for MonteCarlo sampling
		* reticle : (unsigned long) reticle of breadth-first search sampling
        *
		* Returns
		* -------
        * None
		*****************************************************************/
        std::cout<<"MH_SAMPLING"<<std::endl;
		//Initial states
		allocate_more_nodes(sampling_size);

        //Each nU counter
        int nu_poss_value = (this->electrons.up > this->electrons.down) ? \
            this->electrons.up : this->electrons.down;
        nu_poss_value += 1;
        std::vector<sType> nu_state_counter(nu_poss_value, 0);
        //Maximum number of possible state for each nu value
        std::vector<sType> nb_state_per_nu;
        for (int i = 0; i < nu_poss_value; i++) {
            sType n_states = comb_specified(i, this->sys_hubP.n_sites,
                                            this->electrons.up,this->electrons.down);
            nb_state_per_nu.push_back(n_states);
        }

        int filled_nu_layer = -1;

        //StatesArr* currentState = this->clone();
        decltype(this) current_state = new StatesR_T(50);
        //*currentState = StatesArr(arrR.size());

		for (unsigned int i = 0; i < this->arr.size(); i++) {
            sType starting_states = this->arr.at(i).key;
			current_state->add(starting_states);
            int nu_add = Hu(starting_states,this->sys_hubP.n_sites);
            nu_state_counter[nu_add] += 1;
		}

		bool tree_like_sampling = false;
		if(verbose == 13){tree_like_sampling = true;}

		unsigned long MH_size = this->arr.size();

		StateType size_current_step_eval = this->arr.size();
		StateType size_next_step_eval = reticle * size_current_step_eval;
		//StatesArr* nextStepEval = this->clone();
		decltype(this) next_step_eval = new StatesR_T(50);

        //Only used when verbose needs it to be
        std::vector<sType> proposed_array (nu_poss_value, 0);
        std::vector<sType> current_nu_state (nu_poss_value, 0);
        decltype(this) possible_state = new StatesR_T(50);
        //

		unsigned int g = 0;
		auto step1 = std::chrono::high_resolution_clock::now();
		auto step2 = std::chrono::high_resolution_clock::now();
		while (MH_size < sampling_size) {

			size_current_step_eval = current_state->get_length();


			step2 = std::chrono::high_resolution_clock::now();
            if (verbose > 4) {
                for (unsigned long i = 0; i < size_current_step_eval;i++)
                {
                    int this_nu = Hu(current_state->get_at(i),this->sys_hubP.n_sites);
                    current_nu_state[this_nu]++;
                }
                for (unsigned long i = 0; i < possible_state->get_length();i++)
                {
                    int this_nu = Hu(possible_state->get_at(i),this->sys_hubP.n_sites);
                    proposed_array[this_nu]++;
                }

                printf("Sampled(%5ld/%5ld) [",MH_size,sampling_size);
                for (int i = 0; i < nu_poss_value; i++) {
                    printf("%6ld/%ld", nu_state_counter[i],
                           nb_state_per_nu[i]);
                }
                printf("] \033[42m|\033[0m size_ceval:%6ld [", size_current_step_eval);
                for (int i = 0; i < nu_poss_value; i++) {
                    printf("%6ld ",current_nu_state[i]);
                }
                printf("] \033[42m|\033[0m size_prop:%8ld [", possible_state->get_length());
                for (int i = 0; i < nu_poss_value; i++) {
                    printf("%6ld ",proposed_array[i]);
                }
                printf("] \033[42m|\033[0m dt = %s\n",time_formating(step1,step2).c_str());

                proposed_array = std::vector<sType>(nu_poss_value, 0);
                possible_state->remove_all();
                current_nu_state = std::vector<sType>(nu_poss_value, 0);
            }
			step1 = std::chrono::high_resolution_clock::now();


            bool big_sample = reticle * size_current_step_eval > sampling_size;
			size_next_step_eval = big_sample ? sampling_size :
                reticle * size_current_step_eval;

			next_step_eval->remove_all();
			next_step_eval->allocate_more_nodes(size_next_step_eval);
			g++;

			float current_energy;
			std::vector<StateType> possible_new_state;

			for (StateType i = 0; i < size_current_step_eval; i++) {
				//Evolution of Hamiltonian of the current state
                //and energy of the current state
				current_energy = Hu(current_state->get_at(i),
                        this->sys_hubP.n_sites)*this->sys_hubP.u;
				possible_new_state.clear();
				Ht(current_state->get_at(i), &possible_new_state,
                        &this->sys_hubP);

                if (verbose > 4) {
                    //Compute proposed states for print
                    for (unsigned long i = 0; i < possible_new_state.size();i++)
                    {
                        possible_state->add(possible_new_state.at(i));
                    }
                }


				if (possible_new_state.size() != 0) {
					StateType new_state;
					std::vector<StateType> all_accepted_states;
					//Test the breadth algorithm for all possible states found
                    int new_nu;
					for (StateType j = 0; j < possible_new_state.size(); j++)
					{
						new_state = possible_new_state.at(j);
						if (tree_like_sampling) {
							if(this->countains_element(new_state)){
								continue;
							}
						}

						//#Calculates new Energy and accept factor
						new_nu = Hu(new_state,this->sys_hubP.n_sites);
                        float new_energy = new_nu * this->sys_hubP.u;

						float diff_energy = new_energy - current_energy;
						float a = (float)rand() / (float)RAND_MAX;
						bool accepted = exp(-beta * diff_energy) > a;

						//Energy MONTE CARLO Condition
						if (accepted){
                            if (new_nu >= filled_nu_layer) {
                                all_accepted_states.push_back(new_state);
                            }
						}
					}

					//If no states are accepted
					if (all_accepted_states.size() == 0) {
						continue;
					}
					//If the number of accepted states is bigger than reticle
					else if (all_accepted_states.size() > reticle) {
						int* rdm_arr = new int[all_accepted_states.size()];
						for (
                        StateType p = 0; p < all_accepted_states.size(); p++)
						{
							rdm_arr[p] = p;
						}
						std::shuffle(
                            rdm_arr, rdm_arr + all_accepted_states.size(),
                            std::default_random_engine(std::time(NULL)));

						for (StateType q = 0; q < reticle; q++) {
							StateType item = all_accepted_states.at(rdm_arr[q]);

							next_step_eval->add(item);
							if(!this->countains_element(item)){
								add(item);
								MH_size++;
								g=0;
							}
							if (MH_size >= sampling_size) {
								delete[] rdm_arr;
                                //This is a nested break to get out of the for
                                //loop of the evolution and the for loop of the
                                //currentState to evaluate.
								goto nestedBreakForEnoughSampling;
							}
						}
						delete[] rdm_arr;
					}
					//If the number of accepted states is less than reticle
					else {
						for (
                        StateType l = 0; l < all_accepted_states.size(); l++){
							StateType item = all_accepted_states.at(l);


							next_step_eval->add(item);
							if (!this->countains_element(item)) {
								add(item);
								MH_size++;
								g=0;

                                //Add to counter for each type of nu
                                new_nu = Hu(item,this->sys_hubP.n_sites);
                                nu_state_counter[new_nu] += 1;
                                // Prevents already filled sampling
                                if (nu_state_counter[new_nu] \
                                    == nb_state_per_nu[new_nu]) {
                                    filled_nu_layer++;
                                    if (verbose > 4) std::cout<<"FILLED : " <<filled_nu_layer<<std::endl;
                                }
                                //print_vector(nu_state_counter.data(),nu_state_counter.size());
							}
							if (MH_size >= sampling_size) {
                                //This is a nested break to get out of the for
                                //loop of the evolution and the for loop of the
                                //currentState to evaluate.
								goto nestedBreakForEnoughSampling;
							}
						}
					}
				}
				else {//Prevents solo block breaking
					std::cout << "No Evolution\n";
					std::cout << "Current: " << current_state->get_at(i)
                              << "   E : " << current_energy;
					break;
				}
                //std::cout<<"new total"<<this->get_length()<<std::endl;

			}//END OF FOR

            //std::cout<<"g:"<<g<<std::endl;
		nestedBreakForEnoughSampling:
			if (g > PERMISSION) {
				std::cout << "The sample size entered couldn't be met. "
                          << "This can be a result of:\n\t-An unattainable "
                          << "sample size\n\t-A beta value too large\n";
				break;
			}
			//std::swap(currentState,nextStepEval);
			if (next_step_eval->get_length()>0){
				decltype(current_state) temp = current_state;
				current_state = next_step_eval;
				next_step_eval = temp;
			}
			next_step_eval->remove_all();
		}

        delete possible_state;
		delete current_state;
		delete next_step_eval;
	}
	void MC_sampling(unsigned long sampling_size, float beta) {
		/****************************************************************
		* Samples states arround the one given and grows a subspace
        * around of desired size
		*
		* Parameters
		* ----------
		* sampling_size : (unsigned long) size of the subspace desired
		* beta : (float) beta value desired for Monte Carlo sampling
        *
		* Returns
		* -------
        * None
		*****************************************************************/
        std::cout<<"MC_SAMPLING"<<std::endl;
		//Initial states
		allocate_more_nodes(sampling_size);

		sType MH_size = this->arr.size();

        std::vector<sType> index_active;
        std::vector<int> neighbor;
        //Checks how many neighbors are already included in the sample.
        for (sType i = 0; i < MH_size; i++) {
            index_active.push_back(i);

            //Neighbor fixing
			std::vector<StateType> possible_new_state;
            Ht(this->get_at(i), &possible_new_state, &this->sys_hubP);

            int n_neighbor = possible_new_state.size();
            for (uint64_t j = 0; j < possible_new_state.size(); j++) {
                if (this->countains_element(possible_new_state.at(j))) {
                    n_neighbor--;
                }
            }
            neighbor.push_back(n_neighbor);
        }


		unsigned int g = 0;
		while (MH_size < sampling_size) {
			g++;
            //this->show_all_states();
            //print_vector(neighbor.data(), neighbor.size());

            //Generate a random state from the active ones
            sType rdm_element_1 = rdm_number() * (index_active.size()-1);
            sType rdm_index = index_active.at(rdm_element_1);

			float current_energy = Hu(this->get_at(rdm_index),
                            this->sys_hubP.n_sites) * this->sys_hubP.u;

            //Generate a random accessible jumping state from the selected one.
			std::vector<StateType> possible_new_state;
            Ht(this->get_at(rdm_index), &possible_new_state, &this->sys_hubP);
            sType rdm_element_2 = rdm_number() * (possible_new_state.size()-1);

            //Generate the energy of the new state and the difference.
            int new_nu = Hu(possible_new_state.at(rdm_element_2),
                            this->sys_hubP.n_sites);
            float new_energy = new_nu * this->sys_hubP.u;
            float diff_energy = new_energy - current_energy;

            //Generate a random number for the boltzmann condition.
            float rdm_element_3 = rdm_number();
            bool accepted = exp(-beta * diff_energy) > rdm_element_3;

            //Accept condition
            if (accepted){
                sType where_index;
                bool already_countained = this->where_is_element(
                    possible_new_state.at(rdm_element_2), &where_index);

                if(!already_countained){
                    this->add(possible_new_state.at(rdm_element_2));
                    index_active.push_back(MH_size);
                    MH_size++;
                    g=0;

                    // Neighbor reducing
                    std::vector<StateType> neighbor_possible;
                    Ht(possible_new_state.at(rdm_element_2), &neighbor_possible, &this->sys_hubP);
                    // For each alreay existing neighbor will reduce -1 at
                    // possible neighbor and reduce the total of already
                    // existing neighbor to this presente state.
                    int n_neighbor = neighbor_possible.size();
                    int n_already_neighbor = 0;
                    for (int j = 0; j < n_neighbor; j++) {
                        StateType index_neighbor;
                        if (this->where_is_element(neighbor_possible.at(j),
                                                   &index_neighbor)){
                            n_already_neighbor++; //Remove to the current state
                            neighbor[index_neighbor]--; //Remove the neighbor
                            //std::cout<<index_neighbor<<"/"<<neighbor.size()<<std::endl;
                            if (neighbor[index_neighbor] <= 0) {
                                auto remove_index = std::find(index_active.begin(),
                                          index_active.end(),index_neighbor);

                                if (remove_index != index_active.end()) {
                                    sType last_index = index_active.back();
                                    index_active.pop_back();
                                    *remove_index = last_index;
                                    if (verbose > 5) {
                                        printf("Sampled(%5ld/%5ld) | ",MH_size, sampling_size);
                                        printf("Blocking state %5ld | %6ld states possibly generated \n",
                                               this->get_at(index_neighbor),index_active.size());
                                    }

                                }
                            }
                        }
                    }
                    neighbor.push_back(n_neighbor-n_already_neighbor);
                }

            }
            if (verbose > 99) {
                printf("\nState picked (%8ld)  State proposed (%8ld)  Accepted? (%d)",
                       this->get_at(rdm_index),
                       possible_new_state.at(rdm_element_2),
                       accepted);
            }

			if (g > PERMISSION) {
				std::cout << "The sample size entered couldn't be met. "
                          << "This can be a result of:\n\t-An unattainable "
                          << "sample size\n\t-A beta value too large\n";
				break;
			}
		}
	}
public:
	//Constructors
	StatesR_T(uInt reserve = 10) {this->arr.reserve(reserve);}
	StatesR_T(const std::vector<StateType>* array_to_state, hubbardParam hubP){
		this->arr.reserve(array_to_state->size());
		for (uint32_t i = 0; i < array_to_state->size(); i++) {
			add(array_to_state->at(i));
		}

        this->set_hubbard_parameters(hubP);

		this->electrons = find_number_of_electron(this->arr[0].key,
                                                this->sys_hubP.n_sites);
	}
	//Destructors
	~StatesR_T(){}
	//Clone
	StatesR_T* clone(){
		StatesR_T* cloned_sArr = new StatesR_T(10);
		cloned_sArr->set_hubbard_parameters(this->sys_hubP);
		cloned_sArr->set_sampling_parameters(this->sys_sP);
		cloned_sArr->electrons.up = this->electrons.up;
		cloned_sArr->electrons.down = this->electrons.down;
		return cloned_sArr;
	}

	//Samplings
	void sampling(){
		MC_sampling(this->sys_sP.sampling_size,this->sys_sP.beta_MH);
		//MHSamplingOfStates(this->sys_sP.sampling_size,
        //                    this->sys_sP.beta_MH, this->sys_sP.reticle);
	}
	void sampling_least_energy(){
		FermiDiracDistributionSampling();
	}
	void subspace_condition_expanding() {
		if (verbose > 5) std::cout << "nHapply " << this->sys_sP.nHapply <<'\n';
		for(uInt i = 0; i < this->sys_sP.nHapply; i++) {
			Ht_subspace_condition_expanding(this,0,this->arr.size());
		}
	}

	StateType get_at(StateType index) const{//::
		/***************************************
		Gets a specific state in the StatesArr

		Parameters:
		-----------
		index : (int) element desired

		Returns:
		--------
		arr[index] : (States) state at the index 'index'
		****************************************/
		StateType el = 0;
		if (index >= 0 && index < this->get_length()) el = this->arr[index].key;
		else std::cout << "get_at searched out of the array(" << index << "/"
                       << this->get_length() << ")" << std::endl;
		return el;
	}

	//Function overload
	void matrix_creation(double* result_matrix) {
		/*******************************************************************
		Create the matrix of the block of state given

		Parameters
		----------
		bs : (StatesArr)Block of states to matrixify
		result_matrix: (double*) Receptacle to the matrix created
		Returns
		-------
		NONE
		*****************************************************************/
		double mu_value = compute_mu(this->sys_hubP.mu, this->electrons);
		uint64_t cols = this->get_length();

		for (StateType i = 0; i < cols; i++) {
			//Adds diagonal values
			//Hu and Hmu
			result_matrix[i * cols + i]
                += (double)Hu(this->get_at(i), this->sys_hubP.n_sites)
                * (this->sys_hubP.u);
			result_matrix[i * cols + i] += mu_value;

			std::vector<StateType> projHt;
			std::vector<double> jump_energy;
			t_jump_energy(this->get_at(i), &projHt, &jump_energy,
                 &this->sys_hubP);
			for (unsigned int j = 0; j < projHt.size(); j++) {
				StateType index;
				if (!this->where_is_element(projHt.at(j), &index)) continue;
				if (index < i) continue;
				result_matrix[i * cols + index] += jump_energy.at(j);
				result_matrix[index * cols + i] += jump_energy.at(j);
			}
		}
	}

	void H (VectorType* h_phi_n, VectorType* phi_n) {
		/***************************************************************
		* Applies H on the given vector without calculating
        * the H matrix for the r-basis
        *
		* Parameters
		* ----------
		* h_phi_n : (double *) receptacle for the projected vector
        *                       of Lanczos algorithm
		* phi_n   : (double *) current vector of Lanczos algorithm
		* states  : (StatesArr) States used in the subspace
        *
		* Returns
		* --------
		* NONE
		*****************************************************************/
		int elements =  this->get_length();

		//Mu value (const for every state)
		float mu_value = compute_mu(this->sys_hubP.mu, this->electrons);
		#pragma omp parallel for default(none) shared(elements, mu_value,\
            h_phi_n, phi_n, stdout)
		for(int i = 0; i < elements; i++){
			//Hu and hmu
			h_phi_n[i] += (double)(Hu(this->get_at(i),
                            this->sys_hubP.n_sites)*this->sys_hubP.u + mu_value)
                            * phi_n[i];

			//Ht
			std::vector<StateType> proj;
			std::vector<double> energies;
			t_jump_energy(this->get_at(i), &proj, &energies, &this->sys_hubP);

			for (unsigned int j = 0; j < proj.size(); j++) {
				StateType index;

				if (this->where_is_element(proj.at(j), &index)) {
					h_phi_n[i] += energies.at(j) * phi_n[index];
				}
			}
		}
	}


	void add(StateType el){
		/***************************************************************
        * Adds an element to the array
        *
		* Parameters
		* ----------
		* ek : (StateType) State to add
        *
		* Returns
		* --------
		* NONE
		*****************************************************************/
		if (this->arr.capacity() == this->arr.size()) allocate_more_nodes();
		this->arr.push_back(Node(el));
		bool alreadyThere = false;
		insert(this->arr.data(), this->arr.data() + this->arr.size() - 1,
                    &alreadyThere);

		if (alreadyThere && this->arr.size() > 1) this->arr.pop_back();
	}
	virtual bool where_is_element(StateType el, StateType* index) const {
		/***************************************
		* Searches the index of a given state if the array has it
        *
		* Parameters:
		* -----------
		* el : (StateType) state to look for
		* index : (StateType) where is the given state in the array
        *
		* Returns:
		* --------
		* found : (bool) has the element been found
		****************************************/
		bool found = true;
		//Its not there if there is nothing there.
		if (this->arr.size() == 0) {return false;}

		//If Node not present returns a NULL node
		const Node* finding = search(this->arr.data(),el);
		if (finding == NULL) found = false;
		else *index = finding - this->arr.data();

		return found;
	}

	std::string show_all_states_string() const{
		/***************************************
		Writes all the States countained in a string

		Parameters:
		-----------
		NONE

		Returns:
		--------
		allStates : (string) string of all the states
		****************************************/
		//Writing operator
		std::string t, n;
		t = "\t";
		n = "\n";

		std::string allStates = "[";
		for (StateType i = 0; i < this->arr.size(); i++){
			allStates += std::to_string(this->arr[i].key) + t;
		}
		allStates += "]\n";

		return allStates;
	}
};
#endif
