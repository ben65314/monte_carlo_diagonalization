#include "StatesArr.h"

#ifndef __StatesR_T_h__
#define __StatesR_T_h__
template<class StateType, class VectorType> class StatesR_T : public StatesArr<Node,StateType,VectorType>//<Node>
{
private:
	//Allocate MEM
	void allocateMoreNodes(uint64_t memAllocated = 0) {
		//Doubles the MEM allocated if <memAllocated> is 0
		uint64_t newMEMallocated = (memAllocated == 0) ? this->arr.capacity()*2 : memAllocated;

		//If the new MEM allocated is smaller than the current one doesn't change anything
		if (newMEMallocated < this->arr.capacity()) return;
		
		std::vector<StateType> augmented_nodes;
		this->getNumArray(&augmented_nodes);
		this->arr.clear();
		this->arr.reserve(newMEMallocated);
		for (unsigned int i = 0; i < augmented_nodes.size(); i++) {
			add(augmented_nodes.at(i));
		}
	}
	void FermiDiracDistributionSampling() {
		//Reset the StatesArr
		this->removeAll();

		//Number of states per n
		std::vector<uLong> nStates_for_nU;

		//Get the max number of double occupation possible
		uShort max_nU = (this->electrons.up > this->electrons.down) ? this->electrons.down : this->electrons.up;

		float fermi_mu = max_nU;

		uLong remainingStates = this->sys_sP.samplingSize;
		
		//Computes the number of states for each double occupation level
		for (uShort i = 0; i <= max_nU; i++) {
			uLong combinations = combSpecified(i,this->sys_hubP.n_sites,this->electrons.up,this->electrons.down);
			nStates_for_nU.push_back(combinations);

			//Searches for the level that cannot be fully filled according to the samplingSize
			//std::cout<<"REMAINS:"<<remainingStates<<std::endl;
			remainingStates -= combinations;
			//Here we use > because remaindingStates is uLong and once it is < 0 it will become very big
			if (remainingStates > this->sys_sP.samplingSize /*|| remainingStates == 0*/) {
				fermi_mu = i-1;
				break;
			}
		}

		//std::cout<<"FMU:"<<fermi_mu<<std::endl;
		//States stored for the Monte-Carlo sampling
		decltype(this) currentStates = new StatesR_T(50);

		//Creating layers of nU for all the full layers who won't need Monte-Carlo sampling
		for (uShort i = 0; i <= fermi_mu; i++) {
			std::vector<StateType> all_combinations_for_this_nU;
			
			//Calculates all states for a given nUp, nDown, nU
			combinationAll(this->electrons.up,this->electrons.down,this->sys_hubP.n_sites,i, &all_combinations_for_this_nU);
			for (uLong j = 0; j < all_combinations_for_this_nU.size(); j++) {
				//Will take the last filled level has the starting point in the Monte-Carlo sampling
				if (i == (fermi_mu)) currentStates->add(all_combinations_for_this_nU.at(j));

				//Adds the state to the StatesArr
				add(all_combinations_for_this_nU.at(j));
			}
		}

		//Start state if the array is empty
		if (fermi_mu < 0) {
			uLong added_state = createAntiFerro(this->sys_hubP.n_sites,this->electrons.up,this->electrons.down);
			currentStates->add(added_state);
			add(added_state);
		}

		//Complete sampling with Monte-Carlo
		int min_nU = fermi_mu;
		uLong MH_size = this->getLength();

		decltype(this) nextStepEval = new StatesR_T(50);

		while (MH_size < this->sys_sP.samplingSize) {
			for (uLong i = 0; i < currentStates->getLength(); i++) {

				//We dont need to know how many double occupation we have before the jump because we dont set the accepation factor with the differnce of energy but only the energy of the jumped state according to the Fermi-Dirac distribution. But we need it if we use Boltzmann's distribution
				int current_nU = Hu(currentStates->getAt(i),this->sys_hubP.n_sites);

				std::vector<StateType> newStates;
				Ht(currentStates->getAt(i), &newStates, &this->sys_hubP);

				for (uLong j = 0; j < newStates.size(); j++) {
					int new_nU = Hu(newStates.at(j),this->sys_hubP.n_sites);
					if(new_nU < min_nU) continue;

					float randomGen = (float)rand() / (float)RAND_MAX;

					float dE = (new_nU - current_nU)*this->sys_hubP.u;
					if(randomGen < boltzmannDistributionFunction(dE, this->sys_sP.beta_MH)) {
					//if (randomGen <fermiDirac_acceptation_coeficients.at(new_nU)) {
						nextStepEval->add(newStates.at(j));
						if(!this->countainsElement(newStates.at(j))){
							add(newStates.at(j));
							MH_size++;
							//std::cout<<"MH_SIZE:"<<MH_size<<std::endl;
						}
					}
					
					if (MH_size >= this->sys_sP.samplingSize) break;
				}
				if (MH_size >= this->sys_sP.samplingSize) break;
			}

			decltype(this) temp = currentStates;
			currentStates = nextStepEval;
			nextStepEval = temp;
		}
		
		delete currentStates;
		delete nextStepEval;
	}
	std::tuple<unsigned long,unsigned long> MHSamplingOfStates(unsigned long samplingSize, float beta, unsigned long reticle){
		/****************************************************************
		Samples states arround the one given and grows a subspace around of desired size
		
		Parameters
		----------
		samplingSize : (unsigned long) size of the subspace desired
		beta : (float) beta value desired for MonteCarlo sampling
		reticle : (unsigned long) reticle of breadth-first search sampling

		Returns
		-------
		Accept / Evaluated
		*****************************************************************/
		//Initial states
		allocateMoreNodes(samplingSize);
		//StatesArr* currentState = this->clone();
		decltype(this) currentState = new StatesR_T(50);
		//*currentState = StatesArr(arrR.size());

		for (unsigned int i = 0; i < this->arr.size(); i++) {
			currentState->add(this->arr.at(i).key);
		}

		bool TreeLikeSampling = false;
		if(verbose == 10){TreeLikeSampling = true;}

		unsigned long MH_size = this->arr.size();

		unsigned long accept = 0, evaluated = 0;

		StateType size_currentStepEval = this->arr.size();
		StateType size_nextStepEval = reticle * size_currentStepEval;
		//StatesArr* nextStepEval = this->clone();
		decltype(this) nextStepEval = new StatesR_T(50);

		//*nextStepEval = StatesArr(size_nextStepEval);

		unsigned int g = 0;
		auto step1 = std::chrono::high_resolution_clock::now();
		auto step2 = std::chrono::high_resolution_clock::now();
		while (MH_size < samplingSize) {	
			step2 = std::chrono::high_resolution_clock::now();
			if(verbose > 5) std::cout<<"\nCURRENT SAMPLE SIZE:"<<MH_size<<" ("<<(float)MH_size / samplingSize<<std::endl;
			if(verbose > 5) std::cout<<"CurrentStateLen:"<<currentState->getLength()<<std::endl;
			if(verbose > 5) std::cout<<"Time to gather current sample:"<<timeFormating(step1,step2)<<std::endl;
			step1 = std::chrono::high_resolution_clock::now();
			size_currentStepEval = currentState->getLength();
			size_nextStepEval = (reticle * size_currentStepEval > samplingSize) ? samplingSize : reticle * size_currentStepEval;

			if(verbose > 5) std::cout<<"Next step eval size:"<<size_nextStepEval<<std::endl;
			nextStepEval->removeAll();
			nextStepEval->allocateMoreNodes(size_nextStepEval);
			g++;

			float currentEnergy;
			std::vector<StateType> possibleNewState;
			
			for (StateType i = 0; i < size_currentStepEval; i++) {
				//#Evolution of Hamiltonian of the current state and energy of the current state
				currentEnergy = Hu(currentState->getAt(i),this->sys_hubP.n_sites); 
				possibleNewState.clear();
				Ht(currentState->getAt(i), &possibleNewState, &this->sys_hubP);

				if (possibleNewState.size() != 0) {
					StateType newState;
					std::vector<StateType> allPossibleAcceptedStates;
					//Test the breadth algorithm for all possible states found
					for (StateType j = 0; j < possibleNewState.size(); j++)
					{
						newState = possibleNewState.at(j);
						if (TreeLikeSampling) {
							if(this->countainsElement(newState)){
								continue;
							}
						}

						//#Calculates new Energy and accept factor
						float newEnergy;
						newEnergy = Hu(newState,this->sys_hubP.n_sites) * this->sys_hubP.u;

						float diffEnergy = newEnergy - currentEnergy;
						float a = (float)rand() / (float)RAND_MAX;
						bool accepted = exp(-beta * diffEnergy) > a;
						evaluated++;

						//std::cout<<"beta:"<<beta<<"\tret:"<<reticle<<std::endl;
						//Energy MONTE CARLO Condition
						if (accepted) {
							accept++;
							allPossibleAcceptedStates.push_back(newState);
						}
					}

					//Si le nombre d'états voulant être accepté est 0
					if (allPossibleAcceptedStates.size() == 0) {
						continue;
					}
					//Si le nombre d'états voulant être accepté est plus grand que le réticule
					else if (allPossibleAcceptedStates.size() > reticle) {
						//std::cout<<"possAcc:"<<allPossibleAcceptedStates.size()<<"\tret:"<<reticle<<std::endl;
						int* rdmArray = new int[allPossibleAcceptedStates.size()];
						for (StateType p = 0; p < allPossibleAcceptedStates.size(); p++)
						{
							rdmArray[p] = p;
						}
						std::shuffle(rdmArray, rdmArray + allPossibleAcceptedStates.size(), std::default_random_engine(std::time(NULL)));

						for (StateType q = 0; q < reticle; q++) {
							StateType item = allPossibleAcceptedStates.at(rdmArray[q]);
							//ins(item,&nextStepEval);
							nextStepEval->add(item);
							if(!this->countainsElement(item)){
								add(item);
								MH_size++;
								g=0;
							}
							if (MH_size >= samplingSize) {
								delete[] rdmArray;
								goto nestedBreakForEnoughSampling;//This is a nested break to get out of the for loop of the evolution and the for loop of the currentState to evaluate.
							}
						}
						delete[] rdmArray;
					}
					//Si le nombre d'états voulant être accepté est inférieur au réticule
					else {
						for (StateType l = 0; l < allPossibleAcceptedStates.size(); l++){
							StateType item = allPossibleAcceptedStates.at(l);

							//binary_vecSet(item,&nextStepEval);
							nextStepEval->add(item);
							if(!this->countainsElement(item)){
								add(item);
								MH_size++;
								g=0;
							}
							if (MH_size >= samplingSize) {
								goto nestedBreakForEnoughSampling;//This is a nested break to get out of the for loop of the evolution and the for loop of the currentState to evaluate.
							}
						}
					}	
				}
				else {//Prevents solo block breaking
					std::cout << "No Evolution\n";
					std::cout << "Current: " << currentState->getAt(i) << "   E : " << currentEnergy;
					break;
				}     

			}//END OF FOR


		nestedBreakForEnoughSampling:
			if (g > PERMISSION) {
				std::cout << "The sample size entered couldn't be met. This can be a result of:\n\t-An unattainable sample size\n\t-A beta value too large\n";
				break;
			}
			//std::swap(currentState,nextStepEval);
			if (nextStepEval->getLength()>0){
				decltype(currentState) temp = currentState;
				currentState = nextStepEval;
				nextStepEval = temp;
			}
			nextStepEval->removeAll();
		}

		delete currentState;
		delete nextStepEval;
		return {accept, evaluated};
	}
public:
	//Constructors
	StatesR_T(uInt reserve = 10) {this->arr.reserve(reserve);}
	StatesR_T(const std::vector<StateType>* array_to_state){
		this->arr.reserve(array_to_state->size());
		for (uint32_t i = 0; i < array_to_state->size(); i++) {
			add(array_to_state->at(i));
		}
		this->electrons = findNumberOfElectron(this->arr[0].key,this->sys_hubP.n_sites);
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
	void sampling_MH(){
		MHSamplingOfStates(this->sys_sP.samplingSize, this->sys_sP.beta_MH, this->sys_sP.reticle);
	}
	void sampling_least_energy(){
		FermiDiracDistributionSampling();
	} 
	void subspace_condition_expanding() {
		if (verbose > 5) std::cout<<"nHapply "<<this->sys_sP.nHapply<<std::endl;
		for(uInt i = 0; i < this->sys_sP.nHapply; i++) {
			Ht_subspace_condition_expanding(this,0,this->arr.size());
		}
	}

	StateType getAt(StateType index) const{//::
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
		if (index >= 0 && index < this->getLength()) el = this->arr[index].key; 
		else std::cout<< "getAt searched out of the array("<<index<<"/"<<this->getLength()<<")"<<std::endl;
		return el;
	}

	//Function overload
	void matrixCreation(double* result_matrix) {
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
		double muValue = computeMu(this->sys_hubP.mu, this);
		uint64_t cols = this->getLength();

		for (StateType i = 0; i < cols; i++) {
			//Adds diagonal values
			//Hu and Hmu
			result_matrix[i * cols + i] += (double)Hu(this->getAt(i),this->sys_hubP.n_sites) * (this->sys_hubP.u);
			result_matrix[i * cols + i] += muValue;
			
			std::vector<StateType> projHt;
			std::vector<double> jumpEnergy;
			tJumpEnergy(this->getAt(i), &projHt, &jumpEnergy, &this->sys_hubP);
			for (unsigned int j = 0; j < projHt.size(); j++) {
				StateType index;
				if (!this->whereIsElement(projHt.at(j), &index)) continue;
				if (index < i) continue; 
				result_matrix[i * cols + index] += jumpEnergy.at(j);
				result_matrix[index * cols + i] += jumpEnergy.at(j);
			}
		}
	}

	void H (VectorType* h_phi_n, VectorType* phi_n) { 
		/***************************************************************
		Applies H on the given vector without calculating the H matrix for the r-basis

		Parameters
		----------
		h_phi_n : (double *) receptacle for the projected vector of Lanczos algorithm
		phi_n : (double *) current vector of Lanczos algorithm
		states : (StatesArr) States used in the subspace

		Returns
		--------
		NONE
		*****************************************************************/
		int elements =  this->getLength();
		//printf("Elemnt tot = %d",elements);
		//fflush(stdout);

		//Mu value (const for every state)
		float muValue = computeMu(this->sys_hubP.mu,this);
		#pragma omp parallel for default(none) shared(elements,muValue,h_phi_n,phi_n,stdout) 
		for(int i = 0; i < elements; i++){
			//printf("Elemnt i = %d\t, states[%d]:%ld\n",i,i,states->getAt(i));
			//fflush(stdout);
			//Hu and hmu
			h_phi_n[i] += (double)(Hu(this->getAt(i), this->sys_hubP.n_sites)*this->sys_hubP.u + muValue) * phi_n[i];
		
			//printf("HU HMu done without problems\n");
			//fflush(stdout);
			//Ht
			std::vector<StateType> proj;
			std::vector<double> energies;
			tJumpEnergy(this->getAt(i), &proj, &energies, &this->sys_hubP);
			//printf("tJump size:%ld\n",proj.size());
			//fflush(stdout);
			for (unsigned int j = 0; j < proj.size(); j++) {
				StateType index;
				//printf("\tProj:%ld\n",proj.at(j));
				//fflush(stdout);
				if (this->whereIsElement(proj.at(j), &index)) {
					h_phi_n[i] += energies.at(j) * phi_n[index];
				}
			}
		}
	}


	//Pure-Virtual function override
	void add(StateType el){
		if (this->arr.capacity() == this->arr.size()) allocateMoreNodes();
		this->arr.push_back(Node(el));
		bool alreadyThere = false;
		insert(this->arr.data(), this->arr.data() + this->arr.size() - 1, &alreadyThere);

		if (alreadyThere && this->arr.size() > 1) this->arr.pop_back();
	}
	virtual bool whereIsElement(StateType el, StateType* index) const {
		/***************************************
		Searches the index of a given state if the array has it

		Parameters:
		-----------
		el : (StateType) state to look for
		index : (StateType) where is the given state in the array

		Returns:
		--------
		foudn : (bool) has the element been found
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

	std::string showAllStatesString() const{
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

		std::string allStates = "R:[";
		for (StateType i = 0; i < this->arr.size(); i++){
			allStates += std::to_string(this->arr[i].key) + t;
		}
		allStates += "]\n";

		return allStates;
	}
};
#endif
