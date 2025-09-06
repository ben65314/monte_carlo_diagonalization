#pragma once
#include "electronManipulationFunctions.h"

template <class Search, class StateType, class VectorType> class StatesArr {
protected:
	std::vector<Search> arr;

public:
	hubbardParam sys_hubP;
	samplingParam sys_sP;
	virtual ~StatesArr() = default;

	Electrons electrons;
	virtual StatesArr<Search,StateType,VectorType>* clone() = 0;

	//Settings
	void set_hubbard_parameters(hubbardParam hubP) {sys_hubP = hubP;}
	void set_sampling_parameters(samplingParam sP) {sys_sP = sP;}

	//Manipulation
	virtual void add(StateType el) = 0;	
	void removeAll(){
		/***************************************
		* Clear all the states countained
		*
		* Parameters:
		* -----------
		* NONE
		*
		* Returns:
		* --------
		* NONE
		****************************************/
		arr.clear();
	}
	void merge(const StatesArr* addedArray_1){
		/***************************************
		* Merge a StatesArr to the current one
		* Doesn't sort them, if needed call sortStates() right after 
		*
		* Parameters:
		* -----------
		* addedArray_1 : (StatesArr) state to merge
		*
		* Returns:
		* --------
		* this
		****************************************/

		for (unsigned int i = 0; i < addedArray_1->getLength(); i++) {
			add(addedArray_1->getAt(i));
		}
	}	

	//Accessing
	virtual StateType getAt(StateType index) const = 0;
	virtual StateType getLength() const{
		/***************************************
		* Get function of the number of elements in the array
		* 
		* Parameters:
		* -----------
		* NONE
		*
		* Returns:
		* --------
		* arr.size() : (StateType) number of elements in the array
		****************************************/
		
		return arr.size();
	}
	virtual bool countainsElement(StateType el) const{
		/***************************************
		* Checks if a specified state is a part of the array (checks the nature and the number of the state)
		*
		* Parameters:
		* -----------
		* el : (StateType) state to look for
		*
		* Returns:
		* --------
		* countains : (bool) is the state given in the array
		****************************************/
		StateType index = 0;
		bool countains = whereIsElement(el, &index);
		return countains;
	}
	virtual bool whereIsElement(StateType el, StateType* index) const = 0;
	void getNumArray(std::vector<StateType>* receptacle) const{
		/***************************************
		* Returns the vector countaining only the number of the states in the array will not sort the states.
		*
		* Parameters
		* ----------
		* nature	: (char) 
		* receptacle: (vector<StateType>) vector of only the states number
		*
		* Returns
		* -------
		* NONE
		****************************************/
		if (arr.size() > 0) {
			receptacle->reserve(arr.size());
			for (uint64_t i = 0; i <arr.size(); i++) {
				receptacle->push_back(arr.at(i).key);
			}
		}
	}
	
	//Presentation
	void showAllStates() const{
		/***************************************
		* Prints all the States countained
		*
		* Parameters:
		* -----------
		* NONE
		*	
		* Returns:
		* --------
		* NONE
		****************************************/
		std::cout << showAllStatesString();
	}	
	//void showAllStatesInorder() const;
	virtual std::string showAllStatesString() const = 0;

	//Matrices operations
	virtual void matrixCreation(double* result_matrix) = 0;
	virtual void H(VectorType* h_phi_n, VectorType* phi_n) = 0;

	virtual void s_matrix_creation(VectorType* s_matrix) {
		makeIdentity(s_matrix,this->getLength());

	}

	//Sampling
	virtual void sampling_MH() = 0;
	virtual void sampling_least_energy() = 0;

	virtual void subspace_condition_expanding() = 0;

};


template <class A, class B, class C> double computeMu(float mu, StatesArr<A,B,C>* sArr){
	/**************************************************
	* Computes mu interaction used in matrixCreation
	*
	* Parameters
	* ----------
	* mu: (float) mu value
	* state: (States) state affected by the mu relation on the hamiltonian 
	*
	* Returns
	* ------ -
	* result: (double) mu value
	************************************************/

	double result;
	//Counts number of electrons
	long nElectron = sArr->electrons.up + sArr->electrons.down;
	//Computs chemical potential
	result = -nElectron * mu;
	return result;
}


