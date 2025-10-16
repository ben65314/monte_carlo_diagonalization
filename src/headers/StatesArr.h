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
	void remove_all(){
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

		for (unsigned int i = 0; i < addedArray_1->get_length(); i++) {
			add(addedArray_1->get_at(i));
		}
	}	

	//Accessing
	virtual StateType get_at(StateType index) const = 0;
	virtual StateType get_length() const{
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
	virtual bool countains_element(StateType el) const{
		/***************************************
		* Checks if a specified state is a part of the array (checks the nature 
        * and the number of the state)
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
		bool countains = where_is_element(el, &index);
		return countains;
	}
	virtual bool where_is_element(StateType el, StateType* index) const = 0;
	void get_num_array(std::vector<StateType>* receptacle) const{
		/***************************************
		* Returns the vector countaining only the number of the states in 
        * the array will not sort the states.
		*
		* Parameters
		* ----------
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
	void show_all_states() const{
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
		std::cout << show_all_states_string();
	}	
	//void showAllStatesInorder() const;
	virtual std::string show_all_states_string() const = 0;

	//Matrices operations
	virtual void matrix_creation(double* result_matrix) = 0;
	virtual void H(VectorType* h_phi_n, VectorType* phi_n) = 0;

	virtual void s_matrix_creation(VectorType* s_matrix) {
		make_identity(s_matrix,this->get_length());
	}

	//Sampling
	virtual void sampling_MH() = 0;
	virtual void sampling_least_energy() = 0;

	virtual void subspace_condition_expanding() = 0;

};




