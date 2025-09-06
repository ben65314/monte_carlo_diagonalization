#include "hashmap.h"

Hashmap::Hashmap(){}
Hashmap::Hashmap(sType* array_to_hash) {
	hash_array = array_to_hash;
}

Hashmap::~Hashmap(){}

void Hashmap::set_hash_parameters(uChar sites, uChar n_up, uChar n_down, uLong size) {
	/*******************************************************
	* Sets the parameters needed for the hash function
	*
	* Parameters
	* ----------
	* sites : (uChar) number of sites of the system
	* n_up	: (uChar) number of up electrons in the bloc used
	* n_down: (uChar) number of down electrons in the bloc used
	* size	: (uLong) size of the hashmap
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	this->up = n_up;
	this->down = n_down;
	this->sites = sites;

	combinationUp = comb(this->sites,this->up);
	if (size == 0) hash_size = comb(this->sites,this->down) * combinationUp;
	else hash_size = size;

	//Creates the markers for the hash function
	combMarkersUp = std::vector<sType>(sites * up, 0);
	combMarkersDown = std::vector<sType>(sites * down, 0);
	create_comb_markers();
}



sType Hashmap::hash_function(sType item) const {
	/*******************************************************
	* Computes the key for a specific item
	*
	* Parameters
	* ----------
	* item		: (sType) item to add
	*
	* Returns
	* -------
	* hashkey	: (sType) key hashed with the item given
	*******************************************************/
	sType checkDown = 1<<(sites);
	sType checkUp = checkDown<<(sites);
	unsigned char upCount = 0, downCount = 0;
	sType adderUp, adderDown, indexUp = 0, indexDown = 0;

	//Checks where are the up and down electrons to get the markers information to generate the hashkey
	for (unsigned char i = 0; i < sites; i++) {
		checkUp>>=1;
		checkDown>>=1;

		bool locateUp = ((item & checkUp) != 0);
		bool locateDown = ((item & checkDown) != 0);

		//Ups
		if (!locateUp && upCount < up){
			if (up > 0) adderUp = combMarkersUp[upCount * sites + i];
			else adderUp = 1;
			indexUp += adderUp;
		}
		else upCount++;

		//Downs
		if (!locateDown && downCount < down){
			if (down > 0) adderDown = combMarkersDown[downCount * sites + i];
			else adderDown = 1;
			indexDown += adderDown;
		}
		else downCount++;
	}
	sType hashKey = combinationUp * indexDown + indexUp;


	return hashKey;
}

void Hashmap::hash_set(sType item) {
	/*******************************************************
	* Adds an element to the hashmap
	*
	* Parameters
	* ----------
	* item		: (sType) item to add
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	//Get the hashvalue
	sType hashValue = hash_function(item);
	hashValue %= hash_size;

	//Resolves hash conflict
	unsigned int conflict = 0;
	while(hash_array[hashValue] != 0 && hash_array[hashValue] != item) {
		std::cout<<"STUCKED?\t hashValue = "<<hashValue<<"\thash_size = "<<hash_size<<"\titem = "<< item << "\tcurrent Occupant = "<<hash_array[hashValue]<<std::endl;
		hashValue ++;
		conflict ++;
		if(hashValue >= hash_size){
			hashValue %= hash_size;
		}
	}
	averageConflict += conflict;
	if (conflict > maxConflict) maxConflict = conflict;
	//Set item in the hashmap
	hash_array[hashValue] = item;
}

bool Hashmap::hash_find(sType item, sType* index) const {
	/*******************************************************
	* Searches an element in the hashmap
	*
	* Parameters
	* ----------
	* item	: (sType) item searched
	* index	: (sType) the index where the item is placed (if found)
	*
	* Returns
	* -------
	* true	: if the element was found else false 
	*******************************************************/
	//Get hashvalue
	sType hashValue = hash_function(item);
	hashValue %= hash_size;

	//int move= 0;
	unsigned int displacement = 0;
	//Look for possible hashvalue shifting
	while (hash_array[hashValue] != 0) {
		if(hash_array[hashValue] == item){
			/*if(move>1000){std::cout<<"+:"<<move<<std::endl;}*/ 
			*index = hashValue;
			return true;
		}
		displacement ++;
		hashValue ++;
		if(hashValue >= hash_size){
			hashValue %= hash_size;
		}
		if (displacement > maxConflict) break;
		

	}
	//if(move>0){std::cout<<"+:"<<move<<std::endl;}
	return false;
}

std::vector<sType> Hashmap::hashToVector() {
	/*******************************************************
	* Creates a vector of all the elements in the hashmap
	*
	* Parameters
	* ----------
	* NONE
	*
	* Returns
	* -------
	* vecHash : (std::vector<sType>) vector of the elements
	*******************************************************/
	std::vector<sType> vecHash;
	vecHash.reserve(hash_size);

	for	(sType i = 0; i < hash_size; i++){
		if(hash_array[i] != 0){
			vecHash.push_back(hash_array[i]);
		}
	}
	return vecHash;
}

void Hashmap::reset(sType size) {
	/*******************************************************
	* Deletes the hasmap and creates another one of size <size>
	*
	* Parameters
	* ----------
	* size	: (sType) size of the new hashmap
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	delete[] hash_array;
	hash_size = size;
	hash_array = new sType[hash_size]();


}

void Hashmap::create_comb_markers() {
	/*******************************************************
	* Creates the markers used in the hashfunction
	*
	* Parameters
	* ----------
	* NONE
	*
	* Returns
	* -------
	* NONE
	*******************************************************/
	//Up comb markers
	for (unsigned char i = 0; i < up; i++) { //Rows
		for (unsigned char j = 0; j < sites; j++) {  //Columns
			if (sites - j < up - i) {
				combMarkersUp[i * sites + j] = 1;
			}
			combMarkersUp[i * sites + j] = comb(sites - j - 1, up - i - 1);
		}
	}

	//Down comb markers
	for (unsigned char i = 0; i < down; i++) { //Rows
		for (unsigned char j = 0; j < sites; j++) {  //Columns
			if (sites - j < down - i) {
				combMarkersDown[i * sites + j] = 1;
			}
			combMarkersDown[i * sites + j] = comb(sites - j - 1, down - i - 1);
		}
	}
}

