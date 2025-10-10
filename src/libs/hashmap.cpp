#include "hashmap.h"

Hashmap::Hashmap(){}
Hashmap::Hashmap(sType* array_to_hash) {
	this->hash_array = array_to_hash;
}

Hashmap::~Hashmap(){}

void Hashmap::set_hash_parameters(uChar sites, uChar n_up, 
                                  uChar n_down, uLong size) {
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

	this->combinationUp = comb(this->sites, this->up);
	if (size == 0) {
        this->hash_size = comb(this->sites, this->down) * this->combinationUp;
    }
	else {
        this->hash_size = size;
    }

	// Creates the markers for the hash function
	this->combMarkersUp = std::vector<sType>(this->sites * this->up, 0);
	this->combMarkersDown = std::vector<sType>(this->sites * this->down, 0);
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
    // Inspection variables
	sType checkDown = 1<<(this->sites);
	sType checkUp = checkDown<<(this->sites);
	unsigned int upCount = 0, downCount = 0;
	sType adderUp, adderDown, indexUp = 0, indexDown = 0;

	// Checks where are the up and down electrons to get the markers 
    // information to generate the hashkey
	for (unsigned char i = 0; i < this->sites; i++) {
		checkUp>>=1;
		checkDown>>=1;
        
        // Found electron
		bool locateUp = ((item & checkUp) != 0);
		bool locateDown = ((item & checkDown) != 0);

		// Ups
		if (!locateUp && upCount < this->up){
			if (this->up > 0) {
                adderUp = this->combMarkersUp[upCount * this->sites + i];
            }
			else {
                adderUp = 1;
            }
			indexUp += adderUp;
		}
		else upCount++;

		// Downs
		if (!locateDown && downCount < this->down){
			if (this->down > 0) {
                adderDown = this->combMarkersDown[downCount * this->sites + i];
            }
			else {
                adderDown = 1;
            }
			indexDown += adderDown;
		}
		else downCount++;
	}
	sType hashKey = this->combinationUp * indexDown + indexUp;

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
	// Get the hashvalue
	sType hashValue = hash_function(item);
	hashValue %= this->hash_size;

	// Resolves hash conflict
	while(this->hash_array[hashValue] != 0 
          && this->hash_array[hashValue] != item) {
        if (verbose > 10) std::cout << "conflict ";
		hashValue ++;

		if (hashValue >= this->hash_size) {
			hashValue %= this->hash_size;
		}
	}
	// Set item in the hashmap
	this->hash_array[hashValue] = item;
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
	// Get hashvalue
	sType hashValue = hash_function(item);
	hashValue %= this->hash_size;

	// Look for possible hashvalue shifting
	while (this->hash_array[hashValue] != 0) {
		if (this->hash_array[hashValue] == item) {
			*index = hashValue;
			return true;
		}
		hashValue ++;
		if (hashValue >= this->hash_size) {
			hashValue %= this->hash_size;
		}
	}
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
	vecHash.reserve(this->hash_size);

	for	(sType i = 0; i < this->hash_size; i++) {
		if (this->hash_array[i] != 0) {
			vecHash.push_back(this->hash_array[i]);
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
	delete[] this->hash_array;
	this->hash_size = size;
	this->hash_array = new sType[this->hash_size]();


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
	// Up comb markers
	for (unsigned char i = 0; i < this->up; i++) { // Rows
		for (unsigned char j = 0; j < this->sites; j++) {  // Columns
			if (this->sites - j < this->up - i) {
				this->combMarkersUp[i * this->sites + j] = 1;
			}
			this->combMarkersUp[i * this->sites + j] 
                = comb(this->sites - j - 1, this->up - i - 1);
		}
	}

	// Down comb markers
	for (unsigned char i = 0; i < this->down; i++) { // Rows
		for (unsigned char j = 0; j < this->sites; j++) {  // Columns
			if (this->sites - j < this->down - i) {
				this->combMarkersDown[i * this->sites + j] = 1;
			}
			this->combMarkersDown[i * this->sites + j] 
                = comb(this->sites - j - 1, this->down - i - 1);
		}
	}
}

