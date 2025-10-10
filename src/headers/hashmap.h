#pragma once
#include "basicFunctions.h"
class Hashmap {
private :
	bool blankHashmap = false;
	sType* hash_array;
	sType hash_size;

	//Used for the hash function
	unsigned int up, down, sites;
	std::vector<sType> combMarkersUp;
	std::vector<sType> combMarkersDown;
	sType combinationUp;

	void create_comb_markers();
	sType hash_function(sType item) const;
	//unsigned long hashStepFunction(sType item);
public:
	bool hash_setted = false;

	Hashmap(sType* arrayToHash);
	Hashmap();
	~Hashmap();
	
	void set_hash_parameters(uChar sites, uChar nUp, 
                             uChar nDown, uLong size = 0);
	void hash_set(sType item);
	bool hash_find(sType item, sType* index) const;
	std::vector<sType> hashToVector();
	void reset(sType size);
};
