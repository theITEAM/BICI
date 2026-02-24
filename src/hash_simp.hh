#pragma once

#include <vector>

using namespace std;

const auto HASH_SIMP_NUM_DEF = 11;        // Default size for hash code 

struct HashSimpValue
{
	unsigned int num;                        // The number for the index
	unsigned int val;                        // The number which is stored
};

class HashSimp                                 // Stores information about the model
{
public:
	HashSimp();
	
	void off();
	unsigned int find(unsigned int val) const;
	void enlarge();
	void add(unsigned int num, unsigned int val); 
	
	bool on;
	double n;             // Number of elements
	unsigned int size;    // Size of hash table
	vector < vector <HashSimpValue> > table;
};
