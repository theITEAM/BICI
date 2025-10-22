#pragma once

#include <vector>

using namespace std;

const auto HASH_SIMP_NUM = 50;             // How bit the hash code is

struct HashSimpValue
{
	unsigned int num;                        // The number for the index
	unsigned int val;                        // The number which is stored
};

class HashSimp                                 // Stores information about the model
{
public:
	HashSimp();
	unsigned int find(unsigned int val) const;
	void add(unsigned int num, unsigned int val); 
	
	vector < vector <HashSimpValue> > table;
};
