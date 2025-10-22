/// This is for a simpler form of hash table just for unsigned int

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>    
 
using namespace std;

#include "hash_simp.hh"
#include "utils.hh"

HashSimp::HashSimp()
{
	table.resize(HASH_SIMP_NUM);
}

/// Determines if the state is already in the hash table
unsigned int HashSimp::find(unsigned int val) const 
{
	const auto &vec = table[val%HASH_SIMP_NUM];
	
	for(auto j = 0u; j < vec.size(); j++){
		if(vec[j].val == val) return vec[j].num;
	}
	
	return UNSET;
}


/// Adds an element to the hash table
void HashSimp::add(unsigned int num, unsigned int val)
{
	HashSimpValue hsv; hsv.num = num; hsv.val = val;
	
	table[val%HASH_SIMP_NUM].push_back(hsv);
}
