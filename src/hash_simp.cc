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
	size = HASH_SIMP_NUM_DEF;
	table.resize(size);
	n = 0;
	
	on = true;
}


/*
/// Initialises hash table with a given size
void HashSimp::init(unsigned int num)
{
	size = num;
	table.resize(num);
}
*/


/// Turns off the hash table
void HashSimp::off()
{
	on = false;
	
	if(profiling){
		auto max = 0.0; 
		for(const auto &va : table){
			if(va.size() > max) max = va.size();  
		}
		if(max > 30){
			cout << n/size <<  " hash_simp average    max=" << max << endl;
			emsg("Hash simp too large");
		}
	}
	
	table.clear();
}


/// Determines if the state is already in the hash table
unsigned int HashSimp::find(unsigned int val) const 
{
	if(!on) emsg("Hash table turned off");
	const auto &vec = table[(val*HASH_MULT)%size];
	
	for(auto j = 0u; j < vec.size(); j++){
		if(vec[j].val == val) return vec[j].num;
	}
	
	return UNSET;
}


/// Adds an element to the hash table
void HashSimp::add(unsigned int num, unsigned int val)
{
	if(!on) emsg("Hash table turned off");
	HashSimpValue hsv; hsv.num = num; hsv.val = val;
	
	table[(val*HASH_MULT)%size].push_back(hsv);
	
	n++;
	if(n/size > HASH_OCC_THRESH) enlarge();
}


/// Enlarges the hash table
void HashSimp::enlarge()
{  
	size *= 5;
	
	auto table_old = table;
	table.clear();
	table.resize(size);
	for(const auto &row : table_old){
		for(const auto &ele : row){
			table[(ele.val*HASH_MULT)%size].push_back(ele);
		}
	}			
}

