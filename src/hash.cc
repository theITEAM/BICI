/// Implent hash tables for finding items on a list in a faster way

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>    
 
using namespace std;

#include "hash.hh"
#include "utils.hh"

Hash::Hash()
{
	table.resize(HASH_MAX);
	hash_code.resize(HASH_CODE_MAX);
	
	for(auto h = 0u; h < HASH_CODE_MAX; h++){
		hash_code[h] = (unsigned int)(ran()*HASH_MAX);
	}
};	
	
	
/// Determines if the state is already in the hash table
unsigned int Hash::find(string te) const 
{
	auto vec = get_vec_string(te);
	return existing(vec);
}

	
/// Determines if the state is already in the hash table
unsigned int Hash::existing(const vector <unsigned int> &vec) const 
{
	auto num = get_model_hash_number(vec);
		
	for(const auto &ele : table[num]){
		if(equal_vec(ele.vec,vec) == true) return ele.num;
	}
	
	return UNSET;
}
	

/// For equations converts from the equation information to a list of integers
vector <unsigned int> Hash::get_vec_eqn(const string &te, unsigned int type, unsigned int p, unsigned int cl, bool infection_trans, unsigned int ti, const string &c_te) const
{
	vector <unsigned int> vec;
	
	for(auto i = 0u; i < te.length(); i++){
		int ch = te.at(i);
		if(ch < 0) ch += 256;
		vec.push_back(ch);
	}
	
	for(auto i = 0u; i < c_te.length(); i++){
		int ch = c_te.at(i);
		if(ch < 0) ch += 256;
		vec.push_back(ch);
	}
	
	vec.push_back(type);
	vec.push_back(p);
	vec.push_back(cl);
	vec.push_back(ti);
	if(infection_trans == 1) vec.push_back(1); else vec.push_back(0);
	
	return vec;
}


/// COnverts from a text string to a list of integers
vector <unsigned int> Hash::get_vec_string(const string &te) const
{
	vector <unsigned int> vec;
	
	for(auto i = 0u; i < te.length(); i++){
		int ch = te.at(i);
		if(ch < 0) ch += 256;
		vec.push_back(ch);
	}
	
	return vec;
}


/// COnverts from a text string to a list of integers
vector <unsigned int> Hash::get_vec_double(const double val) const
{
	stringstream ss;
	ss << fixed << std::setprecision(16);
	ss << val;
	
	return get_vec_string(ss.str());
}


/// Calculates a hash number from a vector
unsigned int Hash::get_model_hash_number(const vector <unsigned int> &vec) const 
{
	auto num = 0u;
	
	auto i = 0u;
	for(auto j = 0u; j < vec.size(); j++){
		i = (i+256+vec[j])%HASH_CODE_MAX; num += hash_code[i];
	}
	
	return num%HASH_MAX;
}


/// Adds a new element to the hash table
void Hash::add(unsigned int num, const vector <unsigned int> &vec)
{
	HashValue hv; hv.num = num; hv.vec = vec;
	
	table[get_model_hash_number(vec)].push_back(hv);
}


/// Adds a new element to the hash table
void Hash::add(unsigned int num, const string &te)
{
	auto vec = get_vec_string(te);
	
	HashValue hv; hv.num = num; hv.vec = vec;
	
	table[get_model_hash_number(vec)].push_back(hv);
}

/// Creates a hash table
void Hash::create(const vector <string> &v)
{
	for(auto i = 0u; i < v.size(); i++){		
		auto vec = get_vec_string(v[i]);
		add(i,vec);
	}
}
