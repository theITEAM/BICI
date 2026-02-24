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
	on = true;
	size = HASH_MAX_DEF;
	n = 0;
	table.resize(size);
};	
	
	
/// Turns off the hash table
void Hash::off()
{
	on = false;
	
	if(profiling){
		auto max = 0u; 
		for(const auto &va : table){
			if(va.size() > max) max = va.size();
		}
		if(max > 30){
			cout << n/size << " hash average op    max = " << max << endl;
			emsg("Hash too large");
		}
	}
	
	table.clear();
}


/// Determines if the state is already in the hash table
unsigned int Hash::find(string te) const 
{
	if(!on) emsg("Hash table turned off");
	auto vec = get_vec_string(te);
	return existing(vec);
}

	
/// Determines if the state is already in the hash table
unsigned int Hash::existing(const vector <unsigned int> &vec) const 
{
	if(!on) emsg("Hash table turned off");
	auto hash_num = get_model_hash_number(vec);
		
	for(const auto &ele : table[hash_num%size]){
		if(ele.hash_num == hash_num){
			if(equal_vec(ele.vec,vec) == true) return ele.num;
		}
	}
	
	return UNSET;
}
	

/// Determines if the state is already in the hash table (with hash num alread calculated)
unsigned int Hash::existing_with_num(const vector <unsigned int> &vec, unsigned int hash_num) const 
{
	if(!on) emsg("Hash table turned off");
	
	if(false){ 
		auto hash_num_ch = get_model_hash_number(vec);
		if(hash_num != hash_num_ch) emsg("hash number error");
	}
	
	for(const auto &ele : table[hash_num%size]){
		if(ele.hash_num == hash_num){
			if(equal_vec(ele.vec,vec) == true) return ele.num;
		}
	}
	
	return UNSET;
}

/// For equations converts from the equation information to a list of integers
vector <unsigned int> Hash::get_vec_eqn(const string &te, unsigned int type, unsigned int p, unsigned int cl, bool infection_trans, unsigned int ti) const
{
	if(!on) emsg("Hash table turned off");
	
	vector <unsigned int> vec;
	
	for(auto i = 0u; i < te.length(); i++){
		int ch = te.at(i);
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
	if(!on) emsg("Hash table turned off");
	
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
	if(!on) emsg("Hash table turned off");
	
	stringstream ss;
	ss << fixed << std::setprecision(16);
	ss << val;
	
	return get_vec_string(ss.str());
}


/// Calculates a hash number from a vector
unsigned int Hash::get_model_hash_number(const vector <unsigned int> &vec) const 
{
	if(!on) emsg("Hash table turned off");
	
	auto num = 0u;
	for(auto j = 0u; j < vec.size(); j++){
		num = (num+11)*2+vec[j]*(j+97);
		if(num > HASH_NUM_MAX) num = num%HASH_NUM_MAX;
	}		
	
	return num;
}


/// Enlarges the hash table
void Hash::enlarge()
{  
	size *= HASH_ENLARGE_SIZE;
	
	auto table_old = table;
	table.clear();
	table.resize(size);
	for(const auto &row : table_old){
		for(const auto &hv : row){
			table[hv.hash_num%size].push_back(hv);
		}
	}			
}
	
		
/// Adds a new element to the hash table
void Hash::add(unsigned int num, const vector <unsigned int> &vec)
{
	if(!on) emsg("Hash table turned off");
	
	auto hash_num = get_model_hash_number(vec);
	
	HashValue hv; hv.num = num; hv.hash_num = hash_num; hv.vec = vec;
	
	table[hash_num%size].push_back(hv);
	
	n++;
	if(n/size > HASH_OCC_THRESH) enlarge();
}


/// Adds a new element to the hash table
void Hash::add_with_num(unsigned int num, const vector <unsigned int> &vec, unsigned int hash_num)
{
	if(!on) emsg("Hash table turned off");
	
	if(false){
		auto hash_num_ch = get_model_hash_number(vec);
		if(hash_num != hash_num_ch) emsg("hash number error");
	}
	
	HashValue hv; hv.num = num; hv.hash_num = hash_num; hv.vec = vec;
	
	table[hash_num%size].push_back(hv);
	
	n++;
	if(n/size > HASH_OCC_THRESH) enlarge();
}


/// Adds a new element to the hash table
void Hash::add(unsigned int num, const string &te)
{
	if(!on) emsg("Hash table turned off");
	
	auto vec = get_vec_string(te);
	
	auto hash_num = get_model_hash_number(vec);
	
	HashValue hv; hv.num = num; hv.hash_num = hash_num;  hv.vec = vec;
	
	table[hash_num%size].push_back(hv);
	
	n++;
	if(n/size > HASH_OCC_THRESH) enlarge();
}

/// Creates a hash table
void Hash::create(const vector <string> &v)
{
	if(!on) emsg("Hash table turned off");
	
	for(auto i = 0u; i < v.size(); i++){		
		auto vec = get_vec_string(v[i]);
		add(i,vec);
	}
}


/// Gets memory used
double Hash::mem() const
{
	auto sum = 0.0;
	for(auto i = 0u; i < table.size(); i++){
		sum += sizeof(vector <HashValue>);
		for(auto j = 0u; j < table[i].size(); j++){
			sum += sizeof(HashValue)+table[i][j].vec.size()*sizeof(unsigned int);
		}
	}
	
	return sum;
}

