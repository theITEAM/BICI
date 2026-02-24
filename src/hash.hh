#pragma once

#include <vector>

using namespace std;

const unsigned int HASH_MAX_DEF = 11;     // The starting hash number

struct HashValue
{
	unsigned int num;                        // The number (e.g. of the equation)
	unsigned int hash_num;                   // The hash number
	vector <unsigned int> vec;               // The vector specifying the number
};

class Hash                                 // Stores information about the model
{
public:
	Hash();
	void off();
	unsigned int find(string te) const;
	unsigned int existing(const vector <unsigned int> &vec) const;
	unsigned int existing_with_num(const vector <unsigned int> &vec, unsigned int hash_num) const;
	void add(unsigned int num, const vector <unsigned int> &vec);
	void add_with_num(unsigned int num, const vector <unsigned int> &vec, unsigned int hash_num);
	void add(unsigned int num, const string &te);
	vector <unsigned int> get_vec_eqn(const string &te, unsigned int type, unsigned int p, unsigned int cl, bool infection_trans, unsigned int ti) const;
	vector <unsigned int> get_vec_string(const string &te) const;
	vector <unsigned int> get_vec_double(const double val) const;
	void create(const vector <string> &v);
	double mem() const;

private:
	void enlarge();

	bool on;                                 // Determines if hash table turned on

	unsigned int size;                       // The size of the hash table
	
	double n;                                // The total number of elements in the hash table
	
	vector < vector <HashValue> > table;     // For equations stores hash results
	
	unsigned int get_model_hash_number(const vector <unsigned int> &vec) const;
};

