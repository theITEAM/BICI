#pragma once

#include <vector>

using namespace std;

const unsigned int HASH_MAX = 10000;       // The maximum hash number
const unsigned int HASH_CODE_MAX = 10000;  // The maximum number of hash codes

struct HashValue
{
	unsigned int num;                        // The number (e.g. of the equation)
	vector <unsigned int> vec;               // The vector specifying the number
};

class Hash                                 // Stores information about the model
{
public:
	Hash();
	unsigned int find(string te) const;
	unsigned int existing(const vector <unsigned int> &vec) const;
	void add(unsigned int num, const vector <unsigned int> &vec);
	void add(unsigned int num, const string &te);
	vector <unsigned int> get_vec_eqn(const string &te, unsigned int type, unsigned int p, unsigned int cl, bool infection_trans, unsigned int ti) const;
	vector <unsigned int> get_vec_string(const string &te) const;
	vector <unsigned int> get_vec_double(const double val) const;
	void create(const vector <string> &v);

private:
	vector <long> hash_code;                 // Stores codes used to make number

	vector < vector <HashValue> > table;     // For equations stores hash results
	
	unsigned int get_model_hash_number(const vector <unsigned int> &vec) const;
};

