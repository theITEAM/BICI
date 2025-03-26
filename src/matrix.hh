/// Matrix functions 

#pragma once

#include <vector>

using namespace std;

vector < vector <double> > invert_matrix(const vector <vector <double> > &Mst);
vector < vector <double> > invert_matrix_square_root(const vector < vector <double> > &M);
vector < vector <double> > invert_matrix_sparse(vector <vector <double> > M);
void print_matrix(string name, const vector < vector <double> > &M);
void print_vector(string name, const vector <double> &v);
vector <double> matrix_mult(const vector < vector <double> > &M, const vector <double> &vec);
vector <double> matrix_mult(const vector <double> &vec, const vector < vector <double> > &M);
vector < vector <double> > matrix_mult(const vector < vector <double> > &M1, const vector < vector <double> > &M2);
vector < vector <double> > calculate_cholesky(const vector < vector <double> > &M, bool &illegal);
double determinant_fast(const vector < vector <double> > &a);
double determinant_sparse(const vector < vector <double> > &a);
void add_mat_element(unsigned int i, vector <double> &Mjj, vector <unsigned int> &M_elejj,  vector <unsigned int> &M_mapjj);
void rem_mat_element(unsigned int i, vector <double> &Mjj, vector <unsigned int> &M_elejj,  vector <unsigned int> &M_mapjj);
void check_sparse(string name, const vector < vector <double> > &M, const vector < vector <unsigned int> > &M_ele, const vector < vector <unsigned int> > &M_map);
vector < vector <double> > transpose(const vector < vector <double> > &M);
double dot_prod(const vector <double> &vec1, const vector <double> &vec2);
double sparsity(const vector < vector <double> > &a);
vector <double> sample_mvn(const vector < vector <double> > &Z);
void tidy(vector < vector <double> > &M);
vector < vector <double> > create_two_by_two(const vector < vector <double> > &M, unsigned int i, unsigned int j);

