/// Used for synchrosising chains

#pragma once

using namespace std;

#include "chain.hh"
#include "matrix.hh"

void synchronise_proposal(unsigned int s, vector <Chain> &chain, Mpi &mpi);

void print_prop(string te, const vector <Chain> &chain);
string print_prop_info(const vector <PropInfo> &prop_info, const Model &model);