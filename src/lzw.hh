/// Used for compressing and decompressing text strings

#pragma once

using namespace std;

#include "struct.hh"

string encode(string &te);
void decode_lines(vector <string> &lines);
string decode(const vector <unsigned int> &output);
string int_to_UTF8(unsigned int codepoint);
unsigned int UTF8_to_int(unsigned int &i, const string &s);
void test();
