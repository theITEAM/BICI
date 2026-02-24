#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <stdexcept>
#include "stdlib.h"
#include "math.h"
#include <sys/stat.h>
#include <cstring>
#include <signal.h>
#include <math.h>
#include <algorithm>

#include "lzw.hh"
#include "utils.hh"

using namespace std;

const vector <string> dic_list = {" ","!","\"","#","$","%","&","'","(",")","*","+",",","-",".","/","0","1","2","3","4","5","6","7","8","9",":",";","<","=",">","?","@","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","[","\\","]","^","_","`","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","{","|","}","~","\t","\n","\r","α","β","γ","Γ","δ","Δ","ε","ζ","η","Η","θ","Θ","ι","κ","λ","Λ","μ","ν","ξ","Ξ","ο","π","Π","ρ","σ","τ","υ","φ","Φ","χ","ψ","Ψ","ω","Ω","Σ","∫","→","〈","〉"};

const auto COMPRESS_NUM_MAX = 16000u;              // The maximum number used for compression
const auto COMPRESS_DIC_MAX = 512000000u;          // The maximum number of dictionary items           

/// Compresses a string using the LZW algorithm
string encode(string &te)
{
	// Looks to see if lines are repeated in the text
	if(true){
		auto lines = split_basic(te,'\n');
		
		Hash hash;
		for(auto li = 0u; li < lines.size(); li++){
			auto &line = lines[li];

			if(line.length() > 4){
				auto li_last = hash.find(line);
				if(li_last != UNSET){
					stringstream ss; ss << fixed; ss << li-li_last;
					auto lin = "@"+ss.str();
					if(lin.length() < line.length()) line = lin;
				}
				else{
					hash.add(li,line);
				}
			}
		}
		
		te = "";
		for(auto j = 0u; j < lines.size(); j++){
			te += lines[j];
			if(j+1 < lines.size()) te += endli;
		}
	}

	// Generates a map from unicode to dictionary number
	
	vector <unsigned int> dic_num;
	
	auto max = 0u;
	for(const auto &s :	dic_list){
		auto i = 0u;
		auto c = UTF8_to_int(i,s);
		dic_num.push_back(c);
		if(c > max) max = c;
	}
	max++;
	
	vector <unsigned int> map(max,UNSET);
	
	for(auto i = 0u; i < dic_list.size(); i++) map[dic_num[i]] = i;
		
	// Next convert from string to a list of dictionary item
	vector <unsigned int> list;
	{
		auto i = 0u;
		while(i < te.length()){
			auto c = UTF8_to_int(i,te);
			if(c >= max || map[c] == UNSET) emsg("Encoding was not possible");
			list.push_back(map[c]);
		}
	}

	Hash hash;
	vector < vector <unsigned int> > dic;
	
	for(auto i = 0u; i < dic_list.size(); i++){
		vector <unsigned int> vec; 
		vec.push_back(i);
	
		hash.add(dic.size(),vec);
		dic.push_back(vec);
	}
	
	vector <unsigned int> output;
	
	auto imax = list.size();
	auto i = 0u;
	while(i < imax){
		auto jlast = UNSET;
		auto fl = false;
		vector <unsigned int> vec;
		auto j = 0u;
		auto hash_num = 0u;
		for(auto len = 1u; i+len <= imax; len++){ 	
			auto val = list[i+len-1];
			hash_num = (hash_num+11)*2+val*(j+97);
			vec.push_back(val); j++;
			if(hash_num > HASH_NUM_MAX) hash_num = hash_num%HASH_NUM_MAX;
		
			auto j = hash.existing_with_num(vec,hash_num);
			if(j == UNSET){
				if(dic.size() < COMPRESS_DIC_MAX){
					hash.add_with_num(dic.size(),vec,hash_num);
					dic.push_back(vec);
				}
				i += len-1;
				fl = true;
				break;
			}
			
			jlast = j;
		}
		if(fl == false) i = imax;
		
		output.push_back(jlast);
	}

	string out;
	for(auto i = 0u; i < output.size(); i++){
		auto num = output[i]+35;
		if(num >= 92) num++;
		
		if(num >= COMPRESS_NUM_MAX){ // For large numbers unicode is used
			auto mult = (num/COMPRESS_NUM_MAX)-1;
			out += int_to_UTF8(COMPRESS_NUM_MAX + num%COMPRESS_NUM_MAX);
			auto num2 = mult+35; if(num2 >= 92) num2++;
			out += int_to_UTF8(num2);
		}
		else{
			out += int_to_UTF8(num);
		}
	}
	out += endli;
	
	//cout << (unsigned int)((100.0*out.length())/len_start) << "% compression ratio" << endl;
	
	if(true){ // Tests to see if decode the same
		auto te_out = decode(output);
			
		if(te != te_out){
			for(auto i = 0u; i < te.length(); i++){
				if(te.substr(i,1) !=  te_out.substr(i,1)){ emsg("we");}
			}
			emsg("string is wrong");
		}
	}
	
	return out;
}


/// Decodes results from the file
void decode_lines(vector <string> &lines)
{
	if(lines.size() != 1) emsg("Should only be one line");
	
	auto te = lines[0];
	auto len = te.length();
	
	vector <unsigned int> output;
	
	auto i = 0u;
	while(i < len){
		auto c = UTF8_to_int(i,te);
		if(c >= COMPRESS_NUM_MAX){
			auto c2 = UTF8_to_int(i,te);
			if(c2 > 92) c2--;
			c += COMPRESS_NUM_MAX*(c2-35);
		}
		
		if(c != 10){
			if(c > 92) c--;
			output.push_back(c-35);
		}
	}
	
	auto res = decode(output);
	
	lines = split_basic(res,'\n');
	
	// Adds in any repeated lines
	{
		for(auto li = 0u; li < lines.size(); li++){
			if(lines[li].length() > 1){
				if(lines[li].substr(0,1) == "@"){
					auto num = number(lines[li].substr(1));
					if(num != UNSET){
						lines[li] = lines[li-num];
					}
				}
			}
		}		
	}
}


/// Compresses a string using the LZW algorithm
string decode(const vector <unsigned int> &output)
{
	vector < vector <unsigned int> > dic;
	
	for(auto i = 0u; i < dic_list.size(); i++){
		vector <unsigned int> vec; 
		vec.push_back(i);
		dic.push_back(vec);
	}
	
	string te;
	
	auto imax = output.size();
	
	vector <unsigned int> prev;
	for(auto i = 0u; i < imax; i++){
		auto val = output[i];
			
		vector <unsigned int> W;
		if(val < dic.size()) W = dic[val];
		else{ W = prev; W.push_back(prev[0]);}
		
		for(auto j : W) te += dic_list[j];
	
		if(prev.size() != 0){
			auto vec = prev;
			vec.push_back(W[0]);
			
			if(dic.size() < COMPRESS_DIC_MAX){
				dic.push_back(vec);
			}
		}
		prev = W;
	}
	
	return te;
}


/// Converts from an unsigned int to a UTF8
string int_to_UTF8(unsigned int codepoint)
{
	string out;

	if (codepoint <= 0x7f) out.append(1, static_cast<char>(codepoint));
  else if (codepoint <= 0x7ff){
		out.append(1, static_cast<char>(0xc0 | ((codepoint >> 6) & 0x1f)));
		out.append(1, static_cast<char>(0x80 | (codepoint & 0x3f)));
	}
	else if (codepoint <= 0xffff){
		out.append(1, static_cast<char>(0xe0 | ((codepoint >> 12) & 0x0f)));
		out.append(1, static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
		out.append(1, static_cast<char>(0x80 | (codepoint & 0x3f)));
	}
	else{
		out.append(1, static_cast<char>(0xf0 | ((codepoint >> 18) & 0x07)));
		out.append(1, static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f)));
		out.append(1, static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
		out.append(1, static_cast<char>(0x80 | (codepoint & 0x3f)));
	}
	return out;
}


/// Converts from UTF8 to int
unsigned int UTF8_to_int(unsigned int &i, const string &s)
{
	auto n = (int)s[i];
	i++;
	
	if(n >= 0) return n;

	auto n2 = (int)s[i];
	i++;

	if((n & 32) == 0){
		return ((n & 31) << 6) + (n2 & 63);
	}
	else{
		auto n3 = (int)s[i];
		i++;
	
		return ((n & 15) << 12) + ((n2 & 63) << 6)+ (n3 & 63);
	}	

	return UNSET;
}


/// Test to see if working
void test()
{
	if(false){
		for(auto i = 0u; i < 32000; i++){
			auto s = int_to_UTF8(i);

			auto ii = 0u;
			auto j = UTF8_to_int(ii,s);
			if(i != j) emsg("prob");
			cout << i << " " << j << " check" << endl;
		}	
	}
	
	string te = "This is a string";
	
	auto enc = encode(te);
	
	vector <string> lines;
	lines.push_back(enc);
	
	decode_lines(lines);
}
