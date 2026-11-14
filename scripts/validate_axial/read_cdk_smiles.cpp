#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include "chemlarp/smiles.hpp"
#include "chemlarp/cip.hpp"
using namespace chemlarp;
int main(int,char**v){ std::ifstream f(v[1]); std::string l; int n=0,ok=0;
 while(std::getline(f,l)){ std::stringstream s(l); std::string in,lab,w1,w2; std::getline(s,in,'\t'); std::getline(s,lab,'\t'); std::getline(s,w1,'\t'); std::getline(s,w2,'\t');
  for(auto& w: {w1,w2}){ ChemicalSystemState st; Diagnostics d; int sp=parse_smiles(st,w,d); std::string L; if(sp>=0) for(auto&a:cip_axes(st,sp)) L=to_string(a.label); ++n; if(L==lab) ++ok; else printf("DIFF %s ours %s | cdk wrote %s -> %s\n",in.c_str(),lab.c_str(),w.c_str(),L.c_str()); } }
 printf("%d/%d CDK-written SMILES give the same label\n",ok,n); }
