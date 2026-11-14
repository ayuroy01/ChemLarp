// Hindered biaryls: embed, write SDF, read back; prints label and appends the SDF (with label as title).
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include "chemlarp/conformer.hpp"
#include "chemlarp/molio.hpp"
#include "chemlarp/cip.hpp"
using namespace chemlarp;
int main(int argc, char** argv){
  std::vector<std::string> o={"C","Cl","Br","O","OC","N","F","C(F)(F)F","[N+](=O)[O-]","C(C)C","I","C(=O)O"};
  std::vector<std::string> m={"","C","Cl","OC","F"};
  std::mt19937 rng(11); FILE* f=fopen(argc > 1 ? argv[1] : "atropisomers.sdf","w"); int n=0, none=0;
  std::vector<std::string> list={"Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12","Cc1cccc(C)c1-c1c(C)cccc1Cl","CC(C)c1cccc(C(C)C)c1-c1ccccc1[N+](=O)[O-]","Oc1ccccc1-c1c(C)cccc1Br"};
  for(int t=0;t<120;++t){ std::string a=o[rng()%o.size()],b=o[rng()%o.size()],c=o[rng()%o.size()],d= rng()%3? o[rng()%o.size()]:"", x=m[rng()%m.size()];
    list.push_back("c1("+a+")c"+(x.empty()?"":"("+x+")")+"ccc("+b+")c1-c1c("+c+")cccc1"+(d.empty()?"":"("+d+")")); }
  for(auto& s: list){ ChemicalSystemState st; Diagnostics dd; int sp=embed_smiles(st,s,dd); if(sp<0) continue;
    std::string sdf=write_sdf(st,sp); ChemicalSystemState s2; Diagnostics d2; auto r=read_molecules(s2,sdf,d2); if(r.empty()||r[0].species<0) continue;
    std::string lab; for(auto&a:cip_axes(s2,r[0].species)) lab=to_string(a.label); if(lab.empty()){++none; lab="none";}
    sdf.replace(0,sdf.find('\n'),lab+" "+s); fputs(sdf.c_str(),f); ++n; }
  fclose(f); printf("%d labelled, %d without an axis\n",n,none); }
