// Allene test set: input SMILES, ChemLarp label, canonical SMILES, three random spellings.
#include <cstdio>
#include <string>
#include <vector>
#include <random>
#include "chemlarp/smiles.hpp"
#include "chemlarp/canonical.hpp"
#include "chemlarp/cip.hpp"
using namespace chemlarp;
int main(){
  std::vector<std::string> sub={"","C","CC","O","Cl","F","Br","N","C(=O)O","c1ccccc1","OC","C=C","[2H]","C#N"};
  std::mt19937 rng(7); int n=0;
  auto one=[&](const std::string& s){ ChemicalSystemState st; Diagnostics d; int sp=parse_smiles(st,s,d); if(sp<0) return; auto ax=cip_axes(st,sp); if(ax.empty()) return;
    std::string lab=to_string(ax[0].label); if(lab.empty()) return;
    printf("%s\t%s\t%s",s.c_str(),lab.c_str(),st.species[sp].canonical_smiles.c_str());
    for(int k=0;k<3;++k){ SmilesWriteOptions o; o.canonical=false; o.random_seed=1000*n+k; printf("\t%s",write_smiles(st,sp,o).c_str()); } puts(""); ++n; };
  auto grp=[&](const std::string& x){ return x.empty()? std::string(): "("+x+")"; };
  for(int t=0;t<400;++t){ std::string a=sub[rng()%sub.size()],b=sub[rng()%sub.size()],c=sub[rng()%sub.size()],e=sub[rng()%sub.size()];
    if(a.empty()) std::swap(a,b); if(c.empty()) std::swap(c,e); if(a.empty()||c.empty()) continue;
    std::string chi=rng()%2?"@":"@@";
    // forms: a C(b) = [C@] = C(e) c   and with ring-free variants
    one(a+"C"+grp(b)+"=[C"+chi+"]="+"C"+grp(e)+c);
  }
  for(const char* s: {"C1CCCC=[C@]=CCC1","C1CCCC=[C@@]=CCC1","CC1=[C@]=C(C)CCCCC1","C1=[C@]=C(Cl)CCCCC1","OC1CCCCC=[C@]=C1","CC=C=[C@]=C=CC","CC=C=[C@@]=C=C(C)Cl","[C@](=CC)=CF","C(=[C@]=CC)Br","C(F)(=[C@@]=C(O)C)Cl","C[C@H](O)C=[C@]=C[C@@H](C)O","C[C@H](O)C=[C@@]=C[C@H](C)O"}) one(s);
}
