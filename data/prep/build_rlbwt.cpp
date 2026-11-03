// Small/reference-scale construction baseline. Large panels require a scalable builder.
#include <divsufsort64.h>
#include <fstream>
#include <vector>
#include <iterator>
#include <stdexcept>
int main(int argc,char**argv){
 if(argc!=3)return 2;
 std::ifstream in(argv[1],std::ios::binary);
 if(!in)throw std::runtime_error("missing text");
 std::vector<unsigned char> t((std::istreambuf_iterator<char>(in)),{});
 if(t.empty()||t.back()!=10)throw std::runtime_error("text must end in one newline sentinel");
 std::vector<int64_t> sa(t.size());
 if(divsufsort64(t.data(),sa.data(),t.size()))return 1;
 std::ofstream heads(std::string(argv[2])+".bwt.heads",std::ios::binary),lens(std::string(argv[2])+".bwt.len",std::ios::binary);
 unsigned char last=0;uint64_t run=0;
 auto flush=[&](){if(!run)return;heads.put(last);for(int k=0;k<5;k++)lens.put((run>>(8*k))&255);};
 for(auto pos:sa){auto c=t[pos?pos-1:t.size()-1];if(run&&c!=last){flush();run=0;}last=c;++run;}flush();
 return !(heads.good()&&lens.good());
}
