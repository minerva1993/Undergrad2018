#include "tauAnalyzer.h"
R__LOAD_LIBRARY(tauAnalyzer.C+)

void run( const char* input_name, const std::string output_name )
{
   gSystem->CompileMacro("tauAnalyzer.C");

   TChain chain("tree");
   chain.Add(input_name);
   tauAnalyzer t(&chain);
   t.Loop(output_name);
}
