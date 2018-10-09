#include "tauAnalyzer.h"
R__LOAD_LIBRARY(tauAnalyzer.C+)

void run( const char* name )
{
   gSystem->CompileMacro("tauAnalyzer.C");

   TChain chain("tree");
   chain.Add(name);
   tauAnalyzer t(&chain);
   t.Loop();
}
