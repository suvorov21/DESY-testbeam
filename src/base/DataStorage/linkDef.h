#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#include <vector>
#pragma link C++ class vector<int>;
#pragma link C++ class vector<vector<int> > ;
#pragma link C++ class vector<short>;
#pragma link C++ class vector<vector<short> > ;
#pragma link C++ class vector<vector<double> >;

#pragma link C++ class THit+;
#pragma link C++ class std::vector<THit>+;
#pragma link C++ class std::vector<std::vector<THit> >+;
#pragma link C++ class TTrack+;
#pragma link C++ class std::vector<TTrack>+;
#pragma link C++ class TEvent+;

#endif