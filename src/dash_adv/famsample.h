#ifndef FAMSAMPLE_H
#define FAMSAMPLE_H

#include <string>
#include <list>

using namespace std;

/**	Maintains information about a single individual Genotype. */
class FamSample
{
public:
	FamSample( string l ) { line = l; }
	void addCluster( int csid ) { if ( cs.back().first == csid ) cs.back().second = true; else cs.push_back( make_pair(csid,false) ); }
	list< pair<int,bool> > cs;
	string line;
};

#endif
