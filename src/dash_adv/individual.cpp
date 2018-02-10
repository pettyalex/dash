#include "individual.h"

Individual::Individual( size_t ctr , string i , FamSample * fs )
{
	nid = ctr;
	id = i;
	fams = fs;
}

size_t Individual::getNumericID()
{
	return nid;
}

string Individual::getID()
{
	return id;
}

void Individual::addCluster( unsigned int csid )
{
	if ( fams != NULL ) fams->addCluster( csid );
}
