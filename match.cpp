#include "match.h"

Match::Match( size_t i[], unsigned long p[] )
{
	id[0] = i[0];
	id[1] = i[1];
	pos[0] = p[0];
	pos[1] = p[1];
}

unsigned long Match::getLength()
{
	return pos[1] - pos[0];
}

unsigned long Match::getPosition( short p )
{
	return pos[p];
}

size_t Match::getSampleID( short s )
{
	return id[s];
}
