#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "famsample.h"

#include <string>
#include <list>

using namespace std;

/**	Maintains information about a single individual Genotype. */
class Individual
{
public:

	/** Constructor for Genotype
		@param i int identifier
		*/
	Individual( size_t ctr , string i , FamSample * fs );

	/** Returns the string identifier */
	string getID();
	/** Returns the numeric identifier */
	size_t getNumericID();
	/** Adds the Cluster to the FamSample of this individual, if available **/
	void addCluster( unsigned int );
private:
	/// Contains the number id of this individual
	size_t nid;

	/// Contains the string id of this individual
	string id;

	/// Pointer to the FamSample for this individual
	FamSample * fams;
};

#endif
