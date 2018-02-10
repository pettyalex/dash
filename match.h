#ifndef MATCH_H
#define MATCH_H

#include "individual.h"

/** Maintains information about a Match connecting two haplotypes */
class Match
{
public:
	/** Constructor
		@param i points to incident individuals
		@param p contains start and end position of match
		*/
	Match( size_t i[], unsigned long p[] );

	/** Returns one of the two connected nodes
		@param s 0/1 for which vertex
		@return id of the individual
		*/
	size_t getSampleID( short s );

	/** Returns position in match
		@param p indicates the start (0) or end (1) position to return
		@return the physical position
		*/
	unsigned long getPosition( short p );
	
	/** Returns the length of the edge */
	unsigned long getLength();
	
private:
	
	/// Physical position of the underlying match
	unsigned long pos[2];
	
	/// Integer ID of individuals in match
	size_t id[2];
};

#endif
