#include <iostream>
#include <sstream>

#include "main.h"
#include "match.h"
#include "pedigree.h"

#define VERBOSE 0

float MIN_CLUSTER_DENSITY = 0.6f;
float MAX_DIFFERENCE_FOR_PRINT = 0.85f;
unsigned int MIN_GRAPH_SIZE = 4;
unsigned long WINDOW_SIZE = 500000;
unsigned long cur_pos = 0;
/// File name to write cluster output
string CLUSTER_FILE = "dash";
string FAM_FILE = "";
bool HAP_MODE = false;
bool INDEP_MODE = false;
ofstream of_log;

Pedigree MAIN_PED;

/// Analyses parameters and returns weather valid or invalid
bool process_param( int argc , char * argv[] )
{
	bool good = true;
	for( int i = 1; i < argc; i++ )
	{
		if ( string(argv[i]) == "-help" ) { return false; }
		else if ( string(argv[i]) == "-haploid" ) { HAP_MODE = true; }
		else if ( string(argv[i]) == "-indep" ) { INDEP_MODE = true; }
		else if ( string(argv[i]) == "-win" ) { WINDOW_SIZE = atol( argv[++i] ); }
		else if ( string(argv[i]) == "-density" ) { MIN_CLUSTER_DENSITY = (float) atof( argv[++i] ); }
		else if ( string(argv[i]) == "-r2" ) { MAX_DIFFERENCE_FOR_PRINT = (float) atof( argv[++i] ); }
		else if ( string(argv[i]) == "-min" ) { MIN_GRAPH_SIZE = atoi( argv[++i] ); }
		else if ( string(argv[i]) == "-fam" ) { FAM_FILE = argv[++i]; }
		else CLUSTER_FILE = argv[i];
	}

	if ( WINDOW_SIZE <= 0 ) { cerr << "ERROR: -win " << WINDOW_SIZE << ": must be greater than zero" << endl; good = false; }
	if ( MIN_CLUSTER_DENSITY <= 0 || MIN_CLUSTER_DENSITY > 1 ) { cerr << "ERROR: -density " << MIN_CLUSTER_DENSITY << ": must be greater than zero and less than or equal to one" << endl; good = false; }
	if ( MAX_DIFFERENCE_FOR_PRINT < 0 || MAX_DIFFERENCE_FOR_PRINT > 1 ) { cerr << "ERROR: -r2 " << MAX_DIFFERENCE_FOR_PRINT << ": must be greater than or equal to zero and less than or equal to one" << endl; good = false; }
	if ( MIN_GRAPH_SIZE <= 1 ) { cerr << "ERROR: -min " << MIN_GRAPH_SIZE << ": must be greater than one" << endl; good = false; }
	if ( CLUSTER_FILE == "" ) { cerr << "ERROR: No cluster output file provided" << endl; good = false; }
	
	if ( ! good ) cerr << endl;
	return good;
}

/// Comparator for edges; sorts by start position and then end position
bool compare_matches( Match * first , Match * second )
{
	if ( first->getPosition( 0 ) < second->getPosition( 0 ) ) return true;
	else if ( first->getPosition( 0 ) > second->getPosition( 0 ) ) return false;
	else if ( first->getPosition( 1 ) < second->getPosition( 1 ) ) return true;
	else if ( first->getPosition( 1 ) > second->getPosition( 1 ) ) return false;
	else return first < second;
}

int main( int argc , char * argv[] )
{
	if ( argc < 2 || !process_param( argc , argv ) )
	{
		cerr << "Usage: cat [shared segments] | " << argv[0] << " <optional parameters> [cluster output file]" << endl;
		cerr << "\t-help\t\tPrint this output" << endl;
		cerr << "\t-win [VAL]\tSliding window size in bp (default: " << WINDOW_SIZE << ")" << endl;
		cerr << "\t-fam [STRING]\tFilename for list of sample identifies in PLINK format" << endl;
		cerr << "\t-density [VAL]\tMinimum cluster density (default: "<< MIN_CLUSTER_DENSITY << ")" << endl;
		cerr << "\t-r2 [VAL]\tMaximum r^2 for which two haplotypes are considered different and printed, set to 1 to print all (default: " << MAX_DIFFERENCE_FOR_PRINT << ")" << endl;
		cerr << "\t-min [VAL]\tMinimum haplotype size (default: " << MIN_GRAPH_SIZE << ")" << endl;

		return -1;
	}

	ofstream of_cluster( (CLUSTER_FILE + ".hcl").c_str() );
	ofstream of_tped;
	of_log.open( (CLUSTER_FILE + ".log").c_str() );

        of_log << setw(95) << setfill('-') << ' ' << endl << setfill(' ');
        of_log << " Welcome to DASH, a tool for detecing shared haplotypes from" << endl;
        of_log << " pairwise Identity By Descent (IBD) data." << endl;
        of_log << endl;
        of_log << " For more details, please see the web-site at:" << endl;
        of_log << " http://www.cs.columbia.edu/~gusev/dash/" << endl;
        of_log << endl;
        of_log << " DASH was coded by Alexander Gusev and collaborators in " << endl;
        of_log << " Itsik Pe'er's Computational Biology Lab at Columbia University" << endl;
        of_log << setw(95) << setfill('-') << ' ' << endl << setfill(' ');

	if ( !of_cluster )
	{
		cerr << "ERROR! Could not open " << CLUSTER_FILE << ".hcl" << " for writing" << endl;
		return -1;
	}

	string line , str_fid[2] , str_id[2];
	unsigned long pos[2];
	
	Match * cur_match;
	size_t id[2];
	list< Match * > matches;
	list< Match * > active_matches;
	list< Match * > new_matches;

	stringstream ss;

	// load fam file
	if ( FAM_FILE != "" )
	{
		of_log << "Reading FAMily file " << FAM_FILE << endl;
		of_tped.open( (CLUSTER_FILE + ".ped").c_str() );
		ifstream fam_list( FAM_FILE.c_str() );
		while ( getline( fam_list , line ) )
		{
			ss.clear(); ss.str( line );
			ss >> str_fid[0] >> str_id[0];
			MAIN_PED.addFam( str_fid[0] + " " + str_id[0] , line );
		}
		fam_list.close();
		of_log << "Identified " << MAIN_PED.getSize() << " unique individuals in pedigree" << endl;
	}

	// load & sort all matches
	of_log << "Loading IBD shared segments" << endl;
	while( getline( cin , line ) )
	{
		ss.clear(); ss.str( line );
		ss >> str_fid[0] >> str_id[0] >> str_fid[1] >> str_id[1] >> pos[0] >> pos[1];
		id[0] = MAIN_PED.add( str_fid[0] + " " + str_id[0] );
		id[1] = MAIN_PED.add( str_fid[1] + " " + str_id[1] );
		// check if match is long enough
		if ( pos[1] - pos[0] > WINDOW_SIZE )
		{
			cur_match = new Match( id , pos );
			matches.push_back( cur_match );
		}
	}
	of_log << matches.size() << " IBD segments loaded" << endl;
	if ( matches.size() == 0 )
	{
		of_log << "WARNING: No matches were identified, make sure your match file is not empty and the window parameter is small enough to contain the matches" << endl;
		cerr << "WARNING:\tNo matches were identified" << endl << "\t\tMake sure your match file is not empty and the window parameter is small enough to contain the matches" << endl;
	}
	else
	{
		matches.sort(compare_matches);
		of_log << "Sorting complete, range: " << matches.front()->getPosition( 0 ) << " - " << matches.back()->getPosition( 1 ) << endl;
	}
        of_log << setw(95) << setfill('-') << '-' << endl << setfill(' ');
       	of_log << left << setw(25) << "Window" << setw(15) << "CL Count" << setw(15) << "Avg. Size" << setw(15) << "Max. Size" << setw(15) << "Avg. Density" << endl;

	cur_pos = 0;
	for ( list< Match * >::iterator i = matches.begin(); i != matches.end() || !active_matches.empty(); )
	{
		if ( i != matches.end() )
		{
			// catch up to this position
			while ( cur_pos < (*i)->getPosition( 0 ) + WINDOW_SIZE )
			{
				cur_pos += WINDOW_SIZE;
			}

			// remove any inactive edges
			for ( list< Match * >::iterator am = active_matches.begin(); am != active_matches.end(); )
			{
				if ( (*am)->getPosition( 1 ) < cur_pos )
				{
					MAIN_PED.disconnect( *am );
					delete *am;
					active_matches.erase( am++ );
				} else am++;
			}

			// add all matches within this window
			while ( i != matches.end() && (*i)->getPosition( 0 ) <= (cur_pos - WINDOW_SIZE) )
			{
				// add the new edge
				if ( (*i)->getPosition( 1 ) >= cur_pos )
				{
					MAIN_PED.connect( *i );
					active_matches.push_back( *i );
					new_matches.push_back( *i );
				} else 
				{
					// this match does not take up an entire window
					delete *i;
				}
				matches.erase( i++ );
			}
		}
		else
		{
			// find when the next match will be erased
			unsigned long min_end = 0;
			for ( list< Match * >::iterator am = active_matches.begin(); am != active_matches.end(); am++ )
				if ( am == active_matches.begin() || (*am)->getPosition(1) < min_end ) min_end = (*am)->getPosition(1);
			// surpass this window
			while ( min_end >= cur_pos )
			{
				cur_pos += WINDOW_SIZE;
			}

			// remove any inactive edges
			for ( list< Match * >::iterator am = active_matches.begin(); am != active_matches.end(); )
			{
				if ( (*am)->getPosition( 1 ) < cur_pos )
				{
					MAIN_PED.disconnect( *am );
					delete *am;
					active_matches.erase( am++ );
				} else am++;
			}
		}

		cerr << '\r' << (cur_pos - WINDOW_SIZE) << " .. " << cur_pos << " : " << matches.size() << " un-analyzed\t" << flush;

		// Dissolve (make visisble) any known cluster that is now too sparse
		MAIN_PED.dissolveClusters();
		// Add incident non-clustered individuals to clusters
		MAIN_PED.packClusters();
		// Identify (and make invisible) and dense subgraphs
		MAIN_PED.formClusters();
		// Print the clusters
		MAIN_PED.print( of_cluster );
		MAIN_PED.print();
	}
	if ( of_tped ) { MAIN_PED.printPed( of_tped ); of_tped.close(); }
        of_log << setw(95) << setfill('-') << '-' << endl << setfill(' ');
	cerr << endl << "Analysis Complete" << endl;
	of_log << "Analysis Complete" << endl;
	of_cluster.close();
	of_log.close();
	return 0;
}
