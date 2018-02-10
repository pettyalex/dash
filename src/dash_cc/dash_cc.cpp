#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <list>

using namespace std;

// Global Constants
short START = 0;
short END = 1;
class Ind;

// Global Variables
struct Share
{

	Share(long p[2])
	{ 
		overlap[0] = overlap[1] = -1;
		pos[0] = p[0];
		pos[1] = p[1];
	}

	Share(long p0, long p1)
	{
		overlap[0] = overlap[1] = -1;
		pos[0] = p0;
		pos[1] = p1;
	}

	Share(long p[2] , Ind * in)
	{ 
		overlap[0] = overlap[1] = -1;
		pos[0] = p[0];
		pos[1] = p[1];
		connect ( in );
	}
	void	connect( Ind * );
	void	connect();
	void	disconnect();
	void	add ( Ind * );
	void	cut();
	void	merge( Share * );
	void	print(ostream&);
	void	clear();
	void	countPheno();
	void	extend();
	void	setID( int ); 
	// Shrink the overlap START/END markers
	void	addOverlap( long p[2] );
	// Grow the overlap START/END markers
	void	growOverlap( long p[2] );
	bool	equals( Share * );
	Share * split( long );
	
	int csid;
	set<Ind*> i;
	long	pos[2];

	// Track the max/min start/end of any overlapping segments
	long	overlap[2];
	float *	pheno;
};

struct CompareShares
{
	bool operator()(const Share * lhs, const Share * rhs) const { return lhs->pos[0] < rhs->pos[0]; }
};

struct Sample
{
	Sample( string l ) { line = l; }
	void addShare( int id ) { if ( cs.back().first == id ) cs.back().second = true; else cs.push_back( make_pair(id,false) ); }
	list< pair<int,bool> > cs;
	string line;
};

class Ind
{
public:
	Ind( Sample * p , string i ) { parent = p; id = i; }
	void	setParent( Sample * p ) { parent = p; }
	void	addShareID( int id ) { if ( parent != NULL ) parent->addShare( id ); }
	void	addShare(Share * sha){ if ( !s.insert(sha).second ) cerr << "WARNING: New cluster already present: " << sha->pos[0] << ' ' << sha->pos[1] << endl; }
	void	delShare(Share * sha){ s.erase(sha); }
	void	cut( long pos[2] , list<Share*>& ret );
	string	getID() { return id; }
	set< Share* , CompareShares > s;
private:
	Sample * parent;
	string id;
};

list<Share*> all_shares;

string parse_id( string id )
{
	return id.substr( 0, id.length()-2 );
}

void deleteShare ( Share * del )
{
	all_shares.remove( del );
	delete del;
}

void addShare( Share * add )
{
	all_shares.push_back( add );
}

void Share::addOverlap( long p[2] )
{
	if ( overlap[START] == -1 || p[START] > overlap[START] ) overlap[START] = p[START];
	if ( overlap[END] == -1 || p[END] < overlap[END] ) overlap[END] = p[END];
}
void Share::growOverlap( long p[2] )
{
	if ( overlap[START] == -1 || p[START] < overlap[START] ) overlap[START] = p[START];
	if ( overlap[END] == -1 || p[END] > overlap[END] ) overlap[END] = p[END];
}

void Share::setID( int id )
{
	csid = id;
	for ( set< Ind* >::iterator it = i.begin() ; it != i.end() ; it++ ) (*it)->addShareID( id );
}

void Share::extend()
{
	// find this share in the first individual
	set< Share* , CompareShares >::iterator	find = (*i.begin())->s.find( this );
	// move to the next share
	find++;

	if ( find != (*i.begin())->s.end() && 
		 this->pos[END] == (*find)->pos[START] && 
		 this->equals( *find ) )
	{
		// save the second to be deleted
		Share * del = *find;

		// extend the first to the end of the second
		this->pos[END] = (*find)->pos[END];
		this->growOverlap( (*find)->overlap );

		// delete the second
		del->disconnect();
		deleteShare( del );
	}
}

void mergeShares ( list<Share*> sha[2] , Ind * ind[2] , long pos[2] )
{
	list<Share*>::iterator sha_i[2];
	sha_i[0] = sha[0].begin();
	sha_i[1] = sha[1].begin();
	int minor;
	// the farthest position we have examined
	long tail = pos[0];
	Share * new_share;

	if ( sha_i[0] == sha[0].end() && sha_i[1] == sha[1].end() )
	{
		new_share = new Share ( pos );
		new_share->connect(ind[0]);
		new_share->connect(ind[1]);
		new_share->addOverlap( pos );
		addShare( new_share ); 
	}
	else
	{
		// There is a gap at the start
		long min_start;
		if ( sha_i[0] != sha[0].end() )
		{
			min_start = (*sha_i[0])->pos[START];
			if ( sha_i[1] != sha[1].end() && (*sha_i[1])->pos[START] < min_start ) min_start = (*sha_i[1])->pos[START];
		} else min_start = (*sha_i[1])->pos[START];
		// Create a new share for the gap
		if ( tail < min_start )
		{
			new_share = new Share( tail , min_start );
			new_share->connect( ind[0] );
			new_share->connect( ind[1] );
			new_share->addOverlap( pos );
			addShare( new_share );
			tail = min_start;
		}

		while ( sha_i[0] != sha[0].end() && sha_i[1] != sha[1].end() )
		{
			// Equal start and end
			if ( (*sha_i[0])->pos[START] == (*sha_i[1])->pos[START] && (*sha_i[0])->pos[END] == (*sha_i[1])->pos[END] )
			{
				tail = (*sha_i[0])->pos[END];
				
				// merge across both
				if ( *sha_i[0] != *sha_i[1] )
				{
					if ( (*sha_i[0])->i.size() < (*sha_i[1])->i.size() ) minor = 0; else minor = 1;

					// disconnect individuals from the minor
					(*sha_i[ minor ])->disconnect();

					// move them into the major
					(*sha_i[ 1 - minor ])->merge ( *sha_i[ minor ] );
					(*sha_i[1-minor])->addOverlap( pos );

					// delete the minor
					deleteShare( *sha_i[ minor ] );
				}

				// increment both
				sha_i[0]++; sha_i[1]++;
			}
			// Equal starts
			else if ( (*sha_i[0])->pos[START] == (*sha_i[1])->pos[START] )
			{

				if ( (*sha_i[0] )->pos[END] < (*sha_i[1] )->pos[END] ) minor = 0; else minor = 1;
				// we cover until the minor ends
				tail = (*sha_i[minor] )->pos[END];

				// shift the major share to start from the minor end
				(*sha_i[1-minor])->disconnect();
				(*sha_i[1-minor])->pos[START] = (*sha_i[minor])->pos[END];
				(*sha_i[1-minor])->connect();

				// merge major individuals into the minor
				(*sha_i[minor])->merge( *sha_i[1-minor] );
				(*sha_i[minor])->addOverlap( pos );
				// increment the minor
				sha[minor].erase(sha_i[minor]++);
			}
			// Equal ends - SHOULD NOT HAPPEN
/*			else if ( (*sha_i[0])->pos[END] == (*sha_i[1])->pos[END] )
			{
				cerr << "Equal Ends" << endl;
				if ( (*sha_i[0] )->pos[START] < (*sha_i[1] )->pos[START] ) minor = 0; else minor = 1;

				// we cover until the major starts
				tail = (*sha_i[1-minor] )->pos[START];

				// shift the minor share to end from the major start
				(*sha_i[minor])->pos[END] = (*sha_i[1-minor])->pos[START];

				// merge major individuals into the minor
				(*sha_i[minor])->merge( *sha_i[1-minor] );
				(*sha_i[minor])->addOverlap( pos );
				// increment the minor
				sha[minor].erase(sha_i[minor]++);
			}
*/
			// Some overlap
			else if ( (*sha_i[0])->pos[START] < (*sha_i[1])->pos[END] && (*sha_i[0])->pos[END] > (*sha_i[1])->pos[START] )
			{
				if ( (*sha_i[0])->pos[START] < (*sha_i[1])->pos[START] ) minor = 0; else minor = 1;
				// we cover until the major starts
				tail = (*sha_i[ 1 - minor ] )->pos[START];

				// create a new share over the non-overlap
				new_share = new Share( (*sha_i[minor])->pos[START], (*sha_i[1-minor])->pos[START] );
				// disconnect & move the minor share to the end of the non-overlap
				(*sha_i[minor])->disconnect();
				(*sha_i[minor])->pos[START] = (*sha_i[1-minor])->pos[START];
				(*sha_i[minor])->connect();
				// put all of the minor individuals into the non-overlap
				new_share->merge( *sha_i[minor] );
				// also put the major individual into the non-overlap
				new_share->connect( ind[1-minor] );
				// save
				new_share->addOverlap( pos );
				addShare(new_share);
			}
			// No overlap
			else
			{
				if ( (*sha_i[0])->pos[START] < (*sha_i[1])->pos[START] ) minor = 0; else minor = 1;
				// we cover until the minor ends
				tail = (*sha_i[ minor ] )->pos[END];

				// Cover the minor share
				(*sha_i[minor])->connect( ind[1-minor] );
				(*sha_i[minor])->addOverlap( pos );

				// Cover the gap between minor and major
				if ( (*sha_i[minor])->pos[END] < (*sha_i[1-minor])->pos[START] )
				{
					new_share = new Share( (*sha_i[minor])->pos[END] , (*sha_i[1-minor])->pos[START] );
					new_share->connect( ind[1-minor] );
					new_share->addOverlap( pos );
					addShare( new_share );

					// Put the gap as the next "major" element to be processed
					sha_i[1-minor] = sha[1-minor].insert( sha_i[1-minor] , new_share );
				}
				// Move to the next "minor" element
				sha_i[minor]++;
			}
		}
		// One individual still has shares
		if ( sha_i[0] != sha[0].end() || sha_i[1] != sha[1].end() )
		{
			if ( sha_i[0] != sha[0].end() ) minor = 0; else minor = 1;
			while ( sha_i[minor] != sha[minor].end() )
			{
				// connect the remaining matches
				tail = (*sha_i[minor])->pos[END];
				(*sha_i[minor])->connect( ind[ 1-minor ] );
				(*sha_i[minor])->addOverlap( pos );
				sha_i[minor]++;
			}
		}

		// There is a gap at the end
		if ( tail < pos[END] )
		{
			new_share = new Share( tail , pos[END] );
			new_share->connect( ind[0] );
			new_share->connect( ind[1] );
			new_share->addOverlap( pos );
			addShare( new_share );
		}
	}
}

bool Share::equals( Share * cmp )
{
	bool eq = false;
	if ( i.size() == cmp->i.size() )
	{
		eq = true;
		set< Ind *>::iterator iter[2];
		for ( iter[0] = i.begin(), iter[1] = cmp->i.begin();
			  iter[0] != i.end();
			  iter[0]++, iter[1]++ )
		{
			if ( *iter[0] != *iter[1] ) { eq = false; break; }
		}
	}
	return eq;
}

void Share::merge( Share * in )
{
	for ( set< Ind * >::iterator iter = in->i.begin() ; iter != in->i.end(); iter++ ) connect( *iter );
}

Share * Share::split( long pos )
{
	Share * right = new Share(*this);
	this->pos[END] = right->pos[START] = pos;
	addShare( right );
	return right;
}

void Ind::cut( long pos[2] , list<Share*>& ret )
{
	Share * buffer[2]; buffer[0] = buffer[1] = NULL;
	for ( set< Share* , CompareShares >::iterator i = s.begin(); i != s.end(); i++ )
	{
		// End if we've passed the region
		if ( (*i)->pos[START] > pos[END] ) break;
		// Cut if we're overlapping the region
		else if ( pos[START] < (*i)->pos[END] && pos[END] > (*i)->pos[START] )
		{
			Share * cur_share = (*i);
			// Detirmine the overlap start
			if ( pos[START] > cur_share->pos[START] )
			{
				// split the current share at the region start
				// keep the right leftover
				buffer[0] = cur_share = cur_share->split( pos[START] );
			}
			
			// Detirmine the overlap end
			if ( pos[END] < cur_share->pos[END] )
			{
				// split the current share at the region end
				// keep the left leftover
				buffer[1] = cur_share->split( pos[END] );
			}

			// Mark
			ret.push_back( cur_share );
		}
	}
	if ( buffer[0] != NULL ) buffer[0]->connect();
	if ( buffer[1] != NULL ) buffer[1]->connect();
}

void Share::connect( Ind * in )
{
	if ( i.insert( in ).second  )
	{
		in->addShare( this );
	}
}

void Share::connect()
{
	for( set<Ind*>::iterator iter = i.begin(); iter != i.end(); iter++ ) (*iter)->addShare( this );
}

void Share::disconnect()
{
	for( set<Ind*>::iterator iter = i.begin(); iter != i.end(); iter++ ) (*iter)->delShare( this );
}

void Share::add( Ind * in )
{
	i.insert( in );
}

void Share::cut()
{
	for(set<Ind*>::iterator iter = i.begin(); iter != i.end(); iter++) (*iter)->delShare(this);
}

void Share::print(ostream& o)
{
	o << "cs" << csid << '\t' << pos[0] << '\t' << pos[1] << '\t' << overlap[0] << '\t' << overlap[1];
	for(set<Ind*>::iterator iter = i.begin(); iter != i.end(); iter++) o << '\t' << (*iter)->getID();
}

int main (int argc, char* argv[]) 
{
	list< Sample * > fam_list;
	map< string , Sample* > fam;
	map< string , Sample* >::iterator fam_search;

	map< string , Ind* > haps;
	map< string , Ind* >::iterator haps_search[2];

	Ind * id_pointer[2];
	string in , line , id[2] , f_id[2] , pid;
	stringstream ss;
	Ind * new_ind;
	Sample * new_sample;

	if ( argc < 3 )
	{
		cout << "Usage: cat IBD | " << argv[0] << " [fam file] [output prefix]" << endl;
		cout << "IBD format: [fid1] [id1] [fid2] [id2] [match start] [match end]" << endl;
		return 0;
	}

	ifstream f_fam( argv[1] );
	if ( !f_fam ) { cerr << "ERROR: Couldn't open fam file" << endl; return 0; }

	while(getline(f_fam,line))
	{
		ss.clear(); ss.str(line);
		ss >> f_id[0] >> id[0];
		new_sample = new Sample( line );
		fam.insert(make_pair(f_id[0] + " " + id[0] , new_sample));
		fam_list.push_back( new_sample );
	}
	f_fam.close();

	// Process Matches
	long pos[2];
	list<Share*> cur_shares[2];
	bool skip;
	long ctr = 0;

	while(getline(cin,line))
	{
		if ( ctr++ % 100000 == 0 ) { cerr << '\r' << ctr-1 << " segments" << flush; }
		skip = false;
		ss.clear(); ss.str(line);
		ss >> f_id[0] >> id[0];
		ss >> f_id[1] >> id[1];
		ss >> pos[START];
		ss >> pos[END];
		
		for( int i=0; i<2; i++ )
		{
			pid = f_id[i] + " " + id[i];
			if( (haps_search[i] = haps.find( pid )) == haps.end() )
			{
				if( (fam_search = fam.find( parse_id(pid) )) == fam.end() )
				{
					cerr << "WARNING: Individual " << pid << " was not found in fam file, will be skipped" << endl;
					skip = true;
					break;
				} else
				{
					haps_search[i] = haps.insert(make_pair( pid , new Ind( fam_search->second , pid ) )).first;
				}
			}
			
			// get all shares from this individual and this location
			cur_shares[i].clear();
			id_pointer[i] = haps_search[i]->second;
			haps_search[i]->second->cut( pos , cur_shares[i] );
		}

		if ( !skip ) mergeShares( cur_shares , id_pointer , pos );
	}
	cerr << '\r' << "Completed " << ctr << " segments" << endl;

	// merge consecutive & identical shares
	cerr << all_shares.size() << " initial clusters found" << endl;
	for ( list<Share*>::iterator i = all_shares.begin(); i != all_shares.end(); i++ ) (*i)->extend();

	ofstream out_cc( (string(argv[2]) + ".clst").c_str() );

	int csid = 0;
	int csctr;

	// print all the shares
	for( list<Share*>::iterator i = all_shares.begin(); i != all_shares.end(); )
	{
		if ( (*i)->i.size() <= 2 )
		{
			delete (*i);
			all_shares.erase( i++ );
		} else
		{
			(*i)->setID( csid );
			(*i)->print(out_cc); out_cc << endl;
			csid++;
			i++;
		}
	}
	out_cc.close();

	// print ped file
	ofstream out_ped( (string(argv[2]) + ".ped").c_str() );

	cerr << csid << " total clusters" << endl;
	for( list<Sample*>::iterator i = fam_list.begin() ; i != fam_list.end() ; i++ )
	{
		out_ped << (*i)->line;
		list< pair<int,bool> >::iterator cs_i = (*i)->cs.begin();
		for ( csctr = 0 ; csctr < csid ; csctr++ )
		{
			if ( cs_i == (*i)->cs.end() || cs_i->first > csctr )
			{
				out_ped << " 1 1";
			} else if ( cs_i->first == csctr )
			{
				if ( cs_i->second ) out_ped << " 2 2"; else out_ped << " 1 2";
				cs_i++;
			} else
			{
				while ( cs_i != (*i)->cs.end() && cs_i->first < csctr ) cs_i++;
			}
		}
		out_ped << endl;
	}
	out_ped.close();	
	return 1;
}

