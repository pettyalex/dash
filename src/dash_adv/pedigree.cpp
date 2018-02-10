#include "pedigree.h"

Pedigree::Pedigree():filter( get(hidden_t() , main_graph) ) , background_graph(main_graph , keep_all() , filter )
{
	all_weight_map = get(edge_weight , main_graph );
	edge_match_map = get(match_t(), main_graph );
	
    vertex_index_map = get(vertex_index, main_graph);
	vertex_hidden_map = get(hidden_t() , main_graph);
	vertex_pointer_map = get(idpointer_t() , main_graph);
	vertex_cc_map = get( ccomponent_t() , main_graph);
	hcl_num = uniq_cl_num = 0;
}

Pedigree::~Pedigree()
{
	for ( size_t i = 0 ; i < sample_vertex.size() ; i++ )
		delete vertex_pointer_map[ sample_vertex[i] ];
}

string Pedigree::parse_fam_id( string fid )
{
	string sub_id = fid.substr( 0, fid.length()-2 );
	string suff = fid.substr( fid.length()-2 );
	if ( suff != ".0" && suff != ".1" )
	{
		cerr << left << suff << endl;
		cerr << setw(10) << "WARNING:" << "Sample identifier \"" << fid << "\" does not have a .0/.1 suffix" << endl;
		cerr << setw(10) << " " << "When providing a FAM file, individual IBD segment ID's must have suffix designating" << endl;
		cerr << setw(10) << " " << "which haploid copy is involved in the current segment" << endl;
	}
	return sub_id;
}

void Pedigree::addFam( string i , string line )
{
        if ( fam.find( i ) == fam.end() )
        {
			FamSample * new_sample = new FamSample( line );
			fam.insert( make_pair( i , new_sample ) );
        } else of_log << "WARNING: Duplicate individual " << i << " in FAM file" << endl;
}

size_t Pedigree::add(string i)
{
        map< string , size_t >::iterator it;
        if ( ( it = sample_map.find( i )) == sample_map.end() )
        {
			size_t id = sample_vertex.size();
			sample_map.insert( make_pair( i , id ));
			sample_vertex.push_back( add_vertex(main_graph) );

			// Try to find a FamSample for this individual
			FamSample * fs = NULL;
			if ( fam.size() != 0 )
			{
				map< string , FamSample * >::iterator fam_it = fam.find( parse_fam_id(i) );
				if ( fam_it != fam.end() ) fs = fam_it->second;
				else { cerr << "WARNING: " << i << " could not be found in the FAM file, it will be skipped in the PED file" << endl; }
			}

			vertex_pointer_map[ sample_vertex.back() ] = new Individual( id , i , fs );
			vertex_hidden_map[ sample_vertex.back() ] = false;
			vertex_index_map[ sample_vertex.back() ] = id;
			return id;
        } else return it->second;
}

void Pedigree::dissolveClusters()
{
	size_t n;
	bool destroy;
	for ( list< Graph* >::iterator ci = main_graph.m_children.begin() ; ci != main_graph.m_children.end() ; )
	{
		n = num_vertices( **ci );
		destroy = n < MIN_GRAPH_SIZE || (float) (2.0 * num_edges(**ci) / (n * (n-1)) < MIN_CLUSTER_DENSITY );
		// if dense presently, try removing outliers and re-check density
		if ( ! destroy )
		{
			trim( **ci );
			n = num_vertices( **ci );
			destroy = n < MIN_GRAPH_SIZE || (float) (2.0 * num_edges(**ci) / (n * (n-1)) < MIN_CLUSTER_DENSITY );
		}

		if ( destroy )
		{
			dissolve( **ci );
			clearSaved( *ci );
			main_graph.m_children.erase( ci++ );
		} else ci++;
	}
}

void Pedigree::dissolve( Graph& sgraph )
{
	Graph::vertex_iterator vi , vi_end;
	for( tie(vi, vi_end) = vertices( sgraph ); vi != vi_end; ++vi )
		vertex_hidden_map[ sgraph.local_to_global( *vi ) ] = false;
}

void Pedigree::trim( Graph& sgraph )
{
	Graph::vertex_iterator vi , vi_end, next;
	unsigned int min = int((num_vertices( sgraph )-1) * MIN_CLUSTER_DENSITY);
	for( tie(vi, vi_end) = vertices( sgraph ); vi != vi_end; ++vi )
	{
		if ( out_degree( *vi , sgraph ) < min )
		{
			vertex_hidden_map[ sgraph.local_to_global( *vi ) ] = false;
		}
	}
	
	tie(vi, vi_end) = vertices( sgraph );
	for ( next = vi ; vi != vi_end ; vi = next )
	{
		++next;
		if ( ! vertex_hidden_map[ sgraph.local_to_global( *vi ) ] )
		{
			clear_vertex( *vi , sgraph );
			remove_vertex( *vi , sgraph );
		}		
	}
}

void Pedigree::packClusters()
{

	Graph::children_iterator ci, ci_end;
	Graph::edge_iterator ei , ei_end;
	Graph::vertex_iterator vi , vi_end;
	VisibleGraph::adjacency_iterator ai , ai_end;

	map< Vertex , size_t > neighbors;
	map< Vertex , size_t >::iterator n_query;

	unsigned int min;

	// For each cluster:
	// TODO: Order in decreasing cluster size
	for (tie(ci, ci_end) = main_graph.children(); ci != ci_end; ci++ )
	{
		neighbors.clear();
		min = int((num_vertices( *ci )) * MIN_CLUSTER_DENSITY);

		// Construct a list of adjacent visible vertices
		for ( tie( vi , vi_end ) = vertices( *ci ) ; vi != vi_end ; vi++ )
		{
			for ( tie( ai , ai_end ) = adjacent_vertices( ci->local_to_global( *vi ) , background_graph ) ; ai != ai_end ; ai++ )
			{
				n_query = neighbors.find( *ai );
				if ( n_query == neighbors.end() ) neighbors.insert( make_pair( *ai , 1 ) );
				else n_query->second++;
			}
		}
		// Add those nodes that are incident to MIN_CLUSTER_DENSITY nodes
		for ( n_query = neighbors.begin() ; n_query != neighbors.end() ; n_query++ )
		{
			if ( n_query->second >= min )
			{
				vertex_hidden_map[ n_query->first ] = true;
				add_vertex( n_query->first , *ci );
			}
		}
	}
}

void Pedigree::formClusters()
{	
	// Identify the connected components
	int nc = connected_components(background_graph , get(ccomponent_t(), main_graph) );
	if ( nc == 0 ) { return; }

	vector< list< Vertex > > component( nc );
	Graph::vertex_iterator   vi, vi_end;
	for(tie(vi, vi_end) = vertices(main_graph) ; vi != vi_end; ++vi)
	{
		if ( vertex_cc_map[ *vi ] >= 0 )
		{
			component[ vertex_cc_map[*vi] ].push_back( *vi );
			vertex_cc_map[*vi] = -1;
		}
	}

	MinCutDB * mcdb;
	int ncc = 0;
	for ( int i = 0 ; i < nc ; i++ )
	{
		if ( component[i].size() >= MIN_GRAPH_SIZE )
		{
			ncc++;
#if VERBOSE > 0
			cout << "CC:\t";
			for ( list< Vertex >::iterator cc = component[i].begin() ; cc != component[i].end() ; cc++ ) cout << main_graph[ *cc ].id->getID() << ' ';
			cout << endl;
#endif

			mcdb = new MinCutDB( main_graph , all_weight_map , component[i].size() );
			clusterRecurse( component[i] , mcdb );
			delete mcdb;
		}
	}
}

float Pedigree::MinCutDB::getDensity( list< Vertex >& map )
{
	n = map.size();
	unsigned int edge_ctr = 0;
	list< Vertex >::iterator it_i, it_j;

	for( it_i = map.begin(); it_i != map.end(); it_i++ )
	{
		for ( it_j = map.begin() ; it_j != map.end(); it_j++ )
		{
			if ( edge( *it_i , *it_j , main_graph ).second )
				edge_ctr++;
		}
	}
	return (float) edge_ctr / (n * (n-1));
}

float Pedigree::MinCutDB::init( list< Vertex >& map )
{
	n = map.size();
	pair< Edge , bool > p;
	unsigned int edge_ctr = 0;
	size_t i , j;
	list< Vertex >::iterator it_i, it_j;

	for( i = 0 , it_i = map.begin(); it_i != map.end(); i++,it_i++ )
	{
		v[i] = (int) i;
		for ( j = 0, it_j = map.begin() ; it_j != map.end(); j++,it_j++ )
		{
			nodes[i][j] = false;
			p = edge( *it_i , *it_j , main_graph );
			if ( p.second ) 
			{
				g[i][j] = weight_map[ p.first ];
				edge_ctr++;
			}
			else g[i][j] = 0;
		}
		nodes[i][i] = true;
	}
	return (float) edge_ctr / (n * (n-1));
}

void Pedigree::MinCutDB::trim( list< Vertex >& map )
{
	list< Vertex >::iterator it_i;
	unsigned int i , j , ctr , min;
	min = int((n-1) * MIN_CLUSTER_DENSITY);
	for( i = 0 , it_i = map.begin(); it_i != map.end(); i++)
	{
		ctr = 0;
		for ( j = 0 ; j < n ; j++ ) if ( i!=j && g[i][j] != 0 ) ctr++;
		if ( ctr < min )
			map.erase( it_i++ );
		else it_i++;
	}
}

void Pedigree::MinCutDB::cut()
{
	size_t ctr = n;
	long best = -1;
    while( ctr > 1 )
    {
        // initialize the set A and vertex weights
        a[v[0]] = true;
        for( size_t i = 1; i < ctr; i++ )
        {
            a[v[i]] = false;
            w[i] = g[v[0]][v[i]];
        }

        // add the other vertices
        int prev = v[0];
        for( size_t i = 1; i < ctr; i++ )
        {
            // find the most tightly connected non-A vertex
            int zj = -1;
            for( size_t j = 1; j < ctr; j++ )
                if( !a[v[j]] && ( zj < 0 || w[j] > w[zj] ) )
                    zj = (int) j;

            // add it to A
			a[v[zj]] = true;

            // last vertex?
            if( i == ctr - 1 )
            {
        // remember the cut weight
		if ( best == -1 || best > w[zj] )
		{
			best = w[zj];
			for ( size_t i = 0 ; i < n ; i++ ) best_nodes[i] = nodes[zj][i];
		}

                // merge prev and v[zj]
                for( size_t j = 0; j < ctr; j++ )
                    g[v[j]][prev] = g[prev][v[j]] += g[v[zj]][v[j]];
                v[zj] = v[--ctr];
		for ( size_t i = 0 ; i < n ; i++ ) { if ( nodes[prev][i] || nodes[zj][i] ) nodes[prev][i] = nodes[zj][i] = true; }
                break;
            }
            prev = v[zj];

            // update the weights of its neighbours
            for( size_t j = 1; j < ctr; j++ )
		if( !a[v[j]] ) w[j] += g[v[zj]][v[j]];
        }
    }
}

string Pedigree::intToString(int x)
{
	std::ostringstream o;
	o << x;
	return o.str();
}

void Pedigree::clusterRecurse( list< Vertex >& map , MinCutDB * mcdb )
{
	if ( map.size() < MIN_GRAPH_SIZE ) return;
	float d = mcdb->init( map );

	if ( d < MIN_CLUSTER_DENSITY )
	{
		mcdb->cut();
		list< Vertex > graph_cut[2];
		size_t ctr = 0;
		for ( list< Vertex >::iterator i = map.begin() ; i != map.end () ; i++, ctr++ )
		{
			if ( mcdb->best_nodes[ctr] ) graph_cut[0].push_back( *i ); else graph_cut[1].push_back( *i );
		}
		for ( int i = 0 ; i < 2 ; i++ ) clusterRecurse( graph_cut[i] , mcdb );
	} else { 
		bool pass = false;
		if ( d == 1 ) pass = true;
		else 
		{
			mcdb->trim( map );
			if ( map.size() >= MIN_GRAPH_SIZE && mcdb->getDensity( map ) >= MIN_CLUSTER_DENSITY ) pass = true;
		}

		if ( pass )
		{
			for ( list< Vertex >::iterator i = map.begin() ; i != map.end(); i++ )
				vertex_hidden_map[ *i ] = true;
			Graph& sg = main_graph.create_subgraph( map.begin() , map.end() );
		}
	}
}

float Pedigree::compareSaved( Graph * key , Graph& sg )
{
	size_t size = sample_vertex.size();
	dynamic_bitset<> bits( size );
	float dif = 0;
	Graph::vertex_iterator vi , vi_end;

	for ( tie(vi, vi_end) = vertices( sg ) ; vi != vi_end ; vi++ )
	{
		bits.set( vertex_pointer_map[ sg.local_to_global( *vi ) ]->getNumericID() );
	}
	
	map< Graph * , dynamic_bitset<> >::iterator it = graph_safe.find( key );
	if ( it != graph_safe.end() )
	{

		float x11 = ( (float) (bits & it->second).count() ) / size;
		float D = x11 - ( (float) bits.count() / size ) * ( (float) it->second.count() / size );

		/*
		if ( D > 0 )
		{
			float p1q2 = ( (float) bits.count() / size ) * ( (float) ( ~it->second ).count() / size );
			float p2q1 = ( (float) it->second.count() / size ) * ( (float) ( ~bits ).count() / size );
			if ( p1q2 < p2q1 ) D /= p1q2; else D /= p2q1;
		} else if ( D < 0 )
		{
			float p1q1 = -1.0f * ( (float) bits.count() / size ) * ( (float) (it->second).count() / size );
			float p2q2 = -1.0f * ( (float) ( ~bits ).count() / size ) * ( (float) (~it->second).count() / size );
			if ( p1q1 > p2q2 ) D /= p1q1; else D /= p2q2;
		}
		cout << "D-prime: " << D << endl;
		*/

		// Calculate the denominator
		float p = (float) it->second.count() / size;
		p *= (float) bits.count() / size;
		p *= 1 - ((float) it->second.count() / size);
		p *= 1 - ((float) bits.count() / size);
		if ( p != 0 ) dif = D * D / p; else dif = 0;

		// Save the new vector
		if ( dif != 1 ) it->second = bits;

	} else
	{
		graph_safe.insert( make_pair( key , bits ) );
	}

	return dif;
}

void Pedigree::clearSaved( Graph * key )
{
	if ( graph_safe.erase( key ) == 0 ) cout << "WARNING: Cleared haplotype has not been previously saved\t" << key << endl; 
}

void Pedigree::printPed( ofstream& f )
{
        for( map<string,FamSample*>::iterator i = fam.begin() ; i != fam.end() ; i++ )
        {
                f << i->second->line;

                list< pair<int,bool> >::iterator cs_i = i->second->cs.begin();
               	for ( unsigned int csctr = 0 ; csctr < hcl_num ; csctr++ )
                {
                       	if ( cs_i == i->second->cs.end() || cs_i->first > csctr )
                        {
                               	f << " 1 1";
                       	} else if ( cs_i->first == csctr )
                        {
                               	if ( cs_i->second ) f << " 2 2"; else f << " 1 2";
                                cs_i++;
                        } else
                       	{
                               	while ( cs_i != i->second->cs.end() && cs_i->first < csctr ) cs_i++;
                       	}
                }
                f << endl;
       	}
}

void Pedigree::print( ofstream& f )
{
	Graph::children_iterator ci, ci_end;
	Graph::edge_iterator ei , ei_end;
	Graph::vertex_iterator vi , vi_end;
	Match * m;

	unsigned long hap_start , hap_end;
	for (tie(ci, ci_end) = main_graph.children(); ci != ci_end; ci++ )
	{
		// test if this graph should be saved
		if ( compareSaved( *(ci.base()) , (*ci) ) > MAX_DIFFERENCE_FOR_PRINT ) continue;

		// get minimum haplotype boundary
		hap_start = hap_end = 0;
		for ( tie(ei, ei_end) = edges( *ci ) ; ei != ei_end; ei++ )
		{
			m = edge_match_map[ (*ci).local_to_global(*ei) ];
			if ( hap_start == 0 || m->getPosition(0) > hap_start ) hap_start = m->getPosition(0);
			if ( hap_end == 0 || m->getPosition(1) < hap_end ) hap_end = m->getPosition(1);
		}
		f << "cl" << hcl_num << '\t' << hap_start << '\t' << hap_end;
		for ( tie(vi, vi_end) = vertices( *ci ) ; vi != vi_end ; vi++ )
		{
			vertex_pointer_map[ (*ci).local_to_global(*vi) ]->addCluster( hcl_num );
			f << '\t' << vertex_pointer_map[ (*ci).local_to_global(*vi) ]->getID();
		}
		hcl_num++;
		f << endl;
	}
}

void Pedigree::print()
{
	Graph::children_iterator ci, ci_end;
	int num = 0;
	int size , max_size = 0 , avg_size = 0;
	double den , avg_den = 0;
	for (tie(ci, ci_end) = main_graph.children(); ci != ci_end; ++ci)
	{
			num++;
			size = (int) num_vertices( *ci );
			avg_size += size;

			den = ( 2.0f * num_edges( *ci ) ) / ( size * (size - 1) );
			avg_den += den;
			if ( size > max_size ) max_size = size;

			if ( den < MIN_CLUSTER_DENSITY - 0.01 ) cout << "WARNING: " << den << '\t' << num_edges( *ci ) << " edges\t" << size << " vertices" << endl;
	}

	if ( num > 0 ) of_log << setw(25) << cur_pos << setw(15) << num << setw(15) << avg_size / num << setw(15) << max_size << setw(15) << avg_den / num << endl;
}

size_t Pedigree::getSize()
{
	return fam.size();
}

void Pedigree::connect( Match * m )
{
	Edge e; bool inserted;
	tie( e , inserted ) = add_edge( sample_vertex[m->getSampleID(0)] , sample_vertex[m->getSampleID(1)] , main_graph );
	if ( inserted )
	{
		all_weight_map[e] = m->getLength();
		edge_match_map[e] = m;
	}
}

void Pedigree::disconnect( Match * m )
{
	remove_edge( sample_vertex[m->getSampleID(0)] , sample_vertex[m->getSampleID(1)] , main_graph );
}
