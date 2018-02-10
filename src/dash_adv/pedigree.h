#ifndef PEDIGREE_H
#define PEDIGREE_H

#include "main.h"
#include "match.h"
#include "individual.h"
#include "famsample.h"

#include <map>
#include <string>
#include <set>
#include <vector>
#include <iostream>
#include <list>

#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace boost;
using namespace std;

/** Maintains connections between individual IDs and Genotype pointers */
class Pedigree
{
public:
	/** Constructor for Pedigree */
	Pedigree();
	/** Destructor for Pedigree; deletes all individuals and clusters */
	~Pedigree();
	/** Add or get an individual with the parameter ID
		@param i string id of individual to find or add
		@return new Individual or found Individual if ID exists
		*/
	size_t add(string i);
	void addFam(string i , string line);
	
	void formClusters();
	void dissolveClusters();
	void packClusters();

	/** Add a new Match to the graph
		@param m pointer to match to add
		*/
	void connect( Match * m );
	/** Remove an old Match to the graph
		@param m pointer to match to add
		*/
	void disconnect( Match * m );
	/** Output information about the pedigree */
	void print();
	void print( ofstream& f );
	/** Code the clusters as SNPs for a PED file */
	void printPed( ofstream& f );
	/** Return the size of the pedigree */
	size_t getSize();
private:

	/// Maintains Edge properties
	struct match_t { typedef edge_property_tag kind; };
	struct hidden_t { typedef vertex_property_tag kind;	};
	struct idpointer_t { typedef vertex_property_tag kind; };
	struct ccomponent_t { typedef vertex_property_tag kind; };

	/// Main undirected, weighted graph
	typedef subgraph< adjacency_list<setS, setS, undirectedS,
		property<vertex_index_t , size_t,
			property<hidden_t , bool,
				property<idpointer_t , Individual *,
					property<ccomponent_t , int > > > > ,
		property<edge_index_t, size_t,
			property<edge_weight_t, unsigned long,
				property<match_t, Match * > > >
	> > Graph;

	typedef property_map < Graph, edge_weight_t >::type EdgeWeightMap;
	EdgeWeightMap all_weight_map;
	property_map< Graph , match_t >::type edge_match_map;
	
	typedef property_map< Graph , hidden_t >::type VertexHiddenMap;
	VertexHiddenMap vertex_hidden_map;

	property_map< Graph , vertex_index_t >::type vertex_index_map;
	property_map< Graph , idpointer_t >::type vertex_pointer_map;
	property_map< Graph , ccomponent_t >::type vertex_cc_map;

	/// Vertex on main graph
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	/// Struct for detirmining a visible vertex
	template <typename VertexHiddenMap> struct filter_visible_vertex {
		filter_visible_vertex() {}
		filter_visible_vertex( VertexHiddenMap hidden ) :  m_hidden(hidden) {}

		template <typename Vertex> bool operator()(const Vertex& v) const { return ! get( m_hidden , v ); }
		VertexHiddenMap m_hidden;
	};
	
	filter_visible_vertex< VertexHiddenMap > filter;
	/// Filtered graph for visible vertices
	typedef filtered_graph< Graph , keep_all , filter_visible_vertex<VertexHiddenMap> > VisibleGraph;

	/// Main graph containing all vertices
	Graph main_graph;
	/// Filtered graph containing unclustered vertices
	VisibleGraph background_graph;
	/// Mapping between string id and numeric id
	map< string , size_t > sample_map;
	/// Mapping between numeric ID and Vertex pointer
	vector< Vertex > sample_vertex;
	/// Current cluster identifier
	unsigned int hcl_num , uniq_cl_num;

	/// Mapping between string id and fam pointer
	map< string , FamSample * > fam;
	/*** Parses an individual ID into a FamID (stripping suffix) */
	string parse_fam_id( string fid );

	/// Stores and performs the min-cut
	struct MinCutDB {
			int * v;
			long * w;
			bool * a;
			bool ** nodes;
			bool * best_nodes;
			long ** g;
			size_t n;
			
			Graph& main_graph;
			EdgeWeightMap& weight_map;

			MinCutDB( Graph& mgraph , EdgeWeightMap& wmap , size_t init_n ):main_graph(mgraph),weight_map(wmap)
			{
				n = init_n;
				v = new int[n];
				w = new long[n];
				a = new bool[n];
				nodes = new bool * [n];
				g = new long * [n];
				best_nodes = new bool[n];
				for ( unsigned int i = 0 ; i < n ; i++ ) { nodes[i] = new bool[n]; g[i] = new long[n]; }
			}
			~MinCutDB()
			{
				delete[] v;
				delete[] w;
				delete[] a;
				delete[] best_nodes;
				for ( unsigned int i = 0 ; i < n ; i++ ) { delete[] nodes[i]; delete[] g[i]; }
				delete[] g;
				delete[] nodes;
			}
				
			/** Initiates the internal data structures based on the map parameter
				@param map reference to list of Vertices for this cut
				@return density of the current map
				**/
			float init( list< Vertex >& map );
			/** Calculates the density of the current map
				@param map reference to list of Vertices in this graph
				@return density of the current map
				**/
			float getDensity( list< Vertex >& map );
			void cut();
			/// Remove outlier individuals 
			void trim( list< Vertex >& map );
	};

	/// Maintains saved bit vetors of all active clusters for comparison
	map< Graph * , dynamic_bitset<> > graph_safe;

	void clusterRecurse( list< Vertex >& map , MinCutDB * mcdb );
	void trim( Graph& sgraph );
	void dissolve( Graph& sgraph );

	/** Compares the current graph membership vector to the saved one and saves new vector
		@param sgraph reference to sugraph that is to be compared
		@return percentage difference between vectors
		*/
	float compareSaved( Graph * key , Graph& sgraph );

	/** Removes a graph from the saved vectors */
	void clearSaved( Graph * key );

	string intToString(int x);
};

#endif
