//
// Tests the percolation identification algorithm for systems with periodic boundaries of Livraghi et al. (2021).
//


//
// Includes.
//

// Standard includes.
#include <iostream>
#include <array>
#include <set>
#include <map>
#include <cassert>

// The class being tested.
#include "percPerBC.h"


//
// Main.
//
int main( int argc, char **argv )
{
	//
	// Built-in debug methods for the auxiliary methods for spatial dimensions 1, 2 and 3.
	//
	percPerBC<1> test1D;
	percPerBC<2> test2D;
	percPerBC<3> test3D;

	// Test the spanning dimension method, verbosely.
	std::cout << "Running standard tests for spanningDimension() for spatial dimensions 1, 2 and 3." << std::endl;
	test1D.debugSpanningDimension(true);
	test2D.debugSpanningDimension(true);
	test3D.debugSpanningDimension(true);

	// Test the independence method, again verbosely.
	std::cout << "\nRunning standard tests for isIndependent() for spatial dimensions 1, 2 and 3." << std::endl;
	test1D.debugIsIndependent(true);
	test2D.debugIsIndependent(true);
	test3D.debugIsIndependent(true);	// TEMP: Should be true. Also uncomment all other checks.

	//
	// 1D test. This is the same as Fig. 4 of Livraghi et al. (2021) after the bridging links are added.
	//
	std::set<int> nodes1D{0,1,2,3,4,5,6};		// Integer nodes - here, consecutive integers starting from zero.

	// Initial set: Does not percolate.
	std::map< std::array<int,2>, std::array<int,1> > edges1D {
		{ {2,5}, {0} },
		{ {5,6}, {0} },
		{ {6,4}, {0} },
		{ {4,0}, {0} },
		{ {0,1}, {0} },
		{ {0,3}, {0} },
		{ {0,2}, {0} }
	};

	int percDim1D = test1D.percolationDimension( nodes1D, edges1D );

	std::cout << "\nPercolation dimension for 1D example prior to adding 'bridging' links: " << percDim1D << " (should be 0)." << std::endl;
	assert( percDim1D==0 );
	
	// Add final two links such that the cluster perocolates.
	edges1D.emplace( std::make_pair( std::array<int,2>({2,4}), std::array<int,1>({-1}) ) );	// Same as { {4,2}, { 1} }.
	edges1D.emplace( std::make_pair( std::array<int,2>({3,1}), std::array<int,1>({ 1}) ) );	// Same as { {1,3}, {-1} }.

	percDim1D = test1D.percolationDimension( nodes1D, edges1D );

	std::cout << "Percolation dimension for 1D example after adding 'bridging' links: " << percDim1D << " (should be 1)." << std::endl;
	assert( percDim1D==1 );


	//
	// 2D Example. Network arranged as (with all links):
	//
	//    |    5              |  
	// - - - - - - - - - - - - - -
	//    |    |              |  
	//    |    0       2      |  
	//    |      \   /   \    |  
	//  3 | - 6 -  1      3 - | 6 
	//    |        |          |  
	//    |    5 - 4   7 - 8  |  
	//    |    |              | 
	// - - - - - - - - - - - - - -
	//    |    0              |  
	//
	std::set<int> nodes2D{0,1,2,3,4,5,6,7,8};

	// Initial edge set that does not have bridging crosslinks, so should not percolate at all.
	std::map< std::array<int,2>, std::array<int,2> > edges2D {
		{ {0,1}, {0,0} },
		{ {1,2}, {0,0} },
		{ {2,3}, {0,0} },
		{ {1,4}, {0,0} },
		{ {4,5}, {0,0} },
		{ {7,8}, {0,0} },
		{ {1,6}, {0,0} }
	};

	int percDim2D = test2D.percolationDimension( nodes2D, edges2D );

	std::cout << "\nPercolation dimension for 2D example prior to adding 'bridging' links: " << percDim2D << " (should be 0)." << std::endl;
	assert( percDim2D==0 );

	// Add bridging link in one dimension.
	edges2D.emplace( std::make_pair( std::array<int,2>({3,6}), std::array<int,2>({1,0}) ) );	// Same as { {6,3}, {-1,0} }.

	percDim2D = test2D.percolationDimension( nodes2D, edges2D );

	std::cout << "Percolation dimension for 2D example after adding bridging link in one dimension: " << percDim2D << " (should be 1)." << std::endl;
	assert( percDim2D==1 );

	// Add final bridging link in the other dimension.
	edges2D.emplace( std::make_pair( std::array<int,2>({0,5}), std::array<int,2>({0,1}) ) );

	percDim2D = test2D.percolationDimension( nodes2D, edges2D );

	std::cout << "Percolation dimension for 2D example after adding bridging link in the other dimension: " << percDim2D << " (should be 2)." << std::endl;
	assert( percDim2D==2 );	

	// For this final calculation, also output some clustet information, and check a few things.
	int nClust = test2D.numClusters();
	std::cout << " - configuration has " << nClust << " clusters (should be 2)." << std::endl;
	assert( nClust==2 );

	for( auto cIndex=0; cIndex<nClust; cIndex++ )
	{
		int pDim = test2D.clusterPercDim(cIndex);
		std::set<int> cluster = test2D.cluster(cIndex);

		std::cout << " - cluster " << cIndex << ":\t";
		for( auto &n : cluster ) std::cout << n << " ";
		std::cout << "(perc dim=" << pDim << ")" << std::endl;

		assert( (cluster.size()==2&&pDim==0) || (cluster.size()==7&&pDim==2) );		// Two clusters in some order, one dim 0 and one (at this point) dim 2.
	}
	
	//
	// Another 2D example taken stright from the beadSpring code, where it initially failed.
	//
	std::cout << "\nTesting 2D network taken from off-lattice simulation." << std::endl;

	std::set<int> nodesBS{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
	std::map< std::array<int,2>, std::array<int,2> > edgesBS {
		{ { 0, 7}, { 0, 0} },
		{ { 0, 9}, { 0, 0} },
		{ { 0,12}, { 0, 0} },
		{ { 0,16}, { 0, 0} },
		{ { 1, 5}, { 0, 0} },
		{ { 1, 6}, { 0, 0} },
		{ { 2,11}, { 0, 1} },
		{ { 2,14}, { 0, 0} },
		{ { 2,15}, { 0, 1} },
		{ { 3, 8}, { 0, 0} },
		{ { 6, 7}, { 1, 0} },
		{ { 6,13}, { 0, 0} },
		{ { 7,12}, { 0, 0} },
		{ { 7,13}, {-1, 0} },
		{ { 9,10}, {-1, 0} },
		{ { 9,12}, { 0, 0} },
		{ { 9,18}, { 0, 0} },
		{ {10,12}, { 1, 0} },
		{ {10,16}, { 1, 0} },
		{ {10,17}, { 1, 0} },
		{ {11,14}, { 0,-1} },
		{ {16,18}, { 0, 0} },
		{ {16,19}, { 0, 0} }
	};

	// Get the percolation dimension.
	int pDim = test2D.percolationDimension( nodesBS, edgesBS );
	std::cout << "Percolation dimension = " << pDim <<"; expected 0 - " << ( pDim==0 ? "PASS" : "FAIL" ) << std::endl;
	assert( pDim==0 );
	
	// Display all clusters.
	nClust = test2D.numClusters();
	std::cout << " - configuration has " << nClust << " clusters (should be 4)." << std::endl;
	assert( nClust==4 );

	for( auto cIndex=0u; cIndex<nClust; cIndex++ )
	{
		int pDim = test2D.clusterPercDim(cIndex);
		std::set<int> cluster = test2D.cluster(cIndex);

		std::cout << " - cluster " << cIndex << ":\t";
		for( auto &n : cluster ) std::cout << n << " ";
		std::cout << "(perc dim=" << pDim << ")" << std::endl;
	}
}

