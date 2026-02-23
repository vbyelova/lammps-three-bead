/**
 * @file percPerBC.h
 * @author David Head
 * @brief Templated class that implements the percolation identification algorithm of Livraghi et al. (2021).
 */

#ifndef _UD_PERC_PER_BC_H
#define _UD_PERC_PER_BC_H


//
// Includes.
//
#include <iostream>
#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <cassert>
#include <queue>


/**
 * @brief Implememnts the percolation identification algorithm of Livraghi et ai. (2021).
 *
 * Implements the percolation identification algorithm of Livraghi et al., J. Chem. Theory Comput. 17(10) (2021).
 * This includes implementations of the "dim()" and "is_independent()" functions described in the article, renamed
 * here to "spanningDimension()" and "isIndependent()".
 *
 * Templated to the spatial / embedding / maximum dimension, i.e. the dimension of the space in which the objects exist.
 *
 * A Python prototype was developed prior to this C++ version.
 */
template< int dim >
class percPerBC
{

private:

	/// Alias for vector of set-valued vectors.
	using intVecSet = std::multiset< std::array<int,dim> >;

	/// Performs tests on a set of integer-valued vectors for the spanning dimension method.
	void __testSpanningDimension( intVecSet vectors, int answer, bool verbose=false ) const;

	/// Performs tests on a set of integer-valued vectors for the independence method.
	void __testIsIndependent( bool answer, intVecSet vecSet, std::array<int,dim> trailvec, bool verbose=false ) const;

	/// Removes zero rows from the matrix; changes the new number of rows in-place.
	void __removeZeroRows( std::vector<int> &matrix, int nCols, int &nRows ) const;

protected:

	/// Set of clusters and their spanning dimension. Constructed when the algorithm is run.
	std::map< std::set<int>, int > _clusters;

	/**
	 * @brief Returns the dimension spanned by the given set of integer-valued vectors.
	 * 
	 * Uses a modified version of Gaussian elimination to determine the dimensionality of the space spanned
	 * by the given set of integer-valued vectors. This function was called "dim()" and referred to as the
	 * algebraic dimension in the article of Livraghi *et al.*
	 *
	 * @param vecSet `std::multiset<>` of `std::array<int,dim>`'s with `dim` the spatial dimension.
	 * @return Spanning dimension as an integer.
	 */
	int spanningDimension( intVecSet vecSet ) const;

	/**
	 * @brief Returns `true` if the vector is independent of the given set, else `false`.
	 *
	 * Uses the spanningDimension() method to determine if the solitary vector is independent of the
	 * set of vectors also provided. Essentially, this just checks if the spanning dimension increases
	 * when the vector is included in the set.
	 *
	 * @param vecSet `std::multiset<>` of `std::array<int,dim>`'s with `dim` the spatial dimension.
	 * @param trialVec `std::array<int,dim>` of vector to be tested for independence.
	 * @return `true` if not linearly dependent, else `false`.
	 */
	 bool isIndependent( intVecSet vecSet, std::array<int,dim> trialVec ) const;

public:

	// Constructors / destructor.
	percPerBC() {};
	virtual ~percPerBC() {};


	//
	// The percolation identification algorithm.
	//

	/**
	 * @brief Returns the percolation dimension for the given connectivity graph.
	 * 
	 * Calculates the number of dimensions over which the given nodes and connectivity (edges)
	 * percolates, using the algorithm of Livraghi *et al.* (2021). The vertices are assumed to
	 * be integers, and the connectivity is specified by a `std::map` in which the keys are directed
	 * edges between connected nodes as a `std::array<int,2>`, and the value is the periodic shift,
	 * *i.e.* a `std::array<int,dim>` with the shifts in units of system size for that connection.
	 *
	 * After returning, the clusters of connected nodes can also be accessed *via* getter methods.
	 * 
	 * Note that these periodic shifts are called 'labels' in the Livraghi *et al.* paper.
	 *
	 * @param nodes Nodes (vertices) as a `std::set<int>`.
	 * @param edges `std::map` with node pairs as `std::array<int,2>` keys, and periodic shifts as `std::array<int,dim>` values.
	 * @return Integer in the range 0 to `dim` inclusive.
	 */
	int percolationDimension( std::set<int> vertices, std::map< std::array<int,2>, std::array<int,dim> > edges );

	/// Number of clusters after the algorithm has completed.
	int numClusters() const noexcept { return _clusters.size(); }

	/// Set of nodes in a given cluster. Will check index.
	std::set<int> cluster( int index ) const;

	/// Percolation dimension of the given cluster; same index as cluster(), and is also checked.
	int clusterPercDim( int index ) const;

	//
	// Debug methods.
	//

	/**
	 * @brief Tests the 'spanning dimension' function.
	 * 
	 * Runs some tests on the protected method spanningDimension() to check it works, at least for low dimensions.
	 * Will give an assertion error if any test fails. Will also provide verbose outut if the flag is selected.
	 * 
	 * @param verbose Option for verbose output; defaults to false.
	 */
	void debugSpanningDimension( bool verbose=false );
 

	/**
	 * @brief Tests the 'is independent' function.
	 * 
	 * Runs some tests on the protected method isIndependent() to check it works, at least for low dimensions.
	 * Will give an assertion error if any test fails. Will also provide verbose outut if the flag is selected.
	 * 
	 * @param verbose Option for verbose output; defaults to false.
	 */
	void debugIsIndependent( bool verbose=false );
 
};


#pragma mark -
#pragma mark Spanning dimension
#pragma mark -

// Returns the dimension spanned by the given set of integer-valued vectors.
template< int dim >
int percPerBC<dim>::spanningDimension( intVecSet vecSet ) const
{
	// Only retain non-zero vectors. Work on a copy of the passed container.
	intVecSet vSet;
	for( auto &v : vecSet )
		if( std::any_of( v.begin(), v.end(), [](int x) { return x; } ) )
			vSet.insert( v );

	// If no non-zero vectors were passed, return 0.
	if( vSet.empty() ) return 0;

	// Convert to an nxm matrix with n columns (vectors) and m rows (initially m=dim),
	// stored in flattened 1D array and indexed as matrix[row*n+col].
	int nCols = vSet.size(), nRows=dim;
	std::vector<int> matrix( nRows*nCols );

	auto col=0u;
	for( auto &v : vSet )
	{
		for( auto row=0; row<nRows; row++ ) matrix[row*nCols+col] = v[row];
		col++;
	}		

	// Local lambda function to read the flattened matrix by row and column. Cannot be used to write.
	auto mat = [&matrix,nCols]( int row, int col ) { return matrix[row*nCols+col]; };

	// // For debugging, a local lambda function to print the current matrix to stdout.
	// // Uncomment and call when required as e.g. 'printMatrix(nRows)'.
	// auto printMatrix = [&matrix,nCols]( int numRows ) {
	// 	for( auto row=0u; row<numRows; row++ ) {
	// 		for( auto col=0u; col<nCols; col++ ) std::cout << matrix[row*nCols+col] << " ";
	// 		std::cout << std::endl;
	// 	}
	// };

	// Remove any rows consisting entirely of zero values.
	__removeZeroRows( matrix, nCols, nRows );

	// Loop over all diagonal elements by column first, to eliminate values below this diagonal.
	for( auto col=0; col<std::min(nRows,nCols); col++ )
	{
		// If diagonal element is zero, need some form of pivoting.
		while( !mat(col,col) )
		{
			// Try to find the column on the same row with the first non-zero element.
			auto nZeroCol = col;
			while( !mat(col,nZeroCol) )
			{
				nZeroCol++;
				if( nZeroCol==nCols ) break;
			}

			if( nZeroCol<nCols )
			{
				// If exists, swap columns so the diagonal is now non-zero.
				for( auto row=0; row<nRows; row++ )
					std::swap( matrix[row*nCols+col], matrix[row*nCols+nZeroCol] );		// No std::swap for the lambda function.
			}
			else
			{
				break;			// Row of zeros; break loop and remove below.
			}
		}

		// Remove any new zero rows and reduce nRows accordingly; check loop counter as nRows may decrease.
		__removeZeroRows( matrix, nCols, nRows );
		if( col>=nRows ) break;	

		// Now know beta is not equal to zero. Perform the elimination on all rows below the diagonal.
		for( auto row=col+1; row<nRows; row++ )
		{
			if( mat(row,col) )  						// No need to eliminate target element if already zero.
			{
				auto
					alpha = mat(row,col), 				// Non-zero.
					beta  = mat(col,col);  				// Must use current value of the diagonal.

				// Multiply row 'col' by alpha.
				for( auto i=0; i<nCols; i++ ) matrix[col*nCols+i] *= alpha;

				// Multiply row 'row' by beta, and subtract off row 'col'.
				for( auto i=0; i<nCols; i++ ) matrix[row*nCols+i] = beta * mat(row,i) - mat(col,i);
			}
		}
	}
	
	// The spanning dimension is the number of rows that do not consist entirely of zeros.
	// Thought (10/10/25): Why not use __removeZeroRows() and return the number of rows?
	int spanDim = 0;
	for( auto row=0; row<nRows; row++ )
	{
		bool allZero = true;
		for( auto col=0; col<nCols; col++ )
			if( mat(row,col) )
			{
				allZero = false;
				break;
			}
		
		if( !allZero ) spanDim++;
	}

	return spanDim;
}

// Remove all all-zero rows (actually shuffle rows in-place and reduce the nRows accordingly).
// This method - repeating from the start after each zero row is dealth with - is perhaps over-safe and could be optimised.
template< int dim >
void percPerBC<dim>::__removeZeroRows( std::vector<int> &matrix, int nCols, int &nRows ) const
{
	// Local lambda function to read the flattened matrix by row and column. Cannot be used to write.
	auto mat = [&matrix,nCols]( int row, int col ) { return matrix[row*nCols+col]; };

	bool noZeroRows;
	do
	{
		noZeroRows = true;
		for( auto row=0; row<nRows; row++ )
		{
			// Is this a row of all zeros?
			bool zeroRow = true;
			for( auto col=0; col<nCols; col++ )
				if( mat(row,col) )
				{
					zeroRow = false;
					break;
				}
			
			// If so, shuffle the rows and reduce nRows by one.
			if( zeroRow )
			{
				for( auto r=row; r<nRows-1; r++ )
					for( auto c=0; c<nCols; c++ )
						matrix[r*nCols+c] = mat(r+1,c);		// No need to copy zero row to row nRows-1 as it will be ignored anyway once nRows is decremented.

				// Decrement nRows plus sanity check.
				if( --nRows==0 )
					throw std::runtime_error( "percPerBC::spanningDimension(): Somehow had all zero rows after removing all zero vectors from the start!" );

				// Skip the remainin rows and re-check from the beginning.
				noZeroRows = false;
				break;
			}
		}
	} while( !noZeroRows );
}

// Tests one spanning dimension calculation, with an option for verbose output to std::cout.
template< int dim >
void percPerBC<dim>::__testSpanningDimension( intVecSet vectors, int answer, bool verbose ) const
{
	// If verbose output required, output a single line with all of the input vectors.
	if( verbose )
	{
		std::cout << "Testing set";

		for( auto &vec : vectors )
		{
			// Each vector delimited by commas, surrounded by brakets.
			for( auto it=vec.begin(); it!=vec.end(); it++ )
				std::cout << (it==vec.begin() ? " (" : "," ) << *it;

			std::cout << ")";
		}
	}

	// Get the result form spanningDimension().
	int d = spanningDimension( vectors );

	// Verbose output about the result.
	if( verbose )
		std::cout << " - got result " << d << ", expected " << answer << " - " << (d==answer?"PASS":"FAIL") << "." << std::endl;
	
	// Assert the required result.
	assert( d==answer );
}


#pragma mark -
#pragma mark Independence
#pragma mark -

// Checks if the trial vector is independent of the given set of vectors, all assumed to be integer valued.
template< int dim >
bool percPerBC<dim>::isIndependent( intVecSet vecSet, std::array<int,dim> trialVec ) const
{
	// Create a new set with the trial vector added.
	intVecSet extendedSet(vecSet);
	extendedSet.insert( trialVec );

	// Independent if the spanning dimension with the vector is different (one more than) than without.
	return( spanningDimension(vecSet) != spanningDimension(extendedSet) );
}

// Tests one independence calculation, with an option for verbose output to std::cout.
template< int dim >
void percPerBC<dim>::__testIsIndependent( bool answer, intVecSet vecSet, std::array<int,dim> trialVec, bool verbose ) const
{
	// If verbose output required, output a single line with all of the input vectors and the trual vector.
	if( verbose )
	{
		std::cout << "Testing independence of";
		for( auto it=trialVec.begin(); it!=trialVec.end(); it++ )
			std::cout << (it==trialVec.begin() ? " (" : "," ) << *it;

		std::cout << ") to the set";
		for( auto &vec : vecSet )
		{
			// Each vector delimited by commas, surrounded by brakets.
			for( auto it=vec.begin(); it!=vec.end(); it++ )
				std::cout << (it==vec.begin() ? " (" : "," ) << *it;

			std::cout << ")";
		}
	}

	// Get the result form isIndependent().
	bool result = isIndependent( vecSet, trialVec );

	// Verbose output about the result.
	if( verbose )
		std::cout << std::boolalpha << " - got result '" << result << "', expected '" << answer << "' - " << (result==answer?"PASS":"FAIL") << "." << std::endl;
	
	// Assert the required result.
	assert( result==answer );
}


#pragma mark -
#pragma mark Percolation dimension
#pragma mark -

// The main algorithm for determing in how many dimensions a set of nodes percolates, if any.
template< int dim >
int percPerBC<dim>::percolationDimension( std::set<int> nodes, std::map< std::array<int,2>, std::array<int,dim> > edges )
{
	//
	// Sanity check. Check all nodes referenced in the edges are in the set of nodes.
	//
	for( auto &edge : edges )
		if(
			std::find(nodes.begin(),nodes.end(),edge.first[0]) == nodes.end()
			||
			std::find(nodes.begin(),nodes.end(),edge.first[1]) == nodes.end()
		)
			throw std::runtime_error( "FATAL in percPerBC::percolationDimension(): At least one edge had a node that is not in the node set." );

	// A dim-dimensional array of zeros.
	std::array<int,dim> zeroArray;
	for( auto k=0u; k<dim; k++ ) zeroArray[k] = 0;

	//
	// Percolation detection algorithm.
	//

	// Initialise the containers.
	std::set<int> visited;								// Node indices that have been visited.
	std::map< int, std::array<int,dim> > distance;		// 'Distance' per node when added to 'visited'.
	_clusters.clear();									// Clusters and their spanning dimension; persistent.

	// Main loop.
	for( auto &n : nodes )
	{
		// If already visited this node at some point during the search, skip.
		if( std::find(visited.begin(),visited.end(),n) != visited.end() ) continue;

		// Start breadth-first search from this node. Same notation ("B","C","Q") as the SI of the Livraghi et al. article.
		std::set<int> C;
		std::multiset< std::array<int,dim> > B;

		// Initialise the queue with this one node.
		std::queue< std::pair< int, std::array<int,dim> > > Q;
		Q.push( {n,zeroArray} );

		// Breadth-first search expanding outwards from this node.
		while( !Q.empty() )
		{
			// Get next item from the queue.
			auto n1   = Q.front().first;		// 'front' = next in the queue.
			auto dist = Q.front().second;		// Note 'front()' does not remove the item, just accessess it.
			Q.pop();							// Remove the item now; pop() does not return anything.

			// Is this a copy of an already-visited node? 
			if( std::find(visited.begin(),visited.end(),n1) != visited.end() )
			{
				// Calculate 'distance' from the original copy.
				std::array<int,dim> newDist;
				for( auto k=0u; k<dim; k++ ) newDist[k] = distance[n1][k] - dist[k];

				// If independent to the existing distance vectors, add to the set.
				if( isIndependent(B,newDist) ) B.insert( newDist );
			}
			else
			{
				// First time we have encountered this node. Update containers and add connected nodes to queue.
				visited.insert( n1 );
				distance[n1] = dist;
				C.insert( n1 );

				// Add all connected nodes to the queue, i.e. breadth-first search. Note edge.first is the directed
				// edge as a 2-vector, and edge.second is the periodic shift for this edge in the given direction.
				// Check edge in both directions, but flip the periodic shift if traversing in reverse order.
				// (not entirely clear from the paper, but found at least one test case that failes without this).
				for( auto &edge : edges )
				{
					std::array<int,dim> newDist(dist);		// () rather than {} for avoid errors with g++.

					if( edge.first[0]==n1 )
					{
						for( auto k=0u; k<dim; k++ ) newDist[k] += edge.second[k];
						Q.push( {edge.first[1],newDist} );
					}

					if( edge.first[1]==n1 )
					{
						for( auto k=0u; k<dim; k++ ) newDist[k] -= edge.second[k];			// Note sign flip for the periodic shift.
						Q.push( {edge.first[0],newDist} );
					}
				}
			}
		}

		// Have now searched through all vertices starting from n. Add the set to the list of components, with the
		// spanning dimension of all of the distance vectors.
		_clusters[C] = spanningDimension(B);
	}

	// Now have all components, each with a spanning dimension. The percolation dimension is the largest of all of them.
	int percDim = 0;
	for( auto iter=_clusters.begin(); iter!=_clusters.end(); iter++ ) percDim = std::max( percDim, iter->second );

	return percDim;
}

#pragma mark -
#pragma mark Getters
#pragma mark -

// Set of nodes in a cluster, with a range check of the cluster index.
template< int dim >
std::set<int> percPerBC<dim>::cluster( int index ) const
{
	// Range check for the index.
	if( index<0 || index>=_clusters.size() ) throw std::out_of_range( "bad cluster index in percPerBC<>::cluster()" );

	// Iterate through clusters and count as you go.
	int i=0;
	auto iter = _clusters.begin();
	while( i<index ) { iter++; i++; }

	return iter->first;
}

// Percolation dimension of the given cluster, with a range check.
template< int dim >
int percPerBC<dim>::clusterPercDim( int index ) const
{
	// Range check for the index.
	if( index<0 || index>=_clusters.size() ) throw std::out_of_range( "bad cluster index in percPerBC<>::cluster()" );

	// Iterate through clusters and count as you go.
	int i=0;
	auto iter = _clusters.begin();
	while( i<index ) { iter++; i++; }

	return iter->second;}

#endif
