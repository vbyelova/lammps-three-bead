
# Prototype for the percolation detection algorithm of Livragi et al., that is suitable for
# periodic boundaries in off-lattice systems. Call as stand-alone for a test.
#
# Adapted from Livraghi et al., J. Chem. Theory Comput. 17(10), 6449-6457 (2021).
#
# WARNING: Will sometimes fail as the Gaussian elimination has no pivoting.
# This was corrected in the C++ version but not here (which had already served its purpose as a prototype).
#


#
# Imports.
#
import queue
import numpy as np


#
# Main class.
#
class percPerBC:

	def __init__( self, dim ):
		"""Initialise with the spatial / embedding dimension, i.e. the maximum dimension of the cluster."""
		self._dim = dim

	def percolationDimension( self, V, E ):
		"""Performs the percolation detection algorithm as per Livraghi et al., J. Chem. Theory Comput. 17(10), 6449-6457 (2021).
		Returns the largest dimension of all components identified by the algorithm. If this equals the spatial/embedding dimension,
		then by the definition of the article, it has percolated. By returning the maximum dimension we also allow other definitions
		of percolated, in particular, percolated in only 1 dimension (for spatial dimensions >1)."""

		#
		# Sanity check.
		#

		# All elements of V single labels.
		for v in V:
			try:
				if len(v)!=1:
					raise ValueError( "At least one element of V is not size 1." )
			except TypeError:
				pass			# Probably means it was an atomic type, which is fine.

		# All elements of E are size 3.
		for e in E:
			if len(e) != 3:
				raise ValueError( "At least one element of E is not size 3." )

		# The vertex elements of E refer to elements in V.
		for e in E:
			if e[0] not in V or e[1] not in V:		# Note we know by now that each element of e is size 2.
				raise ValueError( "At least one edge refers to a vertex that is not in V." )

		# The 'label' (=normalised periodic shift) elements of E are of the same dimension as the embedding space.
		for e in E:
			if len(e[2]) != self._dim:
				raise ValueError( "At least one edge has a normalised periodic shift that is of a different dimension to the space." )

		# Could also use something like disjointSets.py to confirm all vertices belong to a single cluster.

		#
		# Initialise.
		#
		self._visited    = set()
		self._distance   = { v:np.zeros(self._dim,dtype="int") for v in V }
		self._components = []

		#
		# Main loop.
		#
		for v in V:

			# If already visited, skip.
			if v in self._visited:
				continue

			# Not yet visited. Start breadth-first search. Same notation as the SI of the paper.
			C = set()
			B = []

			Q = queue.SimpleQueue()
			Q.put( [v,np.zeros(self._dim,dtype="int")] )

			while not Q.empty():
				n, dist = Q.get()		# Also removes from the queue.

				# Is this a copy of a vertex? If so, calculate 'distance' from original copy.
				if n in self._visited:
					perShift = self._distance[n] - dist

					if percPerBC.isIndependent( B, perShift ):
						B.append( perShift )

				# First time we have encountered this vertex.
				else:
					self._visited.add( n )			# We have now visited this verted.
					self._distance[n] = dist		# It was found at dist == [0,0,...]
					C.add ( n )						# Also a member of this component.

					# Breadth-first search: Add to queue.
					for e in E:

						if e[0] == n:
							Q.put( [ e[1], dist + np.array(e[2]) ] )

						if e[1] == n:
							Q.put( [ e[0], dist - np.array(e[2]) ] )
			
			# Dealt with all vertices. Add to the list of components.
			self._components.append( [C,percPerBC.spanningDimension(B)] )

		# Now have all of the components. Return the largest dimension of all components.
		return max( [ component[1] for component in self._components ] )

	def __str__( self ):
		"""Informal representatiom of the current state of the algorithm."""

		msg = "percPerBC object for spatial dimension {}.".format(self._dim)

		msg += "\n - current visited set  : {}".format(self._visited)
		msg += "\n - current distance map : {}".format(self._distance) 
		msg += "\n - current components   : {}".format(self._components)

		return msg

	#
	# Class methods.
	#
	@classmethod
	def spanningDimension( cls, setOfVectors ):
		"""Returns the dimensions spanned by the given set of vectors, which are assumed to be 1, 2 or 3 dimensional
		and to be lists/numpy arrays of integers. The method used is essentially a simplified version of Gaussian elimination.
		There is no need to check for floating point issues as all elements are assumed to be integer."""
		
		#
		# Sanity check: All vectors must be the same size, and must be at least one.
		#

		# At least one vector passed. If not, just return zero rather than an error message.
		try:
			if len(setOfVectors)==0:
				return 0
		except Exception as err:
			raise ValueError( "Failure in percPerBC::spanningDimension : Could not determine length of passed set of vectors.\nError: {}".format(err) )

		# Set maximum/spatial/embedding dimension from the first item.
		try:
			dim = len( setOfVectors[0] )
		except Exception as err:
			raise ValueError( "Failure in percPerBC::spanningDimension : Could not determine dimension of vectors.\nError: {}".format(err) )

		# Make sure they are all the same length.
		for v in setOfVectors[1:]:
			if len(v) != dim:
				raise ValueError( "Failure in percPerBC::spanningDimension : Not all vectors have the same dimension." )

		# Only been checked for dims up to 3.
		if dim>3:
			raise NotImplementedError( "Failure in percPerBC::spanningDimension : Only dimensions up to 3 currently checked." )

		#
		# Extract the spanning dimension (called the 'algebraic dimension' in the article).
		#

		# Convert to numpy arrays after removing zero vectors. 
		sVectors = []
		for v in setOfVectors:
			npVec = np.array( v, dtype="int" )
			if np.any( npVec ):
				sVectors.append( npVec )

		# If they are all zero vectors, return a spanning dimension of 0.
		if not len(sVectors):
			return 0

		# Assemble a matrix with the vectors as columns.
		nVec = len(sVectors)
		M = np.zeros( shape=[dim,nVec], dtype="int" )
		for col in range(nVec):
			M[:,col] = sVectors[col]

		# Remove any rows consisting entirely of zero values.
		M = cls.removeZeroRows( M )

		# Loop over all diagonal elements by column first, to eliminate values below this diagonal.
		for col in range( min(M.shape[0],M.shape[1]) ):

			# If diagonal element is zero, need some form of pivoting.
			while( M[col,col]==0 ):

				# Find the column on the same row with the first non-zero element.
				nZeroCol = col
				while( M[col,nZeroCol]==0 ):
					nZeroCol += 1
					if nZeroCol == M.shape[1]:
						break

				if( nZeroCol < M.shape[1] ):

					# If exists, swap columns so the diagonal is now non-zero.
					for row in range(M.shape[0]):
						swap = M[row,col]
						M[row,col] = M[row,nZeroCol]
						M[row,nZeroCol] = swap

				else:
					break

			# Remove any new zero rows; check loop counter as matrix.shape[0] may have just decreased.
			M = cls.removeZeroRows( M )
			if col >= M.shape[0]:
				break

			# Now know 'beta' (see below) is not zero. Perform the elimination on all rows below the diagonal.
			for row in range(col+1,M.shape[0]):
				if( M[row,col] ):

					alpha = M[row,col]
					beta  = M[col,col]

					for i in range(M.shape[1]):
						M[col,i] *= alpha

					for i in range(M.shape[1]):
						M[row,i] = beta*M[row,i] - M[col,i]

		# The spanning dimension is the number of rows that do not consist entirely of zeros.
		M = cls.removeZeroRows( M )

		return M.shape[0]

	@classmethod
	def isIndependent( cls, setOfVectors, testVector ):
		"""Tests whether or not testVector is independent of the provided setOfVectors, by seeing in the spanning dimension is increased
		when testVector is included in the set."""

		# Convenient to call this with no vectors, in which case the test vector must be independent as long as it is non-zero.
		if len(setOfVectors)==0:
			return np.any( testVector )

		return cls.spanningDimension( setOfVectors ) != cls.spanningDimension( setOfVectors+[testVector] )

	@classmethod
	def removeZeroRows( cls, matrix ):
		"""Returns a new matrix that consists only of the rows of the input matrix whose rows have at least one non-zero value."""
		return np.delete( matrix, [ not matrix[row,].any() for row in range(matrix.shape[0]) ], 0 )
	
	@classmethod
	def debugSpanningDimension( cls, verbose=False ):
		"""Checks the spanningDimension() class method for a range of low-dimensional test cases."""

		# Tests and expected result.
		testsAndResults = [
			[	0,	[[0,0],[0,0],[0,0]]			],
			[	1,	[[1],[2]]					],
			[	1,	[[1,2],[2,4],[0,0]]			],
			[	2,	[[1,2],[-1,0]]				],
			[	3,	[[1,0,0],[0,1,0],[0,0,1]]	],
			[	2,	[[1,2,3],[4,5,6],[1,0,-1]]	],
			[	1,	[[1,2,0],[2,4,0],[-3,-6,0]]	]
		]

		# Iterate through all tests.
		for result, test in testsAndResults:

			if verbose:
				print( "Testing on '{}'; expect the result {} ...".format(test,result) )

			if percPerBC.spanningDimension( test ) != result:
				print( "WARNING: spanningDimension() failed for '{}'; did not return {}.".format(test,result) )
			else:
				if verbose:
					print( " - passed." )
	
	@classmethod
	def debugIsIndependent( cls, verbose=False ):
		"""Tests the isIndependent() class method for low-dimensional test cases."""

		# Tests and expected results.
		testsAndResults = [
			[ 	True, 	[[0,1],[0,2]],		[1,0]	],
			[	False,	[[0,1],[1,0]],		[1,0]	],
			[	False,	[[1]],				[2]		],
			[	True,	[[0]],				[-1]	],
			[	True,	[[1,2,3],[1,3,2]],	[1,0,0]	],
			[	False,	[[1,2,3],[1,3,2]],	[1,2,3]	]	
		]

		# Iterate through all of the provided tests.
		for result, set, vec in testsAndResults:

			if verbose:
				print( "Testing on set '{}' and new vector '{}'; expect the result '{}' ...".format(set,vec,result) )
			
			if percPerBC.isIndependent( set, vec ) != result:
				print( "WARNING: isIndependent() failed for set '{}' and vector '{}'.".format(set,vec) )
			else:
				if verbose:
					print( " - passed." )


#
# Test when called directly.
#
if __name__ == "__main__":

	#
	# Test the auxiliary class methods. Will only output messages in case of failure (unless 'verbose' flag set to 'True').
	#
	percPerBC.debugSpanningDimension( verbose=False )
	percPerBC.debugIsIndependent( verbose=False )

	#
	# 1D example. This is Fig. 4 of the paper.
	#
	V = [ 0, 1, 2, 3, 4, 5, 6 ]
	E = [
		[ 2, 5, [ 0] ],		# Format: 'from' vertex, 'to' vertex, periodic shift in units of whole system sizes (confusingly called 'label' in the paper).
		[ 5, 6, [ 0] ],
		[ 6, 4, [ 0] ],
		[ 4, 0, [ 0] ],
		[ 0, 1, [ 0] ],
		[ 0, 3, [ 0] ],
		[ 0, 2, [ 0] ],
		[ 3, 1, [ 1] ],		# Same as [1,3,[-1]].
		[ 2, 4, [-1] ]		# Same as [4,2,[ 1]].
	]

	eg_1d = percPerBC( 1 )
	print( "\n1D example; same example as in the article." )
	print( "Maximum spanning dimension of a component: {}".format(eg_1d.percolationDimension(V,E) ) )
	print( eg_1d )

	#
	# 2D Example. Network arranged as (with all links):
	#
	#    |    5              |  
	# - - - - - - - - - - - - - -
	#    |    |              |  
	#    |    0       2      |  
	#    |      \   /   \    |  
	#  3 | - 6 -  1      3 - | 6 
	#    |        |          |  
	#    |    5 - 4   7 - 8       |  
	#    |    |              | 
	# - - - - - - - - - - - - - -
	#    |    0              |  
	#
	V = [ 0, 1, 2, 3, 4, 5, 6, 7, 8 ]
	E = [
		[ 0, 1, [0,0] ],
		[ 1, 2, [0,0] ],
		[ 2, 3, [0,0] ],
		[ 1, 4, [0,0] ],
		[ 4, 5, [0,0] ],
		[ 7, 8, [0,0] ],
		[ 1, 6, [0,0] ],
		[ 3, 6, [1,0] ],	# Comment-out both of the last two lines to not percolate, even in one dimension.
		[ 0, 5, [0,1] ]		# Comment-out this line only to percolate in only one dimension.
	]

	eg_2d = percPerBC( 2 )
	print( "\n2D example; see code for connectivity." )
	print( "Maximum spanning dimension of a component: {}".format(eg_2d.percolationDimension(V,E) ) )
	print( eg_2d )
