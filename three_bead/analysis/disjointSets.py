#
# Determines disjoint sets given a set of nodes and edges, using union-by-rank and path compression.
# 

# Code by Dr. David Head.

#
# Imports.
#

class disjointSets( object ):

    def __init__( self ):
        """Determines all disjoint sets provided nodes and edges. Uses union by rank and path compression."""
        self._sets = []

    def getDisjointSets( self, nodes, edges, noSingletons=False ):
        """Determines and returns (as a list of sets in no particular order) all disjoint sets from the given
        nodes and edges. If "noSingletons' is True, will not include sets of size 1; defaults to False."""

        # Initialise sets as dictionary of dictionaries of isolated nodes.
        self._sets = { n : {"parent":n,"rank":0} for n in nodes }

        # Loop through all edges.
        for edge in edges:
 
            # Get roots for both clusters containing the nodes at both ends of this edge.
            root1 = self._findRoot( edge[0] )
            root2 = self._findRoot( edge[1] )

            # If already in the same cluster, nothing to do for this edge.
            if root1==root2: continue

            # If on different clusters and have different rank, add (by re-assigning a parent) the one with the
            # lower rank to the one with the higher rank.
            if self._sets[root1]["rank"] > self._sets[root2]["rank"]:
                self._sets[root2]["parent"] = self._sets[root1]["parent"]
                continue

            if self._sets[root2]["rank"] > self._sets[root1]["rank"]:
                self._sets[root1]["parent"] = self._sets[root2]["parent"]
                continue

            # If still here, both have the same rank. Add one to the other and increase the rank.
            self._sets[root1]["parent"]  = root2
            self._sets[root2]["rank"  ] += 1

        # Convert to a dictionary of nodes keyed by the parent.
        completeSets = {}
        for node in self._sets.keys():
            parent = self._findRoot( self._sets[node]["parent"] )

            if parent not in completeSets.keys():
                completeSets[parent] = { node }
            else:
                completeSets[parent].add(node)

        # Return with the (list of) sets, skipping singletons if requested.
        return [ completeSets[n] for n in completeSets.keys() if ( not noSingletons or len(completeSets[n])>1 ) ]

    def _findRoot( self, n ):
        """Finds the root for the given node, applying path compression as it does so."""
        if self._sets[n]["parent"] != n:
            self._sets[n]["parent"] = self._findRoot( self._sets[n]["parent"] )
        
        return self._sets[n]["parent"]

    #
    # Class methods.
    #
    @classmethod
    def check( cls, sets, nodes, edges ):
        """Checks the sets for basic consistency. Uses naive nested loops so will be slow for large problems."""

        # Check all nodes appear exactly once.
        for node in nodes:
            count = 0
            for set in sets:
                if node in set: count += 1
            if count != 1:
                    return false

        # Check nodes at the ends of all edges appear in the same set.
        for edge in edges:
            for set in sets:
                if (edge[0] in set and edge[1] not in set) or (edge[1] in set and edge[0] not in set):
                    return False

        return True

#
# If called from command line: Test.
#
if __name__ == "__main__":
    
    #
    # Example using letters.
    #
    nodes = [ "a", "b", "c", "d", "e", "f" ]
    edges = [ ["a","b"], ["c","d"], ["a","f"] ]

    print( "Letter example:" )
    for i, set in enumerate( disjointSets().getDisjointSets(nodes,edges) ):
        print( "Set {0}: {1}".format(i,set) )
    
    #
    # Random integers.
    #
    import random
    Nnodes, Nedges = 30, 20

    nodes = range(Nnodes)
    edges = [ [random.randrange(Nnodes),random.randrange(Nnodes)] for e in range(Nedges)]

    # Get the disjoint sets.
    print( "\nRandom integer example:" )
    print( "Nodes: {}".format(nodes) )
    print( "Edges: {}".format(edges) )
    for i, set in enumerate( disjointSets().getDisjointSets(nodes,edges) ):
        print( "Set {0}: {1}".format(i,set) )

    #
    # Check.
    #
    if disjointSets.check( disjointSets().getDisjointSets(nodes,edges), nodes, edges ):
        print( "Last example passed test." )
    else:
        print( "Last example FAILED test." )
