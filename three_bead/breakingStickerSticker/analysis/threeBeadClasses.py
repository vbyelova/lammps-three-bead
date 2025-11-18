from collections import defaultdict

class Particle():
    """ A  particle with id, type, molecule id and number of intramolecular
        sticky bonds."""
    def __init__(self):
        self.properties = []

class Molecule():
    """A molecule with an id, list of separations at each timestep and a 
    """
    def __init__(self):
        self.molNum = 0
        self.separation = []
        self.bondedMol = defaultdict(list)
