# Authors: Elena Drake-Knight, Esther Fatehi, Jose Murrieta
# Date: 3/10/2018
# Description: This program reads a file specified in the command line argument.
#   The file contains an instance of the travelling salesman problem. Each line
#   defines a city and each line has three numbers sperated by white space. The
#   first number is the city identifier, the second is the x-coordinate, and the
#   third is the y-coordinate. This program then uses Christofides algorithm to
#   solve the tsp and outputs the solution to a file with '.tour' appended to the
#   input file name.

from sys import argv
from math import sqrt

class City(object):
    def __init__(self, c_id=None, x=None, y=None):
        self.id = c_id
        self.x = x
        self.y = y

def distance(city_a, city_b):
    return round(sqrt((city_a.x - city_b.x)**2 + (city_a.y - city_b.y)**2))

def mst_prim(adj_matrix):
    Q = list(adj_matrix.keys())
    keys = {c_id: float('inf') for c_id in Q}
    keys[Q[0]] = 0
    parents = {c_id: None for c_id in Q}
    T = {}
    while(keys):
        u = min(keys, key = keys.get)
        T[u] = {}
        if parents[u]:
            T[u][parents[u]] = T[parents[u]][u] = keys[u]
        del keys[u]
        for city, dist in adj_matrix[u].items():
            if((city in keys) and adj_matrix[u][city] < keys[city]):
                parents[city] = u
                keys[city] = dist
    return T

def getOddDegVerts(adj_matrix):
    O = set()
    for city in adj_matrix.keys():
        edges = len(adj_matrix[city])
        if edges % 2 != 0:
            O.add(city)
    return O

def minWeightMatching(G, O):
    M = set()
    while O:
        u = O.pop()
        length = float("inf")
        closest = 0
        for v in O:
            if u != v and G[u][v] < length:
                length = G[u][v]
                closest = v
        M.add((u, closest, length))
        O.remove(closest)
    return M

def combineGraphs(T, M):
    H = T
    while M:
        u, v, dist = M.pop()
        H[u][v] = H[v][u] = dist
    return H


# open input file:
fileName = argv[1]
inputFile = open(fileName, 'r')

# create a complete graph G (adjacency matrix) from the file
G = {}
cities = [] 
nextLine = inputFile.readline()
while (nextLine):
    c_id, c_x, c_y = nextLine.split()
    newCity = City(c_id, int(c_x), int(c_y))
    cities.append(newCity)
    G[newCity.id] = {}

    # add distances to new city to adjacency matrix
    for city in cities:
        G[city.id][newCity.id] = G[newCity.id][city.id] = distance(city, newCity)

    nextLine = inputFile.readline()

# close input file
inputFile.close()

# get minimum spanning tree T of G
T = mst_prim(G)

# Let O be the set of vertices with odd degree in T
O = getOddDegVerts(T)

# Find a minimum-weight perfect matching M in the induced subgraph given by the vertices from O
M = minWeightMatching(G, O)

# Combine the edges of M and T to form a connected multigraph H in which each vertex has even degree
H = combineGraphs(T, M)

# Form an Eulerian circuit in H
# TODO

# Make the circuit found in the previous step into a Hamiltonian circuit by skipping repeated vertices
# TODO

# print results in output file
#outputFile = open(fileName + '.tour', 'w')
# TODO
# close output file
#outputFile.close()
