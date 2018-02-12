# Connected-Components-using-MPI

## Objective

Your task in this assignment is to implement efficient C++/MPI function [connected_components](https://en.wikipedia.org/wiki/Connected_component_(graph_theory)). We make several assumptions:

Undirected graph on which we are operating is too large to be represented in the memory of a single compute node.

We have p = q * q ranks available.

The graph has n nodes, and we have that q divides n.

The graph is represented by the adjacency matrix A, in which 1 indicates edge and 0 means no edge.

Adjacency matrix A is 2D-decomposed using q by q row-wise grid of ranks.

## Description

Arguments of the connected_components are as follows:

**A** : adjacency matrix, row-wise block of size n/q by n/q.

**n** : total number of nodes in the graph.

**q** : dimension of the rank grid (p = q * q).

**out** : path to the output file where the assignment of nodes to the connected components should be stored (see below).

**comm** : communicator with p = q * q ranks to work with.

Function must return the total number of connected components found.

