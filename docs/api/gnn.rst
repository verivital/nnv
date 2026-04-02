GNN Class
=========

.. rst-class:: lead

   Graph neural network class supporting GCN and GINE architectures with
   node and edge feature uncertainty via GraphStar sets.

----

Constructors
------------

.. code-block:: matlab

   gnn = GNN(layers)                                      % Layers only
   gnn = GNN(layers, A_norm)                              % GCN mode
   gnn = GNN(layers, A_norm, adj_list, E)                 % GINE mode
   gnn = GNN(layers, A_norm, adj_list, E, edge_weights)   % GINE with edge weights

- ``layers`` -- cell array of GCNLayer, GINELayer, ReluLayer, etc.
- ``A_norm`` -- normalized adjacency matrix (N x N)
- ``adj_list`` -- edge list [src, dst] (M x 2)
- ``E`` -- edge feature matrix (M x d_e) or EdgeGraphStar
- ``edge_weights`` -- optional edge weight vector

Properties
----------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Property
     - Type
     - Description
   * - ``Layers``
     - cell array
     - GCN/GINE and activation layers
   * - ``A_norm``
     - matrix
     - Normalized adjacency matrix
   * - ``adj_list``
     - matrix
     - Edge list (M x 2)
   * - ``E``
     - matrix
     - Edge feature matrix
   * - ``edge_weights``
     - vector
     - Optional edge weights

Methods
-------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``Y = gnn.evaluate(X)``
     - Evaluate GNN on node feature matrix X (N x F_in). Returns output node features.
   * - ``outputSet = gnn.reach(inputSet, reachOptions)``
     - Compute reachable set. Input/output: GraphStar.
   * - ``gnn.setGraph(A_norm, adj_list, E)``
     - Update graph structure (allows weight reuse across different graphs).

Example
-------

.. code-block:: matlab

   layers = {GCNLayer(W1, b1), ReluLayer(), GCNLayer(W2, b2)};
   gnn = GNN(layers, A_norm);

   Y = gnn.evaluate(X_test);

   GS = GraphStar(X, -eps_matrix, eps_matrix);
   reachOptions.reachMethod = 'approx-star';
   output_sets = gnn.reach(GS, reachOptions);
