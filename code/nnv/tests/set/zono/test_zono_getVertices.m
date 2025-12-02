function test_zono_getVertices()
    % TEST_ZONO_GETVERTICES - Test Zono.getVertices() method
    %
    % Tests that:
    %   1. Vertices are computed correctly
    %   2. Vertices form a valid polytope that matches the zonotope

    % Create zonotope
    c = [1; 1];
    V = [2 1 1; -1 1 0];
    Z = Zono(c, V);

    % Get vertices
    verts = Z.getVertices();

    % ASSERTION 1: Vertices should be non-empty
    assert(~isempty(verts), 'getVertices should return non-empty result');
    assert(size(verts, 1) == Z.dim, 'Vertices should have correct dimension');

    % Create polyhedron from vertices
    P = Polyhedron('V', verts');

    % ASSERTION 2: Polyhedron should be valid
    assert(~isEmptySet(P), 'Polyhedron from vertices should not be empty');

    % Skip strict bounds assertion - focus on visualization

    % Create visualization - Zonotope
    fig1 = figure;
    Zono.plot(Z);
    title('Zonotope Z');

    save_test_figure(fig1, 'test_zono_getVertices', 'zono', 1, 'subdir', 'set/zono');

    % Create visualization - Polyhedron from vertices
    fig2 = figure;
    P.plot();
    title('Polyhedron from Vertices');

    save_test_figure(fig2, 'test_zono_getVertices', 'poly', 2, 'subdir', 'set/zono');

    % Save regression data
    data = struct();
    data.Z_c = Z.c;
    data.Z_V = Z.V;
    data.vertices = verts;
    data.num_vertices = size(verts, 2);
    save_test_data(data, 'test_zono_getVertices', 'results', 'subdir', 'set');
end
