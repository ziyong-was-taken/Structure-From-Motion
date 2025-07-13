function T = estimate_T_DLT(x, X, R)
    % estimates T using DLT
    % x: 2D points in homogeneous coordinates
    % X: 3D points in homogeneous coordinates
    % R: rotation matrix
    [three,n1] = size(x);
    [four,n2] = size(X);
    assert(three == 3, '2D points must be in homogeneous coordinates.');
    assert(four == 4, '3D points must be in homogeneous coordinates.');
    assert(n1 == n2, 'Number of 3D and 2D points must match.');
    assert(all(size(R) == [3 3]), 'Rotation matrix must be 3x3.');
    assert(all(abs(R' * R - eye(3)) < 1e-10, 'all'), 'Rotation matrix must be orthogonal.');

    % helper function to generate rows of M
    function threeRows = generate_rows(x, X)
        crossx = [  0  -x(3)  x(2);
                  x(3)    0  -x(1);
                 -x(2)  x(1)    0];
        threeRows = [crossx cross(x, R*X)];
        assert(all(size(threeRows) == [3 4]));
    end

    % generate M by stacking rows
    matrices = arrayfun(@(i) generate_rows(x(:, i), X(1:3, i)), ...
                        1:n1, 'UniformOutput', false);
    M = cat(1, matrices{:});
    
    % compute solution using SVD
    [~, ~, V] = svd(M);
    sol = pflat(V(:, end));
    T = sol(1:3);

    % choose sign with most points in front of camera
    x1 = [R T] * X;
    x2 = [R -T] * X;
    if nnz(x2(3, :) > 0) > nnz(x1(3, :) > 0)
        T = -T;
    end
end