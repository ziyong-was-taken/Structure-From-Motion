function threeDPoint = triangulate_3D_point_DLT(xs, Ps)
    % xs: cell array of 2D points in homogeneous coordinates
    % Ps: cell array of corresponding camera matrices

    % helper function to generate cross product matrix
    cross_matrix = @(x, P) cross(repmat(x, 1, size(P, 2)), P);

    % generate M by stacking cross product matrices
    matrices = cellfun(cross_matrix, xs, Ps, 'UniformOutput', false);
    M = cat(1, matrices{:});

    % compute solution using SVD
    [~, ~, V] = svd(M);
    threeDPoint = V(:, end);
end