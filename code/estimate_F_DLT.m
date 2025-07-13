function [Fn, minSingVal, normMv] = estimate_F_DLT(x1s, x2s)
    % compute sizes
    [m, n] = size(x1s);
    assert(m == 3)
    assert(all(size(x2s) == [3 n]))

    % % initialise M
    % M = zeros(n, 9);

    % % store columnwise Kronecker product (outer product) in rows of M
    % for i = 1:n
    %     M(i, :) = kron(x2s(:, i), x1s(:, i))';
    % end
    matrices = arrayfun(@(i) kron(x2s(:, i), x1s(:, i))', ...
                        1:n, 'UniformOutput', false);
    M = cat(1, matrices{:});

    % compute SVD
    [~, S, V] = svd(M);

    % extract fundamental matrix from solution
    v = V(:, end);
    Fn = reshape(v, [3, 3])';

    % compute smallest singular value and ||Mv||
    singVals = diag(S);
    minSingVal = singVals(end);
    normMv = norm(M * v);

    if n >= 9
        % verify min_i sigma_i = ||Mv|| if n >= 9
        assert(abs(minSingVal - normMv) < 1e-10);
    else
        % verify ||Mv|| = 0 if n < 9
        assert(abs(normMv) < 1e-10);
    end
end