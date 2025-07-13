function [E, S] = enforce_essential(E_approx)
    % set singular values to 0, 1, 1
    [U, S, V] = svd(E_approx);
    E = U * diag([1 1 0]) * V';
end