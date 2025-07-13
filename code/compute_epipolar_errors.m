function [l, distances] = compute_epipolar_errors(F, x1s, x2s)
    % compute epipolar lines
    l = F * x1s;
    l = l ./ vecnorm(l(1:2,:)); % use implicit expansion

    % compute distances
    distances = abs(sum(l .* x2s));
end