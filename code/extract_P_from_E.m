function P2s = extract_P_from_E(E)
    % ensure det(UV') > 0
    [U, ~, V] = svd(E);
    if det(U * V') < 0
        V = -V;
    end

    % compute possible P2s
    W = [0 -1 0;
         1  0 0;
         0  0 1];
    u3 = U(:, 3);
    P2s = {[U * W * V' u3];
           [U * W * V' -u3];
           [U * W' * V' u3];
           [U * W' * V' -u3];};
end