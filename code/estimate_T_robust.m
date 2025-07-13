function bestT = estimate_T_robust(x, X, R, threshold)
    % RANSAC
    bestT = zeros(3, 1);
    largestMask = [];
    for i = 1:5000
        % randomly sample 2 correspondences
        rand2 = randperm(size(x, 2), 2);
        
        % compute candidate translation vector
        T = estimate_T_DLT(x(:, rand2), X(:, rand2), R);
        
        % compute projection using T
        P = [R T];
        xMeas = pflat(x);
        xProj = pflat(P * X);
        
        % compute consensus set mask
        inlierMask = vecnorm(xMeas - xProj) < threshold;

        % update best model and consensus set
        if sum(inlierMask) > sum(largestMask)
            largestMask = inlierMask;
            bestT = T;
        end   
    end
end