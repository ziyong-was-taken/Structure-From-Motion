% use threshold as name instead of eps (eps predefined in MATLAB)
function [bestE, inlierIndices] = estimate_E_robust(x1,x2,threshold)
    bestE = [];
    largestMask = [];
    inlierFraction = 0;
    i = 0;

    % run RANSAC until 
    % - inlier fraction at least 0.2
    % - i at least 5000
    % - early exit if inlier fraction at least 0.8
    while (i < 5000 || inlierFraction < 0.2) && inlierFraction < 0.8
        % randomly sample 8 correspondences
        rand8 = randperm(size(x1, 2), 8);
        
        % compute candidate essential matrix
        [E, ~, ~] = estimate_F_DLT(x1(:, rand8), x2(:, rand8));
        [E, ~] = enforce_essential(E);
        
        % compute errors
        [~, d1] = compute_epipolar_errors(E, x1, x2);
        [~, d2] = compute_epipolar_errors(E', x2, x1);
        
        % compute consensus set mask
        inlierMask = (d1.^2 + d2.^2) / 2 < threshold^2;

        % update best model and consensus set
        if sum(inlierMask) > sum(largestMask)
            inlierFraction = mean(inlierMask);
            largestMask = inlierMask;
            bestE = E;
        end

        % increment counter
        i = i + 1;
    end
    inlierIndices = find(largestMask);
end