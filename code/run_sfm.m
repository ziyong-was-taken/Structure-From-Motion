function run_sfm(dataset)
    % setup
    rng(42);
    vl_setup;
    [K, imgNames, initPair, inlierThreshold] = get_dataset_info(dataset);
    savefile = ['data/' num2str(dataset) '/precomputed.mat'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 0: Precompute keypoints and matches %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Computing SIFT keypoints and matches...");

    % preload all images
    images = cellfun(@imread, imgNames, 'UniformOutput', false);

    % precompute SIFT features and descriptors for each image
    [features, descriptors] = cellfun(@(img) vl_sift(single(rgb2gray(img))), ...
                                      images, 'UniformOutput', false);
    keypoints = cellfun(@(f) f(1:2, :), features, 'UniformOutput', false);
    assert(size(keypoints, 2) == size(descriptors, 2));

    % save(savefile, "keypoints", "descriptors");

    % % load features and descriptors
    % precomputed = load(savefile);
    % keypoints = precomputed.keypoints;
    % descriptors = precomputed.descriptors;

    % store correspondences (normalised matches)
    % i:   images i and i+1
    % end: images i1 and i2
    correspondences = cell(1, length(images));
    for i = 1:length(images)
        curr = i;
        next = i+1;
        if i == length(images)
            curr = initPair(1);
            next = initPair(2);
        end
        
        % extract matching points
        matches = vl_ubcmatch(descriptors{curr}, descriptors{next});
        x1 = keypoints{curr}(:, matches(1, :));
        x2 = keypoints{next}(:, matches(2, :));

        % store matches between initial pair
        if i == length(images)
            matchesInitPair = matches;
        end

        % normalise image points
        x1n = K \ [x1; ones(1, size(x1, 2))];
        x2n = K \ [x2; ones(1, size(x2, 2))];

        % store correspondences
        correspondences{i} = cat(3, x1n, x2n);
    end

    % save(savefile, "correspondences", "matchesInitPair", "-append");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Steps 1-2: Calculate relative orientations & Upgrade rotations %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Calculating relative orientations and upgrading rotations...");

    % % load precomputed data
    % precomputed = load(savefile);
    % correspondences = precomputed.correspondences;

    % helper function to triangulate all points
    function X = triangulate(P1, P2, x1, x2)
        X = arrayfun(@(i) triangulate_3D_point_DLT({x1(:, i), x2(:, i)}, {P1, P2}), ...
                    1:size(x1, 2), 'UniformOutput', false);
        X = pflat(cat(2, X{:}));
    end

    % helper function to calculate
    % - relative orientation
    % - corresponding triangulation
    % - inlier indices
    function [bestP2n, bestX, inliers] = relative_orientation(x1n, x2n)
        % robustly estimate essential matrix
        [E, inliers] = estimate_E_robust(x1n, x2n, inlierThreshold);

        % compute normalised P1 and all possible P2s
        P1n = eye(3, 4);
        P2ns = extract_P_from_E(E);
        assert(numel(P2ns) == 4); % 4 cheirality configurations

        % for each P2n, triangulate 3D points to find P2n with correct cheirality
        bestX = [];
        bestP2n = 0;
        maxInFront = 0;
        for j = 1:4
            P2n = P2ns{j};
            X = triangulate(P1n, P2n, x1n(:, inliers), x2n(:, inliers));

            % count number of points in front of both cameras
            x1 = P1n * X;
            x2 = P2n * X;
            inFront = sum(x1(3, :) > 0 & x2(3, :) > 0);

            % update best values
            if inFront > maxInFront
                maxInFront = inFront;
                bestP2n = P2n;
                bestX = X;
            end
        end
    end

    % initialise cameras
    Pn = zeros(3, 4, length(imgNames));
    Pn(:, 1:3, 1) = eye(3); % R_1 = I

    % calculate relative orientations for each successive image pair
    for i = 1:length(imgNames)-1
        % load corresponding points
        x1n = correspondences{i}(:, :, 1);
        x2n = correspondences{i}(:, :, 2);

        % calculate relative orientation
        [P2n, ~, inliers] = relative_orientation(x1n, x2n);

        % update correspondences
        correspondences{i} = cat(3, x1n(:, inliers), x2n(:, inliers));
        
        % R_i+1         = R_i,i+1     * R_i
        Pn(:, 1:3, i+1) = P2n(:, 1:3) * Pn(:, 1:3, i);
    end

    % save(savefile, "Pn", "correspondences", "-append");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 3: Reconstruct X0 from initial image pair %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Reconstructing X0...");

    % % load corresponding points
    % precomputed = load(savefile);
    % correspondences = precomputed.correspondences;
    % descriptors = precomputed.descriptors;
    % matchesInitPair = precomputed.matchesInitPair;
    % Pn = precomputed.Pn;
    
    % normalise image points
    x1n = correspondences{end}(:, :, 1);
    x2n = correspondences{end}(:, :, 2);

    % calculate relative orientation and triangulation
    [~, X0, inliers] = relative_orientation(x1n, x2n);
    
    % transform X0 to world coordinates (up to a translation)
    XRot = ones(4, size(X0, 2));
    XRot(1:3, :) = Pn(:, 1:3, initPair(1))' * X0(1:3, :);
    
    % extract SIFT descriptors of inliers of i1 (used to compute X0)
    descX = descriptors{initPair(1)}(:, matchesInitPair(1, inliers));
    assert(all(size(XRot) == [4 length(inliers)]));

    % save(savefile, "XRot", "descX", "-append");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 4.1: Precompute correspondences for X0 %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Computing correspondences for X0...");

    % % load precomputed data
    % precomputed = load(savefile);
    % descriptors = precomputed.descriptors;
    % descX = precomputed.descX;

    matches2D3D = arrayfun(@(i) vl_ubcmatch(descriptors{i}, descX), ...
                           1:length(imgNames), 'UniformOutput', false);

    % save(savefile, "matches2D3D", "-append");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 4.2: Robustly calculate translations %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Robustly calculating translations...");

    % % load precomputed data
    % precomputed = load(savefile);
    % keypoints = precomputed.keypoints;
    % matches2D3D = precomputed.matches2D3D;
    % XRot = precomputed.XRot;
    % Pn = precomputed.Pn;

    % for each image i, calculate T_i
    for i = 1:length(imgNames)
        % extract corresponding images points and normalise
        x = keypoints{i}(:, matches2D3D{i}(1, :));
        xn = K \ [x; ones(1, size(x, 2))];

        % extract corresponding 3D points
        X = XRot(:, matches2D3D{i}(2, :));
        
        % update T_i
        R_i = Pn(:, 1:3, i);
        Pn(:,4,i) = estimate_T_robust(xn, X, R_i, 3 * inlierThreshold);
    end
    
    % save(savefile, "Pn", "-append");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 6: Triangulate and plot 3D points and cameras %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Triangulating and plotting 3D points and cameras...");

    % % load precomputed data
    % precomputed = load(savefile);
    % correspondences = precomputed.correspondences;
    % Pn = precomputed.Pn;

    % for each image pair, triangulate 3D points
    X = cell(1, length(imgNames)-1);
    for i = 1:length(imgNames)-1
        x1n = correspondences{i}(:, :, 1);
        x2n = correspondences{i}(:, :, 2);
        X{i} = triangulate(Pn(:,:,i), Pn(:,:,i+1), x1n, x2n);

        % discard points behind camera
        x1nProj = Pn(:,:,i) * X{i};
        x2nProj = Pn(:,:,i+1) * X{i};
        inFront = x1nProj(3, :) > 0 & x2nProj(3, :) > 0;
        X{i} = X{i}(:, inFront);
    end
    X = cat(2, X{:});

    % filter out points far from center of mass
    norms = vecnorm(X - mean(X, 2), 2);
    X = X(:, norms <= quantile(norms, 0.9));

    % unnormalise cameras
    Ps = arrayfun(@(i) K * Pn(:,:,i), ...
                  1:length(imgNames), 'UniformOutput', false);
    
    % plot 3D points and cameras
    figure;
    hold on;
    axis equal;
    plot3d(X, "MarkerSize", 2);
    plotcams(Ps);
    hold off;

    % save(['data/' num2str(dataset) '/results' num2str(dataset) '.mat'], "X", "Ps");
end