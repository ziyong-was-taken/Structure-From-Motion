# 3D Reconstruction

<!-- omit in toc -->
## Table of Contents

- [Usage](#usage)
- [High-Level Overview](#high-level-overview)

## Usage

- set up [VLFeat](https://www.vlfeat.org)
- run the MATLAB function [`code/run_sfm.m`](./code/run_sfm.m)

## High-Level Overview

0. precompute SIFT features and descriptors for all images
1. for each successive image pair ($i$, $i+1$)
   1. match descriptors to get correspondences (normalised matches)
   2. use RANSAC to robustly estimate $E$ using 8-point algorithm
   3. extract all possible $[R\quad T]$ (4 solutions) from $E$
   4. cheirality check: triangulate points using each solution and choose the solution with the most points in front of it
2. upgrade relative rotations to absolute rotations by chaining\
   $R_1 = I$ and $R_{i+1} = R_{i,i+1}R_i$
3. for initial image pair ($i_1$ and $i_2$)
   1. repeat step 1
      1. save best 3D points as $\mathcal X_0$
      2. save descriptors `descX` of inliers of $i_1$ used to triangulate $\mathcal X_0$
   2. rotate $\mathcal X_0$ to (partial) world coordinates $X_{rot} = X_0 - C_{i_1} = R_{i_1}^\top \mathcal X_0$ where camera center $C_{i_1} = -R_{i_1}^\top T_{i_1}$ (unknown since $T_{i_1}$ unknown)
4. use RANSAC to robustly estimate $T_i$ for each image $i$
   1. match descriptors `descriptors{i}` and `descX` to find correspondences (points in image $i$ that correspond to $\mathcal X_0$)
   2. use DLT with 2 point correspondences to estimate $T_i$ (rewrite resection problem with $R$ known)
   3. use reprojection error $\|x_{Meas} - x_{Proj}\|_F^2$ to evaluate $T_i$
5. for each successive image pair ($i$, $i+1$):\
   triangulate points $X_{i,i+1}$ using $[R_i\quad T_i]$ and $[R_{i+1}\quad T_{i+1}]$
6. plot all $X_{i,i+1}$ and $P_i$

For more details, see the [report](./Report.pdf).
