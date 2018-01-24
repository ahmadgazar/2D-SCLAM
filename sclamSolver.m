source "../tools/utilities/geometry_helpers_2d.m"
source "./indices.m"
source "./least_squares_landmarks.m"
source "./least_squares_poses.m"
source "./least_squares_calibration.m"

# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (2xnum_landmarks matrix of landmarks)
#   XK: the initial robot's paramters estimates (Kl, Kr, b)'
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation
#   XK: the paramters obtained by applying the perturbation

function [XR, XL, Xk] = boxPlus(XR, XL, Xk, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index = 1:num_poses)
    pose_matrix_index = poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr = dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index) = v2t(dxr) * XR(:,:,pose_index);
  endfor;
  for(landmark_index = 1:num_landmarks)
    landmark_matrix_index = landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl = dx(landmark_matrix_index:landmark_matrix_index + landmark_dim-1,:);
    XL(:,landmark_index)+= dxl;
  endfor;
  Xk+= dx(end-2:end);
endfunction;

#compute the measurement after obtaining the new estimate of the paramters
function [cl, cr, pose_pose_observations, pose_pose_associations] = compute_pose_pose_measurements(odometry_matrix, Xk, num_poses)

  num_pose_measurements = num_poses - 1;
  cl = zeros(1, num_pose_measurements);
  cr = zeros(1, num_pose_measurements);
  pose_pose_observations = zeros(3,3,num_pose_measurements);
  pose_pose_associations = zeros(2, num_pose_measurements);

  x_prev = 0;
  y_prev = 0;
  theta_prev = 0;
  kl = Xk(1);
  kr = Xk(2);
  b = Xk(3);
  measurement_num = 1;
  for pose_num = 1:num_pose_measurements

      Xi = v2t([x_prev y_prev theta_prev]);
      delta_cl = normalizeAngle(odometry_matrix(2,:,pose_num+1) - odometry_matrix(2,:,pose_num)); #incremental left encoder ticks
      delta_cr = normalizeAngle(odometry_matrix(3,:,pose_num+1) - odometry_matrix(3,:,pose_num)); #incremental right encoder tick

      cl(1, pose_num) = delta_cl;
      cr(1, pose_num) = delta_cr;

      theta = theta_prev + (delta_cr*kr - delta_cl*kl)/b;
      x = x_prev + cos(theta_prev) * (delta_cr*kr + delta_cl*kl)/2;
      y = y_prev + sin(theta_prev) * (delta_cr*kr + delta_cl*kl)/2;

      Xj = v2t([x  y theta]);
      pose_pose_associations(:, measurement_num) = [pose_num, pose_num+1]';
      pose_pose_observations(:,:,measurement_num) = inv(Xi) * Xj;
      measurement_num++;

      x_prev = x;
      y_prev = y;
      theta_prev = theta;
  endfor
endfunction


# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   XK: the initial robot's paramters estimates (Kl, Kr, b)'
#   odometry_matrix:
#   3x1x2507 matrix (:,1,k)= [pose_id, encoder left ticks, encoder right ticks]
#   Zl:  the measurements (2xnum_measurements)
#   landmark associations: 2xnum_measurements.
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   Xk: the parameters after optimization
#   chi_stats_{l,poses, parameters}: array 1:num_iterations, containing evolution of chi2 for landmarks, poses and parameters
#   num_inliers{l,poses, parameters}: array 1:num_iterations, containing evolution of inliers landmarks, poses and parameters


function [XR, XL, Xk, chi_stats_l, num_inliers_l, chi_stats_poses, num_inliers_poses, chi_stats_parameters, num_inliers_parameters, H, b] = doSCLAM(XR, XL, Xk,
       odometry_matrix,
	     Zl, landmark_associations,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold)

  global pose_dim;
  global landmark_dim;
  global k_dim;

  chi_stats_l = zeros(1,num_iterations);
  num_inliers_l = zeros(1,num_iterations);

  chi_stats_poses = zeros(1,num_iterations);
  num_inliers_poses = zeros(1,num_iterations);

  chi_stats_parameters = zeros(1,num_iterations);
  num_inliers_parameters= zeros(1,num_iterations);

  # size of the linear system
  system_size = pose_dim*num_poses + landmark_dim*num_landmarks + k_dim;
  for (iteration = 1:num_iterations)
    H = zeros(system_size, system_size);
    b = zeros(system_size, 1);

    [H_landmarks, b_landmarks, chi_, num_inliers_] = linearizeLandmarks(XR, XL, Zl, landmark_associations, num_poses, num_landmarks, kernel_threshold);
    chi_stats_l(iteration)+= chi_;
    num_inliers_l(iteration) = num_inliers_;

    [cl, cr, Zr, poses_associations] = compute_pose_pose_measurements(odometry_matrix, Xk, num_poses);
    [H_poses, b_poses, Ji, Jj, chi_, num_inliers_] = linearizePoses(XR, Zr, poses_associations, num_poses, num_landmarks, kernel_threshold);
    chi_stats_poses(iteration)+= chi_;
    num_inliers_poses(iteration) = num_inliers_;

    [H_parameters, b_paramters, chi_, num_inliers_] = linearizeParameters(Ji, Jj, XR, Xk, cl, cr, Zr, poses_associations, num_poses, num_landmarks, kernel_threshold);
    chi_stats_parameters(iteration)+= chi_;
    num_inliers_parameters(iteration) = num_inliers_;

    H+= H_landmarks + H_poses + H_parameters;
    b+= b_landmarks + b_poses + b_paramters;

    H+= eye(system_size) * damping;
    dx= zeros(system_size, 1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end) = -(H(pose_dim+1:end,pose_dim+1:end) \ b(pose_dim+1:end,1));
    [XR, XL, Xk] = boxPlus(XR, XL, Xk,  num_poses, num_landmarks, dx);

  endfor
endfunction
