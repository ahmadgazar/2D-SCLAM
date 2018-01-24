source "../tools/utilities/geometry_helpers_2d.m"
source "./indices.m"
source "./least_squares_calibration.m"
source "./compute_pose_pose_measurements.m"
# computes the odometry transition matrix between 2 robot poses
# input:
#   odometry_matrix: 3x1x2507 [:,:,k] = [pose_id left_ticks right_ticks]'
#   Xk: the robot paramters vector [kl kr b]'
# output:
#   cl: 1x2506 vector consisting of left encoder ticks
#   cr: 1x2506 vector consisting of right encoder ticks
#   Zr: 3x3x2506 tensor
#       each 3x3 matrix is the transition between 2 robot poses
#      based on the incremental encoder ticks
#   pose_associations: 2xnum_measurements.
#                 associations(:,k)=[i_idx, j_idx]' means the kth measurement
#                 refers to an observation made from pose i_idx, that
#                 observed the pose j_idx

%{
function [cl, cr, Zr, pose_associations] = computeOdometry(odometry_matrix, Xk, num_poses)

  num_pose_measurements = num_poses - 1;
  Zr = zeros(3,3,num_pose_measurements);
  pose_associations = zeros(2, num_pose_measurements);
  x_prev = 0;
  y_prev = 0;
  theta_prev = 0;

  kl = Xk(1);
  kr = Xk(2);
  b = Xk(3);

  cl = zeros(1,num_pose_measurements);
  cr = zeros(1,num_pose_measurements);

  measurement_num = 1;
  for pose_num = 1:(num_poses - 1)
      Xi = v2t([x_prev y_prev theta_prev]);
      delta_cl = normalizeAngle(odometry_matrix(2,:,pose_num+1) - odometry_matrix(2,:,pose_num)); #incremental left encoder ticks
      delta_cr = normalizeAngle(odometry_matrix(3,:,pose_num+1) - odometry_matrix(3,:,pose_num)); #incremental right encoder ticks

      cr(measurement_num) = delta_cr;
      cl(measurement_num) = delta_cl;

      theta = theta_prev + (delta_cr*kr - delta_cl*kl)/b;
      x = x_prev + cos(theta_prev) * (delta_cr*kr + delta_cl*kl)/2;
      y = y_prev + sin(theta_prev) * (delta_cr*kr + delta_cl*kl)/2;
      Xj = v2t([x  y theta]);
      pose_associations(:, measurement_num) = [pose_num, pose_num+1]';
      Zr(:,:,measurement_num) = inv(Xi) * Xj;
      measurement_num++;
      x_prev = x;
      y_prev = y;
      theta_prev = theta;
  endfor
endfunction
%}

# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   Z:  the measurements (2xnum_measurements)
#   associations: 2xnum_measurements.
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
#   chi_stats_{l,k}: array 1:num_iterations, containing evolution of chi2 for landmarks and platform parameters
#   num_inliers{l,k}: array 1:num_iterations, containing evolution of inliers landmarks and platform

function [XR, Xk, chi_stats_parameters, num_inliers_parameters, H, b] = doCalibration(XR, Xk,
       odometry_matrix,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold)


  global k_dim;

  chi_stats_parameters = zeros(1,num_iterations);
  num_inliers_parameters= zeros(1,num_iterations);

  # size of the linear system
  system_size = k_dim;
  for (iteration = 1:num_iterations)
    H = zeros(system_size, system_size);
    b = zeros(system_size, 1);

    [cl, cr, Zr, poses_associations] = compute_pose_pose_measurements(odometry_matrix, Xk, num_poses);
    [H_parameters, b_paramters, chi_, num_inliers_] = linearizeParameters(XR, Xk, cl, cr, Zr, poses_associations, num_poses, num_landmarks, kernel_threshold);
    chi_stats_parameters(iteration) = chi_;
    num_inliers_parameters(iteration) = num_inliers_;

    H+= H_parameters;
    b+= b_paramters;

    H+= eye(system_size) * damping;
    dx= zeros(system_size, 1);

    dx = -H \ b;
    Xk+= dx

  endfor
endfunction
