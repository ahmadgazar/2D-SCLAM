source "./indices.m"
source "../tools/utilities/geometry_helpers_2d.m"

# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose in world frame (3x3 homogeneous matrix)
#   Xl: the landmark pose (2x1 vector, 3d pose in world frame)
#   z:  measured position of landmark
# output:
#   e: 2x1 is the difference between prediction and measurement
#   Jr: 2x3 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 2x2 derivative w.r.t a the error and a perturbation on the
#       landmark

function [e,Jr,Jl] = landmarkErrorAndJacobian(Xr,Xl,z)
  # inverse transform
  iR = Xr(1:2,1:2)';
  it = -iR * Xr(1:2,3);
  #t = Xr(1:2,3);

  #prediction
  z_hat = iR*Xl + it;
  e = z_hat - z;
  Jr = zeros(2,3);
  Jr(1:2,1:2) = -iR;
  Jr(1:2,3) = iR * [Xl(2) -Xl(1)]';
  Jl = iR;

endfunction;

#linearizes the robot-landmark measurements
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   Zl:  the measurements (2xnum_measurements)
#   associations: 2xnum_measurements.
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   kernel_threshod: robust kernel threshold
# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers

function [H,b, chi_tot, num_inliers] = linearizeLandmarks(XR, XL, Zl, associations, num_poses, num_landmarks, kernel_threshold)
  global pose_dim;
  global landmark_dim;
  global k_dim;

  system_size = pose_dim*num_poses + landmark_dim*num_landmarks + k_dim;
  H = zeros(system_size, system_size);
  b = zeros(system_size, 1);
  chi_tot = 0;
  num_inliers = 0;

  for (measurement_num = 1:size(Zl,2))
    pose_index = associations(1, measurement_num);
    landmark_index = associations(2, measurement_num);
    z = Zl(:,measurement_num);
    Xr = XR(:,:,pose_index);
    Xl = XL(:,landmark_index);


    [e, Jr, Jl] = landmarkErrorAndJacobian(Xr, Xl, z);
    chi = e'*e;
    if (chi > kernel_threshold)
      e*= sqrt(kernel_threshold/chi);
      chi = kernel_threshold;
    else
      num_inliers++;
    endif;
    chi_tot+= chi;

    pose_matrix_index = poseMatrixIndex(pose_index, num_poses, num_landmarks);
    landmark_matrix_index = landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

    H(pose_matrix_index:pose_matrix_index + pose_dim-1,
      pose_matrix_index:pose_matrix_index + pose_dim-1)+= Jr'*Jr;

    H(pose_matrix_index:pose_matrix_index + pose_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+= Jr'*Jl;

    H(landmark_matrix_index:landmark_matrix_index + landmark_dim-1,
      landmark_matrix_index:landmark_matrix_index + landmark_dim-1)+= Jl'*Jl;

    H(landmark_matrix_index:landmark_matrix_index + landmark_dim-1,
      pose_matrix_index:pose_matrix_index + pose_dim-1)+= Jl'*Jr;

    b(pose_matrix_index:pose_matrix_index + pose_dim-1)+= Jr'*e;
    b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+= Jl'*e;
  endfor
endfunction
