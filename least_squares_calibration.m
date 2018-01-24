source "../tools/utilities/geometry_helpers_2d.m"
source "./indices.m"


function [e, Jk] = parametersErrorAndJacobian(Xi, Xj, Xk, cr, cl, Z)

  kl = Xk(1);
  kr = Xk(2);
  b  = Xk(3);

theta = normalizeAngle((Xj)(3) - (Xi)(3));
J_chordal = [1 	0   -0.5*(cr*kr+cl*kl)*sin(theta);
      	     0	1    0.5*(cr*kr+cl*kl)*cos(theta);
      	     0	0   -sin(theta);
      	     0	0    cos(theta);
      	     0  0    -cos(theta);
      	     0  0    -sin(theta)];

J_model = [0.5*cl*cos(theta)    0.5*cr*cos(theta)   0;
           0.5*cl*sin(theta)    0.5*cr*sin(theta)   0;
             -cl/b                 cr/b              (cl*kl - cr*kr)/(b^2)];

Jk = J_chordal * J_model;
Z_hat = inv(Xi) * Xj;
e = flattenIsometryByColumns(Z - Z_hat);
endfunction


function [H,b, chi_tot, num_inliers] = linearizeParameters(Ji, Jj, XR, Xk, cl, cr, Zr, pose_associations, num_poses, num_landmarks, kernel_threshold)
  global pose_dim;
  global landmark_dim;
  global k_dim;

  system_size = pose_dim*num_poses + landmark_dim*num_landmarks + k_dim;
  H = zeros(system_size, system_size);
  b = zeros(system_size, 1);
  chi_tot = 0;
  num_inliers = 0;

  for (measurement_num = 1:size(Zr,3))
    Omega = eye(6);
    #Omega(1:3,1:3)*=1e3; # we need to pimp the rotation  part a little
    pose_i_index = pose_associations(1,measurement_num);
    pose_j_index = pose_associations(2,measurement_num);
    Z = Zr(:,:,measurement_num);
    Xi = XR(:,:,pose_i_index);
    Xj = XR(:,:,pose_j_index);

    [e,Jk] = parametersErrorAndJacobian(Xi, Xj, Xk, cr(1, measurement_num), cl(1, measurement_num), Z);
    chi = e'*Omega*e;
    if (chi > kernel_threshold)
      Omega*= sqrt(kernel_threshold/chi);
      chi = kernel_threshold;
    else
      num_inliers ++;
    endif;
    chi_tot+=chi;


    pose_i_matrix_index = poseMatrixIndex(pose_i_index, num_poses, num_landmarks);
    pose_j_matrix_index = poseMatrixIndex(pose_j_index, num_poses, num_landmarks);

    H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1, end-2:end)+= Ji'*Omega*Jk;
    H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1, end-2:end)+= Jj'*Omega*Jk;

    H(end-2:end, pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+= Jk'*Omega*Ji;
    H(end-2:end, pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+= Jk'*Omega*Jj;
    H(end-2:end, end-2:end)+= Jk'*Omega*Jk;
    b(end-2:end)+= Jk'*Omega*e;

  endfor
endfunction
