close all
clear
clc

addpath '../'
addpath '../tools/g2o_wrapper'
addpath '../tools/visualization'
source "../tools/utilities/geometry_helpers_2d.m"
source "./sclamSolver.m"


kl = 0.1;
kr = 0.1;
b = 0.4;

#load sclam dataset dataset
[landmarks, poses, transitions, observations, odometry] = loadG2o('dataset-sclam.g2o');

#useful dimensions
num_poses = size(poses,2);
num_landmarks = size(landmarks, 2);

#data manipulation

##################### TRUE LANDMARKS############################################
landmarks_matrix = cell2mat(struct2cell(landmarks)); #convert to 2x1x25 matrix
XL_true = zeros(2, num_landmarks);
landmarks_matrix = landmarks_matrix(2:3,:,:);
for i = 1:num_landmarks
  XL_true(:,i) = landmarks_matrix(:,:,i);
endfor

##################### TRUE POSES################################################
poses_matrix = cell2mat(struct2cell(poses));
poses_matrix = poses_matrix(2:4,:,:);#convert to 3x1x2507 matrix
XR_true = zeros(3,3, num_poses);
for i = 1:num_poses
  v = poses_matrix(:,:,i)';
  XR_true(:,:,i) = v2t(v);
endfor


####################### REAL ODOMETRY MEASUREMENTS##############################
odometry_matrix = cell2mat(struct2cell(odometry)); #convert to 3x1x2507 matrix

######################## REAL LANDMARK MEASUREMENTS ############################

disp('Extracting landmarks');

num_landmark_measurements = 0;
pose_id_vec = [];
landmark_id_vec = [];
landmark_x_observation = [];
landmark_y_observation = [];

for i = 1:size(observations, 2)
    dim_curr_obs = size(observations(i).observation, 2);
    pose_id = i; #observations(i).pose_id - 1100; #extract pose id from measurements
    pose_id_vec(1, end+1:end+dim_curr_obs) = pose_id;
  for j = 1:dim_curr_obs
      landmark_id = observations(i).observation(j).id;
      landmark_id_vec(1, end+1) = landmark_id;
      landmark_x_observation(1, end+1) = observations(i).observation(j).x_pose;
      landmark_y_observation(1, end+1) = observations(i).observation(j).y_pose;
  endfor
endfor
pose_landmark_associations = [pose_id_vec; landmark_id_vec];
landmark_observations = [landmark_x_observation; landmark_y_observation];


################# GENERATION OF (WRONG) INITIAL GUESS ##########################

# apply a perturbation to each ideal pose (construct the estimation problem)
disp('Applying perturbations to ideal poses and landmarks!');
pert_deviation = 1;
pert_scale = eye(3)*pert_deviation;
XR_guess = XR_true;
XL_guess = XL_true;

for (pose_num = 2:num_poses)
    xr = rand(3,1) - 0.5;
    dXr = v2t(pert_scale * xr);
    XR_guess(:,:,pose_num) = dXr * XR_guess(:,:,pose_num);
endfor

#apply a perturbation to each landmark
dXl = (rand(landmark_dim, num_landmarks)-0.5) * pert_deviation;
XL_guess+= dXl;

#robot paramters initial guess
Xk_guess = [kl;kr;b];

########################### CALL TRUE SCLAM SOLVER #############################

disp('HERE WE FUCKING GO !');
damping = 0.0001;
kernel_threshold = 1e3;
num_iterations = 10;
 [XR, XL, Xk, chi_stats_l, num_inliers_l, chi_stats_poses, num_inliers_poses, chi_stats_parameters, num_inliers_parameters, H, b] = doSCLAM(XR_guess, XL_guess, Xk_guess,
                              odometry_matrix,
												      landmark_observations, pose_landmark_associations,
												      num_poses,
												      num_landmarks,
												      num_iterations,
												      damping, kernel_threshold);

####################### DISPLAY DATA ###########################################
disp('DISPLAY RESULTS!');

figure(1);
hold on;
grid;

subplot(2,2,1);
title("Landmark Initial Guess");
plot(XL_true(1,:), XL_true(2,:),'b*',"linewidth",2);
hold on;
plot(XL_guess(1,:), XL_guess(2,:),'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;

subplot(2,2,2);
title("Landmark After Optimization");
plot(XL_true(1,:), XL_true(2,:), 'b*',"linewidth",2);
hold on;
plot(XL(1,:), XL(2,:), 'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;

subplot(2,2,3);
title("Poses Initial Guess");
plot(XR_true(1,3,:),XR_true(2,3,:),'b*',"linewidth",2);
hold on;
plot(XR_guess(1,3,:),XR_guess(2,3,:),'ro',"linewidth",2);
legend("Poses True", "Guess");grid;

subplot(2,2,4);
title("Poses After Optimization");
plot(XR_true(1,3,:), XR_true(2,3,:), 'b*',"linewidth",2);
hold on;
plot(XR(1,3,:), XR(2,3,:), 'ro',"linewidth",2);
legend("Poses True", "Guess");grid;


figure(2);
hold on;
grid;
title("chi evolution");


subplot(3,2,1);
plot(chi_stats_parameters, 'r-', "linewidth", 2);
legend("Chi parameters"); grid; xlabel("iterations");
subplot(3,2,2);
plot(num_inliers_parameters, 'b-', "linewidth", 2);
legend("#inliers"); grid; xlabel("iterations");


subplot(3,2,3);
plot(chi_stats_l, 'r-', "linewidth", 2);
legend("Chi Landmark"); grid; xlabel("iterations");
subplot(3,2,4);
plot(num_inliers_l, 'b-', "linewidth", 2);
legend("#inliers"); grid; xlabel("iterations");

subplot(3,2,5);
plot(chi_stats_poses, 'r-', "linewidth", 2);
legend("Chi Poses"); grid; xlabel("iterations");

subplot(3,2,6);
plot(num_inliers_poses, 'b-', "linewidth", 2);
legend("#inliers");grid; xlabel("iterations");
