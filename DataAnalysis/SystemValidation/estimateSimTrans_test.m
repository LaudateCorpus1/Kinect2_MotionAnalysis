% TEST script for similarity transform estimation method
clear
clc

% functions for composing a 3D rotation
Rotx = @(theta) [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Roty = @(theta) [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rotz = @(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

homogenize   = @(x) [x;ones(1,size(x,2))];
dehomogenize = @(x) x(1:end-1,:);

% Initial Data
X1 = randn(3,5);

% pick random rotations uniformly [-pi, pi]
thetaX = -pi + 2*pi*rand;
thetaY = -pi + 2*pi*rand;
thetaZ = -pi + 2*pi*rand;

% Construction random similarity Transform
Rotxyz = Rotx(thetaX)*Roty(thetaY)*Rotz(thetaZ);
Scale  = 0.1 + 99.9*rand;   % [0.1, 100]
Trans  = -10 + 20*rand(3,1);% [-10, 10] in 3D

SimT   = [Scale*Rotxyz Trans; 0 0 0 1];

% Transform X1 into X2
X2     = dehomogenize(SimT*homogenize(X1));

% Estimate Similarity Transform
[Scale_hat,Rotxyz_hat,Trans_hat] = estimateSimTrans(X1,X2);

% Estimated similarity transformation
SimT_hat = [Scale_hat*Rotxyz_hat Trans_hat];
    
% Check errors
fprintf(1,'Scale error          = %6.4e\n',Scale-Scale_hat);
fprintf(1,'Rotation error       = %6.4e\n',norm(Rotxyz-Rotxyz_hat));
fprintf(1,'Translation error    = %6.4e\n',norm(Trans-Trans_hat));
fprintf(1,'Transformation error = %6.4e\n',norm(SimT(1:3,:)-SimT_hat));
