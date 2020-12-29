function [s,R,T] = estimateSimTrans(X1,X2)
% Estimate Similarity Transformation
%
% [s,R,T] = estimateSimTrans(X1,X2)
% 
% X1,X2 ... 3xN matrices with corresponding 3D points
%
% X2 = s*R*X1 + T
% s ... scalar scale
% R ... 3x3 rotation matrix
% T ... 3x1 translation vector
%
% This is done according to the paper:
% "Least-Squares Fitting of Two 3-D Point Sets"
% by K.S. Arun, T. S. Huang and S. D. Blostein
%
% $Id: estsimt.m,v 2.1 2005/05/23 16:24:59 svoboda Exp $
%
% updated 01/17/2015 by John D. Long II, PhD, contact: jlong29@gmail.com

% number of points
N = size(X1,2);

if N ~= size(X2,2)
  error('estimateSimTrans: both sets must contain the same number of points')
end

X1cent = mean(X1,2);
X2cent = mean(X2,2);
% normalize coordinate systems for both set of points
x1 = X1 - repmat(X1cent,1,N);
x2 = X2 - repmat(X2cent,1,N);

% first estimate the scale s by comparing each point's distance from origin
dists1 = sum(sqrt(x1.^2));
dists2 = sum(sqrt(x2.^2));

% Estimating scales according to ratio of paired distances from the origin
scales0 = dists2./dists1;

% Using distances between point pairs; changed by JDL2 on 10/1/2013
C  = combnk(1:N,2);
d1 = zeros(3,size(C,1));
d2 = zeros(3,size(C,1));
for ii = 1:size(C,1)
    d1(:,ii) = x1(:,C(ii,1))-x1(:,C(ii,2));
    d2(:,ii) = x2(:,C(ii,1))-x2(:,C(ii,2));
end
ds1 = sum(d1.^2).^(1/2);
ds2 = sum(d2.^2).^(1/2);

scales1 = ds2./ds1;

s	   = median([scales0 scales1]);
% undo scale
x1s = s*x1;

% finding rotation
H = zeros(3,3);
for i=1:N
  H = H + x1s(:,i)*x2(:,i)';
end

[U,S,V] = svd(H,0);
X = V*U';
if det(X) > 0
  R = X; % estimated rotation matrix
else
  % V = [V(:,1:2),-V(:,3)];
  % R = V*U';
  % s = -s;
  R = X;
end

T = X2cent - s*R*X1cent;
