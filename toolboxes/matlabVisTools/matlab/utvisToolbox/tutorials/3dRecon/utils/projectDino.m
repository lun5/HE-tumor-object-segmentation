function [p, polys, Mint, Mext] = projectDino(f, d, R, sclZ)
  %% [p, f, Mint, Mext] = projectDino(f, d, R, sclZ)
  %%
  %% Return the image positions of points on the Dino model
  %% as viewed from a camera with focal length f, and nodal point at d.  
  %% Optional args:  
  %%  R is the 3x3 rotation matrix, or 
  %     if R =[] or R is omitted, 
  %%    then camera is rotated so that the optical axis passes through the
  %%    mean of Dino's 3D point data.  
  %%  sclZ,  Dino is scaled in its Z direction by
  %%    the factor sclZ (default 1) before projection.
  %% On return, polys is a descriptor for individual polygons in the dino
  %% model.  This is used in showWire.
  %% Also the intrinsic and extrinsic calibration
  %% matrices, Mint and Mext = [R, -R*d], respectively.
  %%
  %% This routine could generate divide by zero warnings depending on the
  %% placement of the camera.
    
if nargin < 3 | size(R,1) == 0
  computeR = 1;
  R = [];
else
  computeR = 0;
end
if nargin < 4
  sclZ = 1;
end

DEBUG = 0;

%% Read  Dino data set
[v polys] = getHalfDino;

%%% Set canonical 3D pose
R0 = [1 0 0; 0 0 -1; 0 -1 0];   %% Rotate and reflect dino (concave away).
mn0 = [0 0 0]';                %% Place Dino's mean at origin
P0 = R0 * v' + repmat(mn0, 1, size(v,1));
if sclZ ~= 1.0
  P0(3,:) = sclZ * P0(3,:);
end

if DEBUG
  figure(1);
  showWire(P0', polys);
  xlabel('Z');ylabel('Z');zlabel('Z');
  title('3D Dino Model');
  pause;
end
%% Build 3D homogeneous coordinates.
P0 = [P0; ones(1,size(P0,2))];

if computeR
  %% Choose rotation to fixate mn0.
  %% That is solve:  R * (mn0 - d) = [0 0 1]';
  R = eye(3);
  R(:,3) = (mn0 - d)/norm(mn0 - d);
  R(:,2) = R(:,2) - (R(:,3)' * R(:,2)) * R(:,3);
  R(:,2) = R(:,2)/norm(R(:,2));
  R(:,1) = R(:,1) - R(:,2:3) * (R(:,2:3)' * R(:,1));
  R(:,1) = R(:,1)/norm(R(:,1));
  R = R';
end

%% Build intrinsic and extrinsic calibration matrices Mint and Mext.
Mint = diag([f, f, 1]);
Mext = [R, -R*d];

%% Construct projection matrix
M = Mint * Mext;

%% Compute the homogeneous image locations
p = M * P0;

%% Convert to image pixels
p = p ./repmat(p(3,:), 3,1);
