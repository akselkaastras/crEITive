
%% Forward problem solver %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;


%% Example: Head&Skull phantom

% Ellipsoid: skull
% Axes
A = eye(3);
% Radii
rs = [0.9,0.95,0.95];
% Conductivity
cs = 0.2;
% Width
width = 0.04;

% Ball: Hemorrhagic stroke
% Center
x = [0.5,0.4,0];
% Radius
rb = 0.15;
% Conductivity
cb = 4;

N = 25;

% Create different conductivity elements for pcc
Spc = PcEllipsoidShell(x,rs,A,cs,width,N);
Bpc = PcBallConductivity(x,rb,cb,N);

% Collect total conductivity data
conddata_pc = MakeConductivityData(Spc,Bpc);

%% Forward map parameters

% Maximal degree of spherical harmonics in the
% computation of the Dirichlet-to-Neumann map.
nd = 25;  

% Mesh used to save the conductivity and q at points.
% Useful to compare the conductivity and reconstruction
% at points of a gmsh mesh.
mesh = 'ball_0p05_3D.msh';  % choose mesh from EITcode/mesh/ folder

savecond = 1;   % {0,1} - save conductivity data or not?
saveq    = 1;   % {0,1} - save q data or not?
parallel = 1;   % {0,1} - parallelize or not? (Requires openMP)

% (optional) Define your own script for running on a cluster
clusterscript = 0; 

% Unique id number for your forward computation
dnmapid  = 11;


%% Run

% Create command file
[complex, commands, log] = WritePcCommandsFile(conddata_pc,nd,dnmapid,...
                                               savecond,saveq,mesh);

% Start
StartDNMAP('pcc', parallel, complex, commands, log, clusterscript)

% Log
LogUpdate('dnmap', dnmapid)

%% plot

% You may choose from {'conductivity','q','qhat'}
type    = 'conductivity'; 

dnmapid = 11;     % include true data by giving dnmapid

v       = [0,0,1];    % normal vector to slice
h       = 0.01;       % stepsize
span    = [0.5,2]; % span of type

[ha, pos] = Plot2D(type, 10, dnmapid, v, h, span, 0,'real');