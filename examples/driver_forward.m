
%% Forward problem solver %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%% Conductivity

% Possible Classes:
        
% Piecewise constant:
%       - PcBallConductivity ( center , radius , amplitude , n )
%       - PcEllipsoidConductivity ( center , radii , axes , amplitude , n )
% Radial:
%       - RadialBallConductivity ( radius , amplitude , innerextent )


%% Example: Hearts and lung phantom

% Orientation
rotang=5*pi/12;

% Prolate spheroid on the right
% Center
xr0=sin(rotang)*0.45; yr0=cos(rotang)*0.45; zr0=0;
% Radius
rr=0.42;
% Conductivity
cr=0.5;

% Prolate spheroid on the left
xl0=-sin(rotang)*0.55; yl0=cos(rotang)*0.55; zl0=0;
% Smallest radius 
rl=0.36;
% Conductivity
cl=0.5;

% Axes
A = [cos(rotang),sin(rotang),0;-sin(rotang),cos(rotang),0;0,0,1];

% Ball
% Center
xb0=-0.09; yb0=-0.55; zb0=0;
% Radius
rb=0.21;
% Conductivity
cb=2;

% Create different conductivity elements for pcc
SRpc = PcEllipsoidConductivity([xr0,yr0,zr0], [1.3*rr, 0.65*rr, 0.65*rr], A', cr, 10);
SLpc = PcEllipsoidConductivity([xl0,yl0,zl0], [1.3*rl, 0.65*rl, 0.65*rl], A, cl, 10);
B1pc = PcBallConductivity([xb0,yb0,zb0],1.3*rb,cb,10);

% Collect total conductivity data
conddata_pc = MakeConductivityData(SRpc,SLpc,B1pc);
%% Plot phantom
PlotPCPhantom(conddata_pc)
%% Forward map parameters

% Maximal degree of spherical harmonics in the
% computation of the Dirichlet-to-Neumann map.
nd = 10;  

% Mesh used to save the conductivity and q at points.
% Useful to compare the conductivity and reconstruction
% at points of a gmsh mesh.
mesh = 'ball_0p1_3D.msh';  % choose mesh from EITcode/mesh/ folder

savecond = 1;   % {0,1} - save conductivity data or not?
saveq    = 1;   % {0,1} - save q data or not?
parallel = 1;   % {0,1} - parallelize or not? (Requires openMP)

% (optional) Define your own script for running on a cluster
clusterscript = 0; 

% Unique id number for your forward computation
dnmapid  = 8;


%% Run

% Create command file
[complex, commands, log] = WritePcCommandsFile(conddata_pc,nd,dnmapid,...
                                               savecond,saveq,mesh);

% Start
StartDNMAP('pcc', parallel, complex, commands, log, clusterscript)

% Log
LogUpdate('dnmap', dnmapid)
