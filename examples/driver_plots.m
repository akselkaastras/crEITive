%%  Plots  %%
%%%%%%%%%%%%%
clc; clear all;

%% Plot 2D

% You may choose from {'conductivity','q','qhat'}
type    = 'conductivity'; 

reconid = [1,2]; % reconstruction id
dnmapid = 0;     % include true data by giving dnmapid

v       = [0,0,1];    % normal vector to slice
h       = 0.01;       % stepsize
span    = [0.65,1.7]; % span of type

[ha, pos] = Plot2D(type, reconid, dnmapid, v, h, span, 0,'real');


%% Plot 1D

% You may choose from {'conductivity','q','qhat'}
type         = 'conductivity';

% comp_reconid is a matrix with columns representing 
% plots, in which the elements are ids to types we 
% want to compare in 1D plot
comp_reconid = cat(2,[1;2],[1;2]);

% include true data by giving dnmapids in a vector
include_true = [1,1];

% comp_what is a cell list length == size(comp_reconid,2) 
% of "comparison"-strings. One "comparison"-string in 
% {'method','ift','ngrid','fixed','pkappa'} for each
% column in comp_reconid.
comp_what    = {'reconmethod','reconmethod'};

% Restriction of the conductivity to a line requires a
% direction v and a point center
v            = [1,0,0];
center       = [0,0,0];

h            = 0.01;      % stepsize
span         = [0.4,1.3]; % span of type

[ha, pos] = Plot1D(type,comp_reconid,include_true,comp_what,v,center,h,span);

%% Plot 3D

% 3D plots supports only conductivity type of data
reconid = 1; % reconstruction id
dnmapid = 0; % include true data by giving dnmapid

% Slices are given as for the Matlab slice function
xslice  = [-0.6,-0.05,0.6]; % slice vector x
yslice  = [];               % slice vector y
zslice  = 0;                % slice vector z

h       = 0.01;      % stepsize
span    = [0.8,1.4]; % span of type

[ha, pos] = Plot3Dc(reconid, dnmapid, xslice, yslice, zslice, h, span);
