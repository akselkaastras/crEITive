%%  Reconstruction demo  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all;
%% Parameters

% Method used for the reconstruction
%   - calderon for a reconstruction with Calderón approximation.
%   - texp for a reconstruction with texp approximation.
%   - t0 for a reconstruction with t0 approximation.
%   - t for a reconstruction with full algorithm.
reconmethod = 'texp';

% Computation of inverse Fourier transform
%   - ifft for an Inverse Fast Fourier Transform
ift = 'ifft';

% Resolution of grid where the scattering transform is computed
ngrid = 11;

% Truncation radius (maximum frequency scattering transform is computed on)
% We must have truncrad <= pi*(ngrid-1)^2/(2*ngrid) by Shannon Sampling
truncrad = 9;

% Complex frequency (zeta) size and type.
%   - |zeta| = pkappa*truncrad/sqrt(2) if fixed = 1.
%   - |zeta| = pkappa*|xi|/sqrt(2) for each xi if fixed = 0.
fixed       = 1; % {0,1} choose zeta fixed (1) or proportional (0)
pkappa      = 1; % |zeta| = p

% Compute conductivity and q on mesh
mesh = 'ball_0p05_3D.msh';  % choose mesh from crEITive/mesh/ folder

dnmapid     = 10; % unique dnmap integer id number of dnmap data
reconid     = 10; % unique reconstruction integer id number
parallel    = 1; % {0,1} - parallelize or not?
cluster     = 0;

%% Run

% Commands
[complex, commands, log] = WriteCommandsFile(reconmethod,ift,ngrid,truncrad,...
                            fixed,pkappa,mesh,reconid,dnmapid);
% Start
StartEIT(parallel, complex, commands, log,cluster);

% Log 
LogUpdate('eit',reconid)


%% Run full algorithm

% Full algorithm
reconmethod = 't';
reconid     = 2;

% Commands
[complex, commands, log] = WriteCommandsFile(reconmethod,ift,ngrid,truncrad,...
                            fixed,pkappa,mesh,reconid,dnmapid);
% Start
StartEIT(parallel, complex, commands, log,cluster);

% Log 
LogUpdate('eit',reconid)
