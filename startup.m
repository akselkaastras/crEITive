% startup.m

clc; clear all;

mydir = pwd;
id  = strfind(mydir,filesep);
folder = mydir(id(end)+1:end);

% navigating to /crEITive/
if strcmpi(folder,'crEITive')
    cd Matlab/interface
    addpath classes
    cd ..
    addpath interface
    addpath interface/cluster
    addpath util
    cd ..
    addpath Matlab
    addpath commands/dnmap_commands
    addpath logs/dnmap_logs
    addpath commands/eit_commands
    addpath logs/eit_logs
    addpath examples
elseif strcmpi(folder,'eitcode')
else
    error('Please, navigate to /crEITive');
end
clear all;