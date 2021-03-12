function [complex, commands, log] = WriteRadialCommandsFile(conddata,nd,h,dnmapid,savecond,saveq,saveeig,radstep)

% Check if dnmap id already exists
filenameshort = sprintf('dnmap/dnmap_%d*',dnmapid);
check = dir(filenameshort);
if ~isempty(check)
    error('Please choose unique dnmapid')
end

% Making folders
folders = {'dnmap_commands','dnmap_logs'};
for i = 1:length(folders)
    if ~exist(folders{i},'dir')
        mkdir(folders{i})
    end
end 

commands = sprintf('dnmap_commands/commands_dnmap_%d_radial_h_%.5f_radstep_%.5f_.dat',dnmapid,h,radstep);
log      = sprintf('dnmap_logs/log_dnmap_%d_radial_h_%.5f_radstep_%.5f_.dat',dnmapid,h,radstep);
commandsfile = fopen(commands,'w');
complex = 0;
logfile = fopen(log,'w');
fclose(logfile);

% Check error in reading
if commandsfile==-1
    error(['Error when opening file ',commandsfile]);
end

% nd check and write
fprintf(commandsfile,'#nd\n');
fprintf(commandsfile,'%d\n',nd);

% dnmap check and write
fprintf(commandsfile,'#dnmap\n');
fprintf(commandsfile,'dnmap_%d.dat\n',dnmapid);

% stepsize check and write
fprintf(commandsfile,'#stepsize\n');
fprintf(commandsfile,'%f\n',h);


% conductivity check and write
for nc = 1:length(conddata)
    if isa(conddata{nc},'RadialBallConductivity')
        fprintf(commandsfile,'#conductivity\n');
        fprintf(commandsfile,'##radius\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.radius);
        fprintf(commandsfile,'##amplitude\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.amplitude);
        fprintf(commandsfile,'##innerextent\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.innerextent);
    else
        error('Please use radially symmetric conductivities')
    end
end

% conductivity saving file check and write
if savecond ~= 0
    fprintf(commandsfile,'#conductivityonradius\n');
    % step size on radius
    fprintf(commandsfile,'%f\n',radstep);
    % conductivity
    fprintf(commandsfile,'conductivity_dnmapid_%d_true.dat\n',dnmapid);
end

% q saving file check and write
if saveq ~= 0
    fprintf(commandsfile,'#qonradius\n');
    % step size on radius
    fprintf(commandsfile,'%f\n',radstep);
    % q
    fprintf(commandsfile,'q_dnmapid_%d_true.dat\n',dnmapid);
end

% eigenvalues saving file check and write
if saveeig ~= 0
    fprintf(commandsfile,'#dnmapeigenvalues\n');
    % eigenvalues
    fprintf(commandsfile,'eig_dnmapid_%d_true.dat\n',dnmapid);
end

fclose(commandsfile);