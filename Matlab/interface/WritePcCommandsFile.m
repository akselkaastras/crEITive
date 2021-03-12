function [complex, commands, log] = WritePcCommandsFile(conddata,nd,dnmapid,savecond,saveq,mesh)

% Check if dnmap id already exists
filenameshort = sprintf('results/dnmap/dnmap_%d*',dnmapid);
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

commands = sprintf('commands/dnmap_commands/commands_dnmap_%d_pcc_nd_%d_%s.dat',dnmapid,nd,mesh(1:(end-4)));
log      = sprintf('logs/dnmap_logs/log_dnmap_%d_pcc_nd_%d_%s.dat',dnmapid,nd,mesh(1:(end-4)));
commandsfile =fopen(commands,'w');
logfile = fopen(log,'w');
fclose(logfile);

% Check error in reading
if commandsfile==-1
    error(['Error when opening file ',commandsfile]);
end

% Complex conductivity?
complex = 0;
for nc = 1: length(conddata)
    if ~isreal(conddata{nc}.amplitude)
        complex = 1;
        break;
    end
end

% nd check and write
fprintf(commandsfile,'#nd\n');
fprintf(commandsfile,'%d\n',nd);

% dnmap check and write
fprintf(commandsfile,'#dnmap\n');
fprintf(commandsfile,'dnmap_%d.dat\n',dnmapid);

% conductivity check and write
for nc = 1:length(conddata)
    fprintf(commandsfile,'#conductivity\n');
    fprintf(commandsfile,'##center\n');
    fprintf(commandsfile,'cartesian\n');
    fprintf(commandsfile,'%f ', conddata{nc}.center);
    fprintf(commandsfile,'\n');

    % ball case
    if isa(conddata{nc},'PcBallConductivity')
        
        % check radius
        fprintf(commandsfile,'##radius\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.radius);
        
    elseif isa(conddata{nc},'PcEllipsoidConductivity')
        
        % check radii
        fprintf(commandsfile,'##radius1\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.radii(1));
        fprintf(commandsfile,'##radius2\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.radii(2));
        fprintf(commandsfile,'##radius3\n');
        fprintf(commandsfile,'%f\n', conddata{nc}.radii(3));
        
        % check axes
        fprintf(commandsfile,'##axis1\n');
        fprintf(commandsfile,'cartesian\n');
        fprintf(commandsfile,'%f ', conddata{nc}.axes(1,:));
        fprintf(commandsfile,'\n');
        fprintf(commandsfile,'##axis2\n');
        fprintf(commandsfile,'cartesian\n');
        fprintf(commandsfile,'%f ', conddata{nc}.axes(2,:));
        fprintf(commandsfile,'\n');
        fprintf(commandsfile,'##axis3\n');
        fprintf(commandsfile,'cartesian\n');
        fprintf(commandsfile,'%f ', conddata{nc}.axes(3,:));
        fprintf(commandsfile,'\n');
    else
        error('Please, use piecewise constant conductivities');
    end
    
    % check amplitude
    fprintf(commandsfile,'##amplitude\n');

    if complex % at least one component is complex
        amplitudenc=conddata{nc}.amplitude;
        fprintf(commandsfile,'%f %f\n',real(amplitudenc),imag(amplitudenc));
    else
        fprintf(commandsfile,'%f\n',conddata{nc}.amplitude);
    end

    % check n
    fprintf(commandsfile,'##n\n');
    fprintf(commandsfile,'%d\n',conddata{nc}.n);
end


% conductivity saving file check and write
if savecond ~= 0
    fprintf(commandsfile,'#conductivityonmesh\n');
    % mesh
    fprintf(commandsfile,'%s\n',mesh);
    % conductivity
    fprintf(commandsfile,'conductivity_dnmapid_%d_true.dat\n',dnmapid);
end

% q saving file check and write
if saveq ~= 0
    fprintf(commandsfile,'#qonmesh\n');
    % mesh
    fprintf(commandsfile,'%s\n',mesh);
    % q
    fprintf(commandsfile,'q_dnmapid_%d_true.dat\n',dnmapid);
end

fclose(commandsfile);