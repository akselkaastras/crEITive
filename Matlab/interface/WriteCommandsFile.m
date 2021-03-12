function [complex, commands, log] = WriteCommandsFile(reconmethod,...
                                    ift,ngrid,truncrad,fixed,pkappa,mesh,...
                                    reconid,dnmapid)

% Check if eit id already exists
filenameshort = sprintf('results/conductivity/conductivity_%d*',reconid);
check = dir(filenameshort);
if ~isempty(check)
    error('Please choose unique eitid')
end

if ~(truncrad <= pi*(ngrid-1)^2/(2*ngrid))
    error(['Please choose truncation radius <=', num2str(pi*(ngrid-1)^2/(2*ngrid))])
end

% Making folders
folders = {'eit_commands','eit_logs'};
for i = 1:length(folders)
    if ~exist(folders{i},'dir')
        mkdir(folders{i})
    end
end           

if fixed == 1
    fix = 'fixed';
elseif fixed == 0
    fix = 'proportional';
else
    error('Please choose fixed (1) or proportional (0)');
end

% rows and cols of dnmap
[~,complex,nd] = read_dnmap(sprintf('results/dnmap/dnmap_%d.dat',dnmapid));
                                
commands = sprintf('commands/eit_commands/commands_eit_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);
log      = sprintf('logs/eit_logs/log_eit_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);
commandsfile = fopen(commands,'w');
logfile = fopen(log,'w');
fclose(logfile);

% Check error in reading
if commandsfile==-1
    error(['Error when opening file ',commandsfile]);
end

% dnmap check and write
fprintf(commandsfile,'#dnmap\n');
fprintf(commandsfile,'dnmap_%d.dat\n',dnmapid);


fprintf(commandsfile,'#nd\n');
fprintf(commandsfile,'%d\n',nd);

% method write
fprintf(commandsfile,'#method\n');
fprintf(commandsfile,'%s\n',reconmethod);

% ift write
fprintf(commandsfile,'#ift\n');
fprintf(commandsfile,'%s\n',ift);

% ngrid check and write
fprintf(commandsfile,'#ngrid\n');
fprintf(commandsfile,'%d\n',ngrid);

% truncrad check and write

fprintf(commandsfile,'#truncrad\n');
fprintf(commandsfile,'%d\n',truncrad);

% zeta write
fprintf(commandsfile,'#zeta\n');
if fixed == 1
    fprintf(commandsfile,'fixed\n');
elseif fixed == 0
    fprintf(commandsfile,'proportional\n');
end

% pkappa check and write
fprintf(commandsfile,'#pkappa\n');
fprintf(commandsfile,'%d\n',pkappa);

% mesh check and write
fprintf(commandsfile,'#mesh\n');
fprintf(commandsfile,'%s\n',mesh);

% cgo traces check and write
fprintf(commandsfile,'#cgostraces\n');
fprintf(commandsfile,'cgo_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat\n',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);

% q check and write
fprintf(commandsfile,'#q\n');
fprintf(commandsfile,'q_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat\n',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);
% qhat check and write
fprintf(commandsfile,'#qhat\n');
fprintf(commandsfile,'qhat_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat\n',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);


% conductivity check and write
fprintf(commandsfile,'#conductivity\n');
fprintf(commandsfile,'conductivity_%d_%s_nd_%d_%s_%s_ngrid_%d_kap_%d_dnmap_%d.dat\n',...
                    reconid,reconmethod,nd,fix,ift,ngrid,pkappa,dnmapid);

nd = sprintf('%02d', nd);
% dnmap0 check and write
% if it does not exists for nd, then we make one!
% and if it does exist we load it
fprintf(commandsfile,'#dnmap0\n');
dnmap0name = sprintf('dnmap0_nd%s.dat',nd);
if isfile(['results/dnmap0/', dnmap0name])
    fprintf(commandsfile,'read\n');
    fprintf(commandsfile,[dnmap0name,'\n']);
else 
    fprintf(commandsfile,'write\n');
    fprintf(commandsfile,[dnmap0name,'\n']);
end

% s0 check and write
% if it does not exists for nd, then we make one!
% and if it does exist we load it
fprintf(commandsfile,'#s0\n');
s0name = sprintf('s0_nd%s.dat',nd);
if isfile(['results/s0/', s0name])
    fprintf(commandsfile,'read\n');
    fprintf(commandsfile,[s0name,'\n']);
else 
    fprintf(commandsfile,'write\n');
    fprintf(commandsfile,[s0name,'\n']);
end

fclose(commandsfile);