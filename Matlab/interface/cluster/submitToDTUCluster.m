function process = submitToDTUCluster(jobid)

jobname = ['job_', jobid];

memcore = 500;
maxmem  = 800;
email   = 's153569@student.dtu.dk';

% update jobscript string
str = '#!/bin/sh\n\n#BSUB -q hpc\n';
str = append(str,['#BSUB -J ', jobname, '\n']);
str = append(str,['#BSUB -n ', num2str(ncores), '\n#BSUB -R "span[hosts=1]"\n']);
str = append(str,['#BSUB -R "rusage[mem=', num2str(memcore),'MB]"\n']);
str = append(str,['#BSUB -M ', num2str(maxmem), 'MB\n#BSUB -W 48:00\n']);
str = append(str,['#BSUB -u ', email, '\n#BSUB -N\n#BSUB -o OutputCluster.out\n#BSUB -e ErrorCluster.err\n']);
str = append(str,'export OMP_NUM_THREADS=$LSB_DJOB_NUMPROC\n');
str = append(str,['module load suitesparse\n',cmd]);

% name
jobscript = ['submitToCluster_',jobid,'.sh'];
fileid = fopen(jobscript,'w');
fprintf(fileid,str);
fclose(fileid);

% Submit job
cmdsub = ['bsub < ', jobscript];
[~, cmdout] = system(cmdsub);
% Get process
id = regexp(cmdout, '\d*','match');
process = id{end};
% Update current processes
fileid = fopen('Matlab/CallerCode/current_processes_dnmap.dat','a+');
fprintf(fileid,[process, ' dnmap', dnmapid,'\n']);
fclose(fileid);
% Delete jobscript
delete(jobscript);