function submitToDTUCluster(jobid,cmd,email)

jobname = ['job_', jobid];

ncores = 12;
memcore = 500;
maxmem  = 800;


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
% Delete jobscript
delete(jobscript);