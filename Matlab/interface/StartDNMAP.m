function process = StartDNMAP(method, parallelize, complex, commands, log, cluster)



dnmapid = regexp(log,'\d*','once','match');
if strcmpi(method,'radial')
    if parallelize
        cmd = ['./bin/fast_radsym ',commands,' ',log];
    else
        cmd = ['./bin/radsym ',commands,' ',log];
    end
else
    if parallelize
        if complex
            cmd = ['./bin/fast_pcc_cmplx ',commands,' ',log];
        else
            cmd = ['./bin/fast_pcc ',commands,' ',log];
        end
    else
        if complex
            cmd = ['./bin/pcc_cmplx ',commands,' ',log];
        else
            cmd = ['./bin/pcc ',commands,' ',log];
        end
    end 
end
if ~cluster
    if strcmpi(computer,'GLNXA64')
        s = input('Do you want me to load suitesparse via "module load"? (y/n): ','s');
        
        if isempty(s)
            s = 'n';
        end
        if strcmpi(s,'y')
            cmd = ['module load suitesparse;', cmd];
        end
        s = input(['Run the program from the terminal? (y/n): '], 's');
        if strcmpi(s,'y')
            disp('Program started ...');
        else
            error('Program aborted.');
        end
        [status, cmdout] = system([cmd,'& echo $!'],'-echo');
        %disp(cmdout);
        if status
            error('Something went wrong. Please make sure libraries are installed and loaded correctly');
        end
    else
        s = input(['Run the program from the terminal? (y/n): '], 's');
        if strcmpi(s,'y')
            disp('Program started ...');
        else
            error('Program aborted.');
        end
        [status, cmdout] = system([cmd,'& echo $!'],'-echo'); % run command in terminal 
        disp(cmdout);
        if ~status
            error('Something went wrong. Please make sure libraries are installed and loaded correctly');
        end
        process = str2double(cmdout); % get process name on pc

        fileid = fopen('processes.txt','a+');
        fprintf(fileid,[num2str(process), ' dnmap', dnmapid,'\n']); % save processname
        fclose(fileid);
    end
else
    cluster(jobid); % Run clusterscript
end
end
