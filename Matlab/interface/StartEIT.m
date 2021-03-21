function StartEIT(parallelize, complex, commands, log, cluster)
eitid = regexp(log,'\d*','once','match');
if parallelize
    if complex
        cmd = ['./bin/fast_eit_cmplx ',commands,' ',log];
    else
        cmd = ['./bin/fast_eit ',commands,' ',log];
    end
else
    if complex
        cmd = ['./bin/eit_cmplx ',commands,' ',log];
    else
        cmd = ['./bin/eit ',commands,' ',log];
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
        %disp(cmdout);
        if ~status
            error('Something went wrong. Please make sure libraries are installed and loaded correctly');
        end
        process = str2double(cmdout); % get process name on pc

        fileid = fopen('processes.txt','a+');
        fprintf(fileid,[num2str(process), ' eit', eitid,'\n']); % save processname
        fclose(fileid);
    end
    
    
    %[~, cmdout] = system([cmd,'& echo $!']); % run command in terminal 
    %process = str2double(cmdout); % get process name on pc

    %fileid = fopen('processes.txt','a+');
    %fprintf(fileid,[num2str(process), ' eit', eitid,'\n']); % save processname
    %fclose(fileid);
else
    cluster(jobid); % Run clusterscript
end
end