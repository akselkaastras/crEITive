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
    [~, cmdout] = system([cmd,'& echo $!']); % run command in terminal 
    process = str2double(cmdout); % get process name on pc

    fileid = fopen('processes.txt','a+');
    fprintf(fileid,[num2str(process), ' eit', eitid,'\n']); % save processname
    fclose(fileid);
else
    cluster(jobid); % Run clusterscript
end
end