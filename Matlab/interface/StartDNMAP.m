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
    [~, cmdout] = system([cmd,'& echo $!']); % run command in terminal 
    process = str2double(cmdout); % get process name on pc

    fileid = fopen('processes.txt','a+');
    fprintf(fileid,[num2str(process), ' dnmap', dnmapid,'\n']); % save processname
    fclose(fileid);
else
    cluster(jobid); % Run clusterscript
end
end
