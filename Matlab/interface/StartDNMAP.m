function StartDNMAP(method, parallelize, complex, commands, log, cluster,email)



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

if ~isa(cluster,'function_handle')
    disp('Please provide clusterscript')
    return;
else
    cluster(dnmapid,cmd,email); % Run clusterscript
end
end
