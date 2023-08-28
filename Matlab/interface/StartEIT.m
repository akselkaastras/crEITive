function StartEIT(parallelize, complex, commands, log, cluster, email)
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


if ~isa(cluster,'function_handle')
    disp('Please provide clusterscript')
    return;
else
    cluster(eitid,cmd,email); % Run clusterscript
end
end