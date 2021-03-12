function LogUpdate(type, ids)
% creates waitbar and prints log outputs of program ids entered

wb = waitbar(0,'Initializing...');
pause(0.5) % wait for commandsfile and logfile to be made



 % import process
K = length(ids);
% make cell of logfile names
for i = K:-1:1
    filenameshort = sprintf('logs/%s_logs/log_%s_%d*',type,type,ids(i));
    matches = struct2cell(dir(filenameshort));
    if isempty(matches)
        error('No log files of current process');
    else
        logfiles{i} = sprintf('logs/%s_logs/%s',type,matches{1});
        dirlog = dir(logfiles{i});
        logfiledate{i} = dirlog.date;
    end
end

% determine dnmaptype if type is dnmap
for i = K:-1:1
    dnmaptype{i} = type;
    if strcmp(type,'dnmap')
    tmp = regexp(logfiles{i},'_','split');
    dnmaptype{i} = tmp{5};
    end
end

% initialize
finished = zeros(K,1);
last = zeros(K,1);
new = ones(K,1);
% don't be too eager!

firsttime = 0;
pause(0.3)
waitbarstr = 'Running...';
while sum(finished)~=K
    for i = 1:K
        doc = fileread(logfiles{i});
        if new(i)
            % print log
            fprintf([doc,'\n']);
        end
        
        % detect if finished 
        if ~isempty(regexp(doc, 'Elapsed time','once')) && ~last(i)
            finished(i) = 1;
            doc = fileread(logfiles{i});
            % print log one last time
            fprintf([doc,'\n']);
            % output name id in case of successful job
            fprintf( ['Job done: ', num2str(ids(i)), '\n']);
            last(i) = 1;

            if ~firsttime 
                firsttime = firsttime+1;
                waitbarstr = 'Jobs done: ';
            end
            waitbarstr = append(waitbarstr, [num2str(ids(i)), ' ']);
        end
        pause(1)
        % detect change in logfile
        dirlog = dir(logfiles{i});
        new(i) = ~strcmp(logfiledate{i},dirlog.date);
        logfiledate{i} = dirlog.date;
        
        % detect user stopping
        
       
        %update waitbar
        waitbar(sum(finished)/K,wb,waitbarstr);
        
    end
    
end


end

