function [ha, pos] = Plot1D(type,comp_reconid,include_true,comp_what,v,center,h,span)
% Mandatory input:
%   type         - {'conductivity', 'qhat', 'q'}
%   comp_reconid - matrix with columns representing plots, in which the
%                  elements are ids to types we want to compare in 1D plot
%   include_true - vector of dnmapids of true conductivity, if element 
%                - is zero, we don't plot the truth
%   comp_what    - vector of "comparison"-strings of length ==
%                  size(comp_reconid,2). One "comparison"-string in 
%                  {'method','ift','ngrid','fixed','pkappa'} for each
%                  column in comp_reconid.
%   v            - restriction line through [0,0,0] is given by vector v =
%                  [v1,v2,v3]
%   h            - fineness of grid on which type is interpolated
%   span         - [a,b] determines the y-window
%
% Output:
%   ha      - array of handles of the axes objects starting from upper 
%             left corner, going row-wise as in subplot
%   pos     - positions of the axes objects


if size(v,1)~=1
    v=v';
end

% check length of include_true and comp_reconid
if size(comp_reconid,2)~=length(include_true)
    error('Please choose comp_reconid and include_true of same length')
end

% Comparable parameters
parset = {'reconmethod', 'nd', 'zeta', 'ift', 'ngrid', 'pkappa'};


% Start and end of plot interval 
%a = -1-h;
%b = 1+h;


% number of plots
M = size(comp_reconid,2);
for i = M:-1:1
    if include_true(i)
        [filename,radial] = getfilename_dnmapid(include_true(i),type);
        if radial % real radially symmetric conductivity
            if strcmp(type,'conductivity')
                [r,c] = read_conductivity(filename);
            elseif strcmp(type,'qhat')
                [r,c] = read_qhat(filename,'meshgrid');
            elseif strcmp(type,'q')
                [r,c] = read_q(filename);
            else
                error('Please choose a type');
            end
            % 1D mesh
            x = [-flip(r);r(2:end)];
            v = [flip(c);c(2:end)];          
            xq = a:h:b;
            vq = interp1(x,v,xq);

            % insert to struct:
            truecond(i).x = xq;
            truecond(i).radial = 1;
            truecond(i).c = vq;
            truecond(i).true = 1;

        else
            if strcmp(type,'conductivity')
                [x,y,z,c] = read_conductivity(filename);
            elseif strcmp(type,'qhat')
                [x,y,z,c] = read_qhat(filename,'meshgrid');
            elseif strcmp(type,'q')
                [x,y,z,c] = read_q(filename);
            else
                error('Please choose a type');
            end
            % insert to struct:
            truecond(i).x = x;
            truecond(i).y = y;
            truecond(i).z = z;
            truecond(i).c = c;
            truecond(i).radial = 0;
            truecond(i).true = 1;
        end
    end
    
    % Go through all elements in column
    K = size(comp_reconid,1);
    for j = 1:K
        if ~comp_reconid(j,i)
            continue;
        end
        [filename,info] = getfilename_reconid(comp_reconid(j,i),type);

        % read file depending on type
        if strcmp(type,'conductivity')
            [x,y,z,c] = read_conductivity(filename);
        elseif strcmp(type,'qhat')
            [x,y,z,c] = read_qhat(filename,'meshgrid');
        elseif strcmp(type,'q')
            [x,y,z,c] = read_q(filename);
        else
            error('Please choose a type');
        end
        
        % index k
        k = (i-1)*K+j;
        % insert to struct:
        cond(k).x = x;
        cond(k).y = y;
        cond(k).z = z;
        cond(k).c = c;
        cond(k).radial = 0;
        cond(k).true = 0;
        % Get parameter values
        cond(k).reconid = info{2};
        cond(k).reconmethod = info{3};
        cond(k).nd = info{5};
        cond(k).zeta = info{6};
        cond(k).ift = info{7};
        cond(k).ngrid = info{9};
        cond(k).pkappa = info{11};
        dnmapdat = split(info{13},'.');
        cond(k).dnmapid = dnmapdat{1};
    end
end

% Determine arrangement and size of plots
f=figure;
if M == 1
    Nh = 1; % horizontal
    Nw = 1; % vertical
    width = 12;
    height = 11;
    gap = [0, 0];
    margh = [0.15, 0.1];
    margw = [0.1, 0.1];
elseif rem(M,2) == 0
    Nh = M/2; % horizontal
    Nw = 2; % vertical
    width = 22;
    height = ceil(M/2)*10;
    gap = [0.1, 0.1];
    margh = [0.1, 0.1];
    margw = [0.05, 0.05];
else
    error('Please choose {1,2,4,6,...} plots')
end

% Tight subplot
[ha,pos] = tight_subplot(Nh, Nw,gap,margh,margw);


%plot
for ii = 1:M
    axes(ha(ii));
    
    hold on
    if include_true(ii)
        % plot true type 
        if truecond(ii).radial
            plot(truecond(ii).x, truecond(ii).c);
        else
            if strcmp(type,'conductivity')
                % cross line through 0
                [xs,ys,zs,P] = cross_line(-1,1,h,v,center);

                % interpolate
                Vq = scatteredInterpolant(truecond(ii).x,truecond(ii).y,truecond(ii).z,truecond(ii).c,'linear','none');
                Vq = Vq(xs,ys,zs);

                % if complex
                if ~isreal(Vq)
                Vq = abs(real(Vq))+abs(imag(Vq));
                end

                % plot
                plot(P,Vq,'black','LineWidth',2);

            elseif strcmp(type,'qhat')

                % extract ximax
                xix = truecond(ii).x;
                ximax = -xix(1,1,1);

                % cross sectional plane
                [xs,ys,zs,P] = cross_line(-ximax,ximax,h,v,center);

                % interpolate
                Vq = interp3(truecond(ii).x,truecond(ii).y,truecond(ii).z,truecond(ii).c,xs,ys,zs);


                % if complex
                if ~isreal(Vq)
                Vq = abs(real(Vq))+abs(imag(Vq));
                end

                % plot
                plot(P,Vq,'LineWidth',1.5);


            elseif strcmp(type,'q')
                % cross sectional plane
                [xs,ys,zs,P] = cross_line(-1,1,h,v,center);

                % interpolate
                Vq = scatteredInterpolant(truecond(ii).x,truecond(ii).y,truecond(ii).z,truecond(ii).c,'linear','none');
                Vq = Vq(xs,ys,zs);

                % if complex
                if ~isreal(Vq)
                Vq = abs(real(Vq))+abs(imag(Vq));
                end

                % plot
                plot(P,Vq,'LineWidth',1.5);

            end      
        end
    end
    
    % plot reconstructions {conductivity, qhat, q}
    for jj = 1:K 
        if ~comp_reconid(jj,ii)
            continue;
        end
        %index 
        kk = (ii-1)*K+jj;
        
        % Plot depending on type
        if strcmp(type,'conductivity')
            % cross line through 0
            [xs,ys,zs,P] = cross_line(-1,1,h,v,center);

            % interpolate
            Vq = scatteredInterpolant(cond(kk).x,cond(kk).y,cond(kk).z,cond(kk).c,'linear','none');
            Vq = Vq(xs,ys,zs);

            % if complex
            if ~isreal(Vq)
            Vq = abs(real(Vq))+abs(imag(Vq));
            end
            linelist = {'-','--','-.',':'};
            % plot
            plot(P,Vq,linelist{jj},'LineWidth',1.5);

        elseif strcmp(type,'qhat')

            % extract ximax
            xix = cond(kk).x;
            ximax = -xix(1,1,1);

            % cross sectional plane
            [xs,ys,zs,P] = cross_line(-ximax,ximax,h,v,center);

            % interpolate
            Vq = interp3(cond(kk).x,cond(kk).y,cond(kk).z,cond(kk).c,xs,ys,zs);


            % if complex
            if ~isreal(Vq)
            Vq = abs(real(Vq))+abs(imag(Vq));
            end

            % plot
            plot(P,Vq,'LineWidth',1.5);
 


        elseif strcmp(type,'q')
            % cross sectional plane
            [xs,ys,zs,P] = cross_line(-1,1,h,v,center);

            % interpolate
            Vq = scatteredInterpolant(cond(kk).x,cond(kk).y,cond(kk).z,cond(kk).c,'linear','none');
            Vq = Vq(xs,ys,zs);

            % if complex
            if ~isreal(Vq)
            Vq = abs(real(Vq))+abs(imag(Vq));
            end
            linelist = {'-','--','-.',':'};
            % plot
            plot(P,Vq,linelist{jj},'LineWidth',1.5); 
                
            ax = gca;
        end
        
        % legend for each plot
        leg{jj} = cond(kk).(comp_what{ii});
    end
    
    % legend
    if include_true(ii)
        leg = ['true', leg];
    end
    legend(leg);
    clearvars leg
    parsetdiff = setdiff(parset,(comp_what{ii}));
    titletext = maketitle(cond((ii-1)*K+1),parsetdiff);

    if strcmp(type,'conductivity')
        titletext = ['$$\gamma$$', titletext]; 
        title(titletext, 'Interpreter', 'latex','FontWeight','normal');
    elseif strcmp(type,'qhat')
        titletext = ['$$\hat{q}$$', titletext]; 
        title(titletext, 'Interpreter', 'latex','FontWeight','normal')
    elseif strcmp(type,'q')
        titletext = ['$$q$$', titletext]; 
        title(titletext, 'Interpreter', 'latex','FontWeight','normal')
    end
    ax = gca;
    figsettings(ax,span);

end
    

% resize figure
set(gcf,'units','centimeters','position',[10,10,width,height])


end

function str = maketitle(cond,parsetdiff)
    str = '';
    for h = 1:length(parsetdiff)
        if strcmp(parsetdiff{h},'reconmethod')
            str = append(str,[', method: $$', cond.reconmethod, '$$']);
        elseif strcmp(parsetdiff{h},'nd')
            str = append(str,[', $$nd = ', cond.nd, '$$']);
        elseif strcmp(parsetdiff{h},'zeta')
            str = append(str, [', ', cond.zeta, ' $$|\zeta|$$']);
        elseif strcmp(parsetdiff{h},'ift')
            str = append(str, [', ', cond.ift]);
        elseif strcmp(parsetdiff{h},'ngrid')
            str = append(str,[', $$ngrid = ', cond.ngrid, '$$']);
        elseif strcmp(parsetdiff{h},'pkappa')
            str = append(str,[', $$\kappa = ', cond.pkappa','$$']);
        end
    end
end
function figsettings(ax,span)
     
    grid off
    ax.FontSize = 12; 
    box on;
    ax.LineWidth = 1.5;
    ylim(span)
end
function [filename,info] = getfilename_reconid(reconidid,type)
filenameshort = sprintf('results/%s/%s_%d_*',type,type,reconidid);
matches = struct2cell(dir(filenameshort));
s2 = size(matches,2);
% check if reconid is unique
if s2 > 1
    error('Nonunique reconstruction id, = %d. Please choose a unique id', reconidid);
elseif s2 == 0
    error('No such reconstruction id, reconid = %d, exists', reconidid);
end
% get info
info = regexp(matches{1},'_','split');
filename = sprintf('results/%s/%s',type,matches{1});
end
function [filename,radial] = getfilename_dnmapid(dnmapidid,type)

filenameshort = sprintf('results/%s/%s_dnmapid_%d_*',type,type,dnmapidid);
matches = struct2cell(dir(filenameshort));
filename = sprintf('results/%s/%s',type,matches{1});
s2 = size(matches,2);
% check if reconid is unique
if s2 > 1
    error('Nonunique dnmap id, = %d. Please choose a unique id', dnmapidid);
elseif s2 == 0
    error('No such dnmap id, = %d, exists', dnmapidid);
end
% check if radial
fid = fopen(filename,'r');

% Check error
if fid == -1
    error(['Error when opening file ',filename]);
end

% Read the first line and remove leading and trailing white space
firstline=strtrim(fgetl(fid));
radial = strcmp(firstline,'radial');
end
function [xs,ys,zs,P] = cross_line(a,b,h,v,center)
% Crossectional line determined by v
v = v/norm(v);
if size(v,2) > 1
    v=v';
end
P = a:h:b;
xs = v(1,1)*P+center(1); % x-coordinates of line
ys = v(2,1)*P+center(2); % y-coordinates of line
zs = v(3,1)*P+center(3); % z-coordinates of line
end
