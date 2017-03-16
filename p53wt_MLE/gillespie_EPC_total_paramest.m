% disp('GILLESPIE MODELLING OF EPITHELIUM DYNAMICS');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Constant parameters:
% dens = 0.5; % density of proliferating (stem) cells in the basal layer
% m = 1.5; % ratio of nucleated suprabasal cells to basal cells -> arbitrary value in this case
% lambda = 3; % week-1
% r = 0.05;
% gamma = dens*lambda/(1-dens); % week-1 // assumption of homeostasis
% mu = dens*lambda/m; % week-1 // assumption of homeostasis

function [nx_basal,nx_total,ntime] = gillespie_EPC_total_paramest(rtime,dens,lambda,r,gamma,mu,m,indiv)

    tic
    % General initial parameters:
    if (nargin < 8)
        indiv=100000;
    end
    timelim=rtime(1,end); %round(365/7); % weeks
    timeErrTol = 0.05; % weeks % should be smaller than the difference between consecutive time points!
    nval=size(rtime,2)-1; % no. of registered time points (excluding first one)
%     parcialtime=timelim./nval;
    ntime = ones(indiv,1) * rtime;%ones(indiv,1) * [0:nval].*parcialtime;
    nx=zeros(indiv,nval+1,3);

    nx_basal=zeros(indiv,nval+1);
    nx_total=zeros(indiv,nval+1);
    nxss=zeros(indiv,3);

    % ITERATION FOR DIFFERENT INDIVIDUAL CLONES
    for it=1:indiv
        
        % Initial cell content:
        x=[1  0  0]; % only single-cell clones with just one stem cell.

        % Initial variables:
        time=0;
        p=zeros(5,1);
        multiplo=0;

%         % Save the populations of cells at time 0:
%         nx(it,multiplo+1,:) = x(:);
%         ntime(it,multiplo+1)=time;
%         multiplo=multiplo+1;

        % ITERATION FOR EACH SINGLE CLONE:
        while (time <= timelim)

            x_pre = x;
            
            % Calculation of single-event probabilities:
            p(1) = lambda*r*x(1); % a -> a + a
            p(2) = lambda*(1-2*r)*x(1); % a -> a + b
            p(3) = lambda*r*x(1); % a -> b + b
            p(4) = gamma*x(2); % b -> c
            p(5) = mu*x(3); % c -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Random number generator:
            r1=rand;
            r2=rand;
            tau=-(1./pt)*log(r1);
            r2n=r2*pt;

            % Event selection:
            event = find(cumsum(p)>r2n,1);
            if (event==1)
                x(1)=x(1)+1;
            elseif (event==2)
                x(2)=x(2)+1;
            elseif (event==3)
                x(1)=x(1)-1;
                x(2)=x(2)+2;
            elseif (event==4)
                x(2)=x(2)-1;
                x(3)=x(3)+1;
            elseif (event==5)
                x(3)=x(3)-1;
            end

            % Calculate time to that event:
            time=time+tau;

            % Save the populations of cells at certain time points:
            if (multiplo < (nval+1)) && (time >= rtime(1,multiplo+1))
                if (time < (rtime(1,multiplo+1)+timeErrTol))
                    nx(it,multiplo+1,:) = x(:);
                    ntime(it,multiplo+1)=time;
                    multiplo=multiplo+1;
                elseif (time == Inf)
                    while (multiplo < (nval+1))
                        nx(it,multiplo+1,:) = x(:);
                        ntime(it,multiplo+1)=rtime(1,multiplo+1);
                        multiplo=multiplo+1;
                    end
                else
                    stop_iter = multiplo;
                    while ((multiplo < (nval)) && (stop_iter < (nval)))
                        if (time >= rtime(1,multiplo+2))
                            nx(it,multiplo+1,:) = x_pre(:);
                            ntime(it,multiplo+1)= rtime(1,multiplo+1);
                            multiplo=multiplo+1;
                        else
                            stop_iter = nval;
                        end
                    end
                    if (time < (rtime(1,multiplo+1)+timeErrTol))
                        nx(it,multiplo+1,:) = x(:);
                        ntime(it,multiplo+1) = time; 
                    else
                        nx(it,multiplo+1,:) = x_pre(:);
                        ntime(it,multiplo+1) = rtime(1,multiplo+1); 
                    end
                    multiplo=multiplo+1;
                end
            end

        end
        %multiplo

        % Final values of the variables:
        nxss(it,:)=x(:);
        
%         [it nxss(it,:)]
    end

    % Sum both types of basal cells to get basal-layer clone sizes:
    nx_basal = nx(:,:,1)+nx(:,:,2);
    % Sum the 3 types of cells to get total clone sizes:
    nx_total = nx(:,:,1)+nx(:,:,2)+nx(:,:,3);
    toc
end
