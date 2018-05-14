%% OVERALL GROWTH DYNAMICS OF MUTANT/WT POPULATIONS: deterministic simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL PARAMETERS:
% Numerical integration parameters:
t_ini = 0; % weeks
t_step = 0.1; % weeks
t_end = 100; % weeks

% CONSTANT PARAMETERS OF THE WT BACK SKIN (MLE VALUES)
lambda = 1.13; % week-1
r = 0.0506;
gamma = 2.8624; % week-1
mu = 0.8195; % week-1
dens = gamma / (lambda + gamma); % assumption of homeostasis %dens = 0.7170; 
m = lambda * dens / mu; % assumption of homeostasis %m = 0.9886;


%% SIMULATE TIME EVOLUTION OF OVERALL WT CELL POPULATION:
delta = 0;

% Initial populations in WT: (size of the "tissue unit block")
% we assume that the three populations are balanced in homeostasis:
x0=[1 0 0]; % populations of a,b,c = basal proliferating | basal differentiated | suprabasal cells
x0(2) = (x0(1)-dens*x0(1))./dens; % balanced with dynamics
x0(3) = m.*(x0(1)+x0(2));
x0_WT = x0./sum(x0); % relative populations

% Numerical integration
ode=@(t,x) Competition_det_eq(t,x,delta,r,lambda,gamma,mu);
[t,u_WT]=ode45(ode,[t_ini:t_step:t_end],x0_WT);
u_WT_all = sum(u_WT,2);

% Plotting population time course:
figure(1)
subplot(1,2,1)
plot(t,u_WT_all)
title('WT (populations)')
subplot(1,2,2)
hold on; plot(t,u_WT(:,1)); plot(t,u_WT(:,2)); plot(t,u_WT(:,3));
title('WT (populations)')


%% SIMULATE TIME EVOLUTION OF OVERALL MUTANT CELL POPULATION UNDER DIFFERENT PARAMETER CONDITIONS:
%delta = let it variable
frec_mutant = 0.01; % initial fraction of p53 mutant basal cells

% Numerical integration
% we loop for different parameter values:
param2loop = 1; % 1=mu | 2=gamma | 3=r
mu_all = mu.*10.^(-[0:5/20:5]);
gamma_all = gamma.*2.^([0:5/20:5]);
r_all = [0:0.5/20:0.5];
delta_all = [0:1/20:1];
u_Mut = [];
u_Mut_all = [];
for buc_mu = 1:size(mu_all,2)
    
    switch param2loop %provided the 3 parameter vectors have a size equal as mu_all
        case 1
            mu = mu_all(1,buc_mu);
        case 2
            gamma = gamma_all(1,buc_mu);
        case 3
            r = r_all(1,buc_mu);
    end
    
    for buc_delta = 1:size(delta_all,2)

        delta = delta_all(1,buc_delta); %delta = 0.05/r;

        % Initial population of Mutant cells:
        dens_mut = 1; % we assume that only "a" mutant cells are present at the beginning
        x0_Mut(1) = dens_mut.*(x0_WT(1) + x0_WT(2)).*(frec_mutant/(1-frec_mutant));
        x0_Mut(2) = (x0_Mut(1)-dens_mut*x0_Mut(1))./dens_mut;
        m_mut = 0;
        x0_Mut(3) = m_mut.*(x0_Mut(1)+x0_Mut(2));

        ode=@(t,x) Competition_det_eq(t,x,delta,r,lambda,gamma,mu);
        [t,u]=ode45(ode,[t_ini:t_step:t_end],x0_Mut);

        u_Mut(:,:,buc_delta,buc_mu) = u;
        u_Mut_all(:,buc_delta,buc_mu) = sum(u,2);
    end

end

% Plotting the population time course under different delta-parameter conditions:
% Select subdata of interest:
subdata = 1;
switch param2loop
    case 1
        mu = mu_all(1,subdata);
    case 2
        gamma = gamma_all(1,subdata);
    case 3
        subdata = 3;
        r = r_all(1,subdata);
end

figure(2)
plot(t,u_Mut_all(:,:,subdata));
set(gca,'YScale','log')
title('Mutant (overall population)')


%% RECONSTRUCTING GLOBAL (WT+MUTANT) POPULATION TIME COURSES AND CHARACETERISTIC MEASURMENTS IN THE TISSUE:
% Check overall population ratios over time:
u_Tot_sb = []; u_Tot_b = []; u_Mut1 = []; u_Mut2 = []; u_Mut3 = [];
u_Tot_sb(:,:,:) = u_WT(:,3) + u_Mut(:,3,:,:);
u_Tot_b(:,:,:) = u_WT(:,1) + u_WT(:,2) + u_Mut(:,1,:,:) + u_Mut(:,2,:,:);
u_Mut1(:,:,:) = u_Mut(:,1,:,:); u_Mut2(:,:,:) = u_Mut(:,2,:,:); u_Mut3(:,:,:) = u_Mut(:,3,:,:);
% overall SB/B ratio:
m_Tot = u_Tot_sb ./ u_Tot_b;
% %p53+ basal cells - estimated as (a_Mut + b_Mut) / (a_Tot + b_Tot)
% (density of basal cells is assummed constant over time and independent of whether cells are p53 or WT)
p53_b(:,:,:) = (u_Mut1 + u_Mut2) ./ u_Tot_b(:,:,:) .* 100;
% %p53+ projected area - estimated as c_Mut / c_Tot
% (density of sb cells changes over time but is assumed independent of whether cells are p53 or WT)
p53_sb(:,:,:) = u_Mut3 ./ u_Tot_sb(:,:,:) .* 100;
% %EdU+ basal cells
EdU_b(:,:,:) = (u_Mut1 + u_WT(:,1)) ./ u_Tot_b(:,:,:) .* 100;

% Calculating gradients in tissue parameters/measurement changes at 6 months vs 0 months:
timepoint1 = find(t>=0,1); % 0 weeks = 0 months
timepoint2 = find(t>=(365/2/7),1); % 26 weeks = 6 months
% Gradient in m:
m_step(:,:) = m_Tot(timepoint2,:,:) ./ m;
% Gradient in %EdU+ basal cells:
EdU_b_step(:,:) = EdU_b(timepoint2,:,:) ./ (dens.*100);
% Gradient in %p53+ basal cells:
p53_b_step(:,:) = p53_b(timepoint2,:,:) - p53_b(timepoint1,:,:);
% Gradient in %p53+ suprabasal cells:
p53_sb_step(:,:) = p53_sb(timepoint2,:,:) - p53_sb(timepoint1,:,:);
% Gradient in thickness:
% we assume that size of SB cells is phi times the size of B cells
phi = [3.5 3.5]; % 3.5 is the normalized lateral surface area of each SB cell in p53* at 0-3 months (assumed as in WT) and 6 months
thick_step(:,:) = (1.1.*u_Tot_b(timepoint2,:,:) + phi(2).*u_Tot_sb(timepoint2,:,:)) ./ (u_Tot_b(timepoint1,:,:) + phi(1).*u_Tot_sb(timepoint1,:,:));


%% PLOTTING INFERRED CHANGES IN TISSUE PROPERTIES (i.e. in the B and SB compartments and in tissue thickness):
figure(3)

% Proportion of p53 mutant basal cells
ax1 = subplot(3,4,((param2loop-1)*4)+1);
imagesc(p53_b_step',[0 20])
xlabel('\Delta'); set(gca,'XTick',[1:4:length(delta_all)]); set(gca,'XTickLabel',delta_all(1:4:length(delta_all)));
switch param2loop
    case 1; ylabel('\mu (relative to WT)'); set(gca,'YTick',[1:4:length(mu_all)]); set(gca,'YTickLabel',mu_all(1:4:length(mu_all))./mu);
    case 2; ylabel('\gamma (relative to WT)'); set(gca,'YTick',[1:4:length(gamma_all)]); set(gca,'YTickLabel',gamma_all(1:4:length(gamma_all))./gamma);
    case 3; ylabel('r'); set(gca,'YTick',[1:4:length(r_all)]); set(gca,'YTickLabel',r_all(1:4:length(r_all)));
end
colormap(ax1,'parula')
colorbar
title('%p53* B cells, 6 months')

% Proportion of p53 mutant suprabasal cells
ax2 = subplot(3,4,((param2loop-1)*4)+2);
imagesc(p53_sb_step',[0 70])
xlabel('\Delta'); set(gca,'XTick',[1:4:length(delta_all)]); set(gca,'XTickLabel',delta_all(1:4:length(delta_all)));
switch param2loop
    case 1; ylabel('\mu (relative to WT)'); set(gca,'YTick',[1:4:length(mu_all)]); set(gca,'YTickLabel',mu_all(1:4:length(mu_all))./mu);
    case 2; ylabel('\gamma (relative to WT)'); set(gca,'YTick',[1:4:length(gamma_all)]); set(gca,'YTickLabel',gamma_all(1:4:length(gamma_all))./gamma);
    case 3; ylabel('r'); set(gca,'YTick',[1:4:length(r_all)]); set(gca,'YTickLabel',r_all(1:4:length(r_all)));
end
colormap(ax2,'parula')
colorbar
title('%p53* SB cells, 6 months')

% Tissue thickness
ax3 = subplot(3,4,((param2loop-1)*4)+3);
imagesc(thick_step',[1 2.5])
xlabel('\Delta'); set(gca,'XTick',[1:4:length(delta_all)]); set(gca,'XTickLabel',delta_all(1:4:length(delta_all)));
switch param2loop
    case 1; ylabel('\mu (relative to WT)'); set(gca,'YTick',[1:4:length(mu_all)]); set(gca,'YTickLabel',mu_all(1:4:length(mu_all))./mu);
    case 2; ylabel('\gamma (relative to WT)'); set(gca,'YTick',[1:4:length(gamma_all)]); set(gca,'YTickLabel',gamma_all(1:4:length(gamma_all))./gamma);
    case 3; ylabel('r'); set(gca,'YTick',[1:4:length(r_all)]); set(gca,'YTickLabel',r_all(1:4:length(r_all)));
end
colormap(ax3,'parula')
colorbar
title('Thickness, 6 months (rel. to WT)')

% SB/B ratio
ax4 = subplot(3,4,((param2loop-1)*4)+4);
imagesc(m_step',[0 2])
xlabel('\Delta'); set(gca,'XTick',[1:4:length(delta_all)]); set(gca,'XTickLabel',delta_all(1:4:length(delta_all)));
switch param2loop
    case 1; ylabel('\mu (relative to WT)'); set(gca,'YTick',[1:4:length(mu_all)]); set(gca,'YTickLabel',mu_all(1:4:length(mu_all))./mu);
    case 2; ylabel('\gamma (relative to WT)'); set(gca,'YTick',[1:4:length(gamma_all)]); set(gca,'YTickLabel',gamma_all(1:4:length(gamma_all))./gamma);
    case 3; ylabel('r'); set(gca,'YTick',[1:4:length(r_all)]); set(gca,'YTickLabel',r_all(1:4:length(r_all)));
end
mycolmap = [[0:0.02:1]', [0:0.02:1]', ones(51,1); ones(50,1), fliplr([0:0.02:0.98])' fliplr([0:0.02:0.98])'];
colormap(ax4,mycolmap)
colorbar
title('SB/B ratio, 6 months (rel. to WT)')
