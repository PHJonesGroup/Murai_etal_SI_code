%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading experimental data on total clone sizes:
load 'TotalCloneSizes_raw_data.mat'
% rtime (weeks)
% rx_basal_time#_indiv#
% rx_total_time#_indiv#
nmice = [4 4 3 5 5 4 3]; %depending on time point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GROUPING INPUT DATA:
rx_basal = {rx_basal_t1_indiv1 rx_basal_t2_indiv1 rx_basal_t3_indiv1 rx_basal_t4_indiv1 rx_basal_t5_indiv1 rx_basal_t6_indiv1 rx_basal_t7_indiv1; rx_basal_t1_indiv2 rx_basal_t2_indiv2 rx_basal_t3_indiv2 rx_basal_t4_indiv2 rx_basal_t5_indiv2 rx_basal_t6_indiv2 rx_basal_t7_indiv2; rx_basal_t1_indiv3 rx_basal_t2_indiv3 rx_basal_t3_indiv3 rx_basal_t4_indiv3 rx_basal_t5_indiv3 rx_basal_t6_indiv3 rx_basal_t7_indiv3; rx_basal_t1_indiv4 rx_basal_t2_indiv4 [] rx_basal_t4_indiv4 rx_basal_t5_indiv4 rx_basal_t6_indiv4 []; [] [] [] rx_basal_t4_indiv5 rx_basal_t5_indiv5 [] []};
rx_total = {rx_total_t1_indiv1 rx_total_t2_indiv1 rx_total_t3_indiv1 rx_total_t4_indiv1 rx_total_t5_indiv1 rx_total_t6_indiv1 rx_total_t7_indiv1; rx_total_t1_indiv2 rx_total_t2_indiv2 rx_total_t3_indiv2 rx_total_t4_indiv2 rx_total_t5_indiv2 rx_total_t6_indiv2 rx_total_t7_indiv2; rx_total_t1_indiv3 rx_total_t2_indiv3 rx_total_t3_indiv3 rx_total_t4_indiv3 rx_total_t5_indiv3 rx_total_t6_indiv3 rx_total_t7_indiv3; rx_total_t1_indiv4 rx_total_t2_indiv4 [] rx_total_t4_indiv4 rx_total_t5_indiv4 rx_total_t6_indiv4 []; [] [] [] rx_total_t4_indiv5 rx_total_t5_indiv5 [] []};

rx_basal_all = {};
rx_total_all = {};
for aa = 1:size(rtime,2)
    rx_basal_all{1,aa} = [];
    rx_total_all{1,aa} = [];
    for ae = 1:max(nmice)
        rx_basal_all{1,aa} = [rx_basal_all{1,aa}; rx_basal{ae,aa}];
        rx_total_all{1,aa} = [rx_total_all{1,aa}; rx_total{ae,aa}];
    end
end

%% VALUE FOR THE LINEAR SCALING SLOPE (tau^{-1}) OBTAINED FROM THE AVERAGE BASAL-LAYER CLONE SIZES OVER TIME:
est_tau_inv_mean = 0.0797;
est_tau_inv_min95ci = 0.0671;
est_tau_inv_max95ci = 0.0923;

%% VALUE FOR THE DIVISION RATE (lambda), AS TAKEN FROM HISTONE-DILUTION EXPERIMENTS:
lambda_mean = 1.16;
lambda_min95ci = lambda_mean-0.06;
lambda_max95ci = lambda_mean+0.06;

%% VALUE FOR THE SB/B CELL RATIO (m), ESTIMATED FROM CLONES CONTAINING AT LEAST (4 BASAL & 4 SUPRABASAL) CELLS:
m_mean = 1.0720;
m_min95ci = 0.9053;
m_max95ci = 1.2387;

%% PRUNING OF EXTREMELY LARGE CLONES AT LATE TIME POINTS (OUTLIERS > 99th PERCENTILE FROM A LOG-NORMAL DISTRIBUTION)
rx_basal_all_raw = rx_basal_all;
rx_total_all_raw = rx_total_all;

time4pruning = [6,7];
[rx_total_all, rx_basal_all] = pruning_outlierClones(rx_total_all,rx_basal_all,time4pruning);

%% OBTAINING EXPERIMENTAL CLONE SIZE FREQUENCIES (BINNED IN POWERS OF 2):
[rfreq_tot_all, rfreq_tot_all_rel] = size2freq(rx_total_all,rtime,1,rx_total_all,2,rx_basal_all,1);
[rfreq_tot_dim_all, rfreq_tot_dim_all_rel, dim_label] = size2freqbinned(rfreq_tot_all,rx_total_all,rtime,1);

%% CALCULATION OF LIKELIHOOD VALUES WITH DIFFERENT PARAMETER SETS:
% GRID SEARCH ON UNKNOWN PARAMETER VALUES WITHIN THE CONSTRAINTS GIVEN BY THE HOMEOSTATIC CONDITION/CONFIDENCE BOUNDS PROVIDED ABOVE:
lambda_all = lambda_min95ci:((lambda_max95ci-lambda_min95ci)/4):lambda_max95ci;
est_tau_inv_all = est_tau_inv_min95ci:((est_tau_inv_max95ci-est_tau_inv_min95ci)/4):est_tau_inv_max95ci;
m_all = m_min95ci:((m_max95ci-m_min95ci)/4):m_max95ci;

loop_param = 0;
all_Ltotal = [];
all_r_all = [];
all_pre_param = [];
for lup1 = 1:size(lambda_all,2)    
    lambda = lambda_all(lup1);    
    for lup2 = 1:size(est_tau_inv_all,2)        
        est_tau_inv = est_tau_inv_all(lup2);        
        for lup3 = 1:size(m_all,2)
            m = m_all(lup3);
            [lambda est_tau_inv m]
            loop_param = loop_param+1

            %% GRID SEARCH ON PARAMETER r:
            r_ini = 0;
            r_end = est_tau_inv/lambda-0.001; % r = est_tau_inv*dens/lambda... the value of DENS cannot be higher than 1
            r_nsteps = 11;
            r = r_ini;

            Ltotal = [];
            r_all = [];
            for buc1 = 1:(r_nsteps+1)

                [r]
                % Constant parameters:
                lambda; % week-1
                dens = r*lambda/est_tau_inv; % week-1 // assumption of homeostasis
                gamma = dens*lambda/(1-dens); % week-1 // assumption of homeostasis
                m;
                mu = dens*lambda/m; % week-1 // assumption of homeostasis

                % GILLESPIE SIMULATION OF THE PDF over time:
                [nx_basal,nx_total,ntime] = gillespie_EPC_total_paramest(rtime,dens,lambda,r,gamma,mu,m);
                
                % COLLECTING PROBABILITY VALUES FOR EACH CLONE SIZE (BINNED IN POWERS OF 2), FROM THE SIMULATION OUTCOME:
                [nfreq_tot, nfreq_tot_rel] = size2freq(nx_total,ntime,2,nx_total,2,nx_basal,1);
                [nfreq_tot_dim, nfreq_tot_dim_rel, dim_label] = size2freqbinned(nfreq_tot,nx_total,ntime,2);

                % CALCULATION OF LIKELIHOOD FUNCTION (v1 & v2):
                [Ltotal(buc1,1)] = logLike_calc(rfreq_tot_dim_all_rel,nfreq_tot_dim_rel,rtime);

                r = r + (r_end/r_nsteps);
                
            end
            r_all = r_ini:r_end/r_nsteps:r_end;
            
            all_Ltotal(:,loop_param) = Ltotal;
            all_r_all(loop_param,:) = r_all;
            all_pre_param(loop_param,:) = [lambda est_tau_inv m];
        end        
    end
end

%% RETRIEVING THE MAXIMUM LIKELIHOOD ESTIMATE OF THE PARAMETER VALUES:
% PARAMETER r: MLE
[xloc_opti,yloc_opti] = find(all_Ltotal == max(max(all_Ltotal)));
mle_mean.r      = all_r_all (yloc_opti,xloc_opti);

% PARAMETER r: CI
all_xloc_min95ci = [];
all_xloc_max95ci = [];
all_r_min95ci = [];
all_r_max95ci = [];
for buz = 1:size(all_Ltotal,2)
    [xloc_95ci,yloc_95ci] = find(all_Ltotal(:,buz) >= (max(max(all_Ltotal)) - 1.92));
    all_xloc_min95ci(1,buz) = min(xloc_95ci);
    all_xloc_max95ci(1,buz) = max(xloc_95ci);
    all_r_min95ci(1,buz) = all_r_all(buz,all_xloc_min95ci(1,buz));
    all_r_max95ci(1,buz) = all_r_all(buz,all_xloc_max95ci(1,buz));
end
mle_min95ci.r   = min(all_r_min95ci);
mle_max95ci.r   = max(all_r_max95ci);

% PARAMETER lambda:
lambda_mean;
lambda_min95ci;
lambda_max95ci;

% PARAMETER est_tau_inv:
est_tau_inv_mean;
est_tau_inv_min95ci;
est_tau_inv_max95ci;

% PARAMETER dens:   (= r*lambda/est_tau_inv)
mle_mean.dens       = mle_mean.r * all_pre_param(yloc_opti,1) / all_pre_param(yloc_opti,2);
mle_min95ci.dens    = mle_min95ci.r * lambda_min95ci / est_tau_inv_max95ci;
mle_max95ci.dens    = mle_max95ci.r * lambda_max95ci / est_tau_inv_min95ci;
if mle_min95ci.dens < 0; mle_min95ci.dens = 0; end
if mle_max95ci.dens > 1; mle_max95ci.dens = 1; end

% PARAMETER gamma:  (= dens*lambda/(1-dens)
mle_mean.gamma      = mle_mean.dens * all_pre_param(yloc_opti,1) / (1-mle_mean.dens);
mle_min95ci.gamma   = mle_min95ci.dens * lambda_min95ci / (1-mle_min95ci.dens);
mle_max95ci.gamma   = mle_max95ci.dens * lambda_max95ci / (1-mle_max95ci.dens);

% PARAMETER m:
m_mean;
m_min95ci;
m_max95ci;

% PARAMETER mu:     (= dens*lambda/m)
mle_mean.mu         = mle_mean.dens * all_pre_param(yloc_opti,1) / all_pre_param(yloc_opti,3);
mle_min95ci.mu      = mle_min95ci.dens * lambda_min95ci / m_max95ci;
mle_max95ci.mu      = mle_max95ci.dens * lambda_max95ci / m_min95ci;

% SUMMARY:
mle_mean
mle_min95ci
mle_max95ci