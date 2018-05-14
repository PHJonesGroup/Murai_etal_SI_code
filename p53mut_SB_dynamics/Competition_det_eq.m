%% ODEs GOVERNING MUTANT/WT POPULATION DYNAMICS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f]=Competition_det_eq(t,x,delta,r,lambda,gamma,mu)
    %--------------------------------------------------------------------------
    % Parameter values:
    %already provided by user

    %--------------------------------------------------------------------------
    % ODE set describing 'unbalanced single-progenitor model' dynamics:
    f=zeros(3,1); 
    f(1) = 2*delta*r*lambda*x(1);
    f(2) = lambda*x(1) - 2*delta*r*lambda*x(1) - gamma*x(2);
    f(3) = gamma*x(2) - mu*x(3);
end