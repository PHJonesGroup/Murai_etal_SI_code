%% LOG-LIKELIHOOD CALCULATION OF THE MATCH PDF vs EMPIRICAL DISTRIBUTION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% rfreqs: cell array of clone size frequencies {:,timepoints}
% myPDF: matrix of clone size frequencies (:,timepoints)
% timepoints: row-vector containing time points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LogLike] = logLike_calc(rfreqs,myPDF,timepoints)

LogLike_t = [];
LogLike = [];
for baz = 1:size(timepoints,2)
    if (size(myPDF,1) < size(rfreqs{:,baz},1))
        scale_up = size(rfreqs{:,baz},1);
        myPDF = [myPDF; zeros(scale_up-size(myPDF,1),size(myPDF,2))];
    end
    LogLike_t(1,baz) = sum ( rfreqs{:,baz}(find(rfreqs{:,baz}~=0),1) .* log(myPDF(find(rfreqs{:,baz}~=0),baz)) );
end
LogLike = sum(LogLike_t,2);
