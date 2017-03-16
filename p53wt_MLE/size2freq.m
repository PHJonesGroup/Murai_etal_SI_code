%% CLONE-SIZES-TO-FREQUENCIES CONVERTER:
% COLLECT FREQUENCIES OF CLONE SIZES, FROM THE EXPERIMENTAL DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% clonesizes: cell array {:,timepoints} or matrix [:,timepoints] containing clone sizes
% timepoints: row-vector containing time points
% vartype: type of variable that clonesizes represent: 1=cell array | 2=matrix
% clonesizes_ref: cell array or matrix of clone sizes {:,timepoints} used for cutoff1 (usually the TOTAL clone sizes)
% cutoff1: minimum clone size considered for the distribution
% clonesizes_ref2: cell array or matrix of clone sizes {:,timepoints} used for cutoff2 (usually the BASAL clone sizes)
% cutoff2: minimum clone size considered for the distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Freq, Freq_rel] = size2freq(clonesizes,timepoints,vartype,clonesizes_ref,cutoff1,clonesizes_ref2,cutoff2)

if (nargin < 4)
    clonesizes_ref = clonesizes;
    cutoff1 = 0;
end

switch vartype
    
    case 1
        Freq = cell(1,size(timepoints,2));
        Freq_rel = cell(1,size(timepoints,2));
        for aa = 1:size(timepoints,2)
            loc_prolif = find(clonesizes_ref{:,aa}>=cutoff1);
            if (nargin > 5)
                loc_prolif = loc_prolif(find(clonesizes_ref2{:,aa}(loc_prolif,1)>=cutoff2));
            end
            Freq{:,aa} = histc(clonesizes{:,aa}(loc_prolif,1),[0:1:max(clonesizes{:,aa})]);
            Freq_rel{:,aa} = Freq{:,aa}./sum(Freq{:,aa});
        end
        
    case 2      
        Freq = zeros(size([0:max(max(clonesizes))],2),size(timepoints,2));
        Freq_rel = zeros(size([0:max(max(clonesizes))],2),size(timepoints,2));
        for at = 1:size(timepoints,2) % We exclude clones with less than 2 cells, since in the experiments they may come from labelling of already differentiated cells.
            loc_prolif = find(clonesizes_ref(:,at)>=cutoff1);
            if (nargin > 5)
                loc_prolif = loc_prolif(find(clonesizes_ref2(loc_prolif,at)>=cutoff2));
            end
            Freq(:,at) = histc(clonesizes(loc_prolif,at),[0:1:max(max(clonesizes))]);
            Freq_rel(:,at) = Freq(:,at)./sum(Freq(:,at));
        end

end
