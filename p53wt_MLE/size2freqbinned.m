%% CLONE-SIZES-TO-BINNED FREQUENCIES CONVERTER:
% COLLECT FREQUENCIES OF CLONE SIZES (BINNED IN POWERS OF TWO), FROM THE EXPERIMENTAL DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% freqs: cell array  {:,timepoints} or matrix [:,timepoints] of clone size frequencies
% clonesizes: cell array {:,timepoints} or matrix [:,timepoints] of clone sizes
% timepoints: row-vector containing time points
% vartype: type of variable that clonesizes and freqs do both represent: 1=cell array | 2=matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Freq_binned, Freq_binned_rel, Bin_labels] = size2freqbinned(freqs,clonesizes,timepoints,vartype)

% We categorize the cell No. into ranges increasing in size in powers of 2
switch vartype
    
    case 1
        Freq_binned = cell(1,size(timepoints,2));
        Freq_binned_rel = cell(1,size(timepoints,2));
        for aa = 1:size(timepoints,2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Binned groups:
            max_dim = ceil(log2(max(clonesizes{:,aa})));
            % Bin_labels = cell(1,max_dim+2);
            Freq_binned{:,aa}(1,1) = freqs{:,aa}(1,1); % freq of clones with 0 cells
            Bin_labels(1) = cellstr(num2str(0));
            for ot = 0:max_dim
                if ot == 0  % freq of single-cell clones, i.e. No_cells = 2^0 = 1
                    Freq_binned{:,aa}(ot+2,1) = freqs{:,aa}(2^ot+1,1); % these are stored in 2nd position in rfreq
                else  % e.g. freq of (3,4)-cell clones -> ot=2
                    Freq_binned{:,aa}(ot+2,1) = 0;
                    for ut = 1:(2^(ot)-2^(ot-1)) % 2 values contained
                        if (2^(ot-1)+1+ut) <= max(clonesizes{:,aa})+1
                            Freq_binned{:,aa}(ot+2,1) = Freq_binned{:,aa}(ot+2,1) + freqs{:,aa}(2^(ot-1)+1+ut,1); % these are stored in 4th and 5th position in nfreq
                        end
                    end
                end

                % Making the labels of the categories (unnecessary):
                if (2^(ot)-2^(ot-1)) > 1 % when more than 1 value, print the range
                    Bin_labels(ot+2) = cellstr([num2str(2^(ot-1)+1) '-' num2str(2^ot)]);
                else % when only 1 value contained, print the single value
                    Bin_labels(ot+2) = cellstr(num2str(2^ot));
                end
            end
            Freq_binned_rel{:,aa} = Freq_binned{:,aa}./sum(Freq_binned{:,aa});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

    case 2
        max_dim = ceil(log2(max(max(clonesizes))));
        Freq_binned = zeros(max_dim+2,size(freqs,2)); % (we will include 0-cells clones, therefore is max_dim+2)
        Freq_binned_rel = zeros(max_dim+2,size(freqs,2)); % (we will include 0-cells clones, therefore is max_dim+2)
        Bin_labels = cell(1,max_dim+2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Binned groups:
        Freq_binned(1,:) = freqs(1,:); % freq of clones with 0 cells
        Bin_labels(1) = cellstr(num2str(0));
        for ot = 0:max_dim
            if ot == 0  % freq of single-cell clones, i.e. No_cells = 2^0 = 1
                Freq_binned(ot+2,:) = freqs(2^ot+1,:); % these are stored in 2nd position in rfreq
            else  % e.g. freq of (3,4)-cell clones -> ot=2
                for ut = 1:(2^(ot)-2^(ot-1)) % 2 values contained
                    if (2^(ot-1)+1+ut) <= max(max(clonesizes))+1
                        Freq_binned(ot+2,:) = Freq_binned(ot+2,:) + freqs(2^(ot-1)+1+ut,:); % these are stored in 4th and 5th position in nfreq
                    end
                end
            end

            % Making the labels of the categories:
            if (2^(ot)-2^(ot-1)) > 1 % when more than 1 value, print the range
                Bin_labels(ot+2) = cellstr([num2str(2^(ot-1)+1) '-' num2str(2^ot)]);
            else % when only 1 value contained, print the single value
                Bin_labels(ot+2) = cellstr(num2str(2^ot));
            end
        end
        Freq_binned_rel = Freq_binned./sum(Freq_binned,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
end
