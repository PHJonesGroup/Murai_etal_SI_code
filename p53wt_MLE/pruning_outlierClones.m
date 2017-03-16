function [clonesizes1, clonesizes2] = pruning_outlierClones(clonesizes1,clonesizes2,time4pruning)

    for ptime = 1:size(time4pruning,2)
        buc_prone = time4pruning(ptime);
        size(clonesizes1{:,buc_prone})
        mean_logx = mean(clonesizes1{:,buc_prone});
        sd_logx = std(clonesizes1{:,buc_prone},0);
        logx_99th = 2.33*sd_logx+mean_logx;
        clonesizes1{:,buc_prone} = clonesizes1{:,buc_prone}(find(clonesizes1{:,buc_prone}<logx_99th),1);
        clonesizes2{:,buc_prone} = clonesizes2{:,buc_prone}(find(clonesizes1{:,buc_prone}<logx_99th),1);
        size(clonesizes1{:,buc_prone})
    end

end
