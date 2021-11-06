function [net, sing] = chartCalcs(cellObj, NumberofPerm)
    rewNet  = [];
    nRewNet = [];
    envR    = [];
    envNR   = [];
    
    for i = 1:4
        % Align matrices in proper direction
        if size(cellObj{i,1},1) == 1
            cellObj{i,1} = cellObj{i,1}';
        end
        if size(cellObj{i+1,2},1) == 1
            cellObj{i+1,2} = cellObj{i+1,2}';
        end
        
        % Calculations
        sing.n(i,1)   = length(cellObj{i,1});
        sing.n(i+1,2) = length(cellObj{i+1,2});
        
        nRewNet  = cat(1, nRewNet, cellObj{i,1});
        envNR    = cat(1, envNR, repmat(i, length(cellObj{i,1}), 1));
        rewNet   = cat(1, rewNet, cellObj{i+1,2});
        envR     = cat(1, envR, repmat(i+1, length(cellObj{i+1,2}), 1));
    end
    
    for i = 1:3
        sing.p_rank(i,1)   = ranksum(cellObj{i,1,1}, cellObj{i+1,1,1});
        sing.p_rank(i+1,2) = ranksum(cellObj{i+1,2,1}, cellObj{i+2,2,1});
    end
    
    [net.p(1), net.rho(1)] = permutation_pair_test_fast(nRewNet,envNR,NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
    [net.p(2), net.rho(2)] = permutation_pair_test_fast(rewNet,envR,NumberofPerm,'rankcorr'); %rankcorr for non parametric rank based; corr is linear
end








