%% progress bar local function
function pe_counter(itIdx,nIterations)
    if nIterations > 1
        if itIdx==1; fprintf(['iterations (of ' num2str(nIterations) '): ']); end
        if any(itIdx==round(logspace(log10(1),log10(nIterations),40)))
            fprintf([num2str(itIdx) ' '])
        end
        if itIdx==nIterations; fprintf(' ...iterations completed \n'); end
    end
end

