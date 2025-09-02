%% progress bar local function
function pe_iterationCounter(counter,total)
    if total > 1
        
        % display message
        if counter==1
            fprintf(['counting (of ' num2str(total) ' total): '])
        end

        % display the count
        if any(counter==round(logspace(log10(1),log10(total),40)))
            fprintf([num2str(counter) ' '])
        end


        % display end of count
        if counter==total; fprintf(' ...completed \n'); end
    end
end

