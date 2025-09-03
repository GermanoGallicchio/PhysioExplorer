%% progress bar local function
function pe_counter(counter,total)
    persistent counting startTime

    if total > 1
        % set persistent variables
        

        % display message
        if counter==1
            fprintf(['counting (of ' num2str(total) ' total): '])
        end

        % start counting 
        if counter==1  &&  total>=1000
            counting = true;
            startTime = tic;
        end

        % display ETA
        if counter==100
            if counting
                elapsedTime = toc(startTime);
                remainingTime = (elapsedTime*total/counter) - elapsedTime;
                remainingMinutes = floor(remainingTime / 60);
                remainingSeconds = mod(remainingTime, 60);
                fprintf([ '     ' char(9203) ' ETA ' num2str(remainingMinutes) 'min, ' num2str(round(remainingSeconds)) 'sec ' char(9203) '     '])
                counting = false;
            end
        end

        % display the count
        if any(counter==round(logspace(log10(1),log10(total),40)))
            fprintf([num2str(counter) ' '])
        end


        % display end of count
        if counter==total; fprintf(' ...completed \n'); end
    end
end

