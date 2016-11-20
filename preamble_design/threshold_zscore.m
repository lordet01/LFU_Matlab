function [signals,avgFilter,stdFilter] = threshold_zscore(y,lag,threshold,influence)
% Initialise signal results
signals = zeros(length(y),1);
% Initialise filtered series
filteredY = y(1:lag+1);
% Initialise filters
avgFilter(lag+1,1) = mean(y(1:lag+1));
stdFilter(lag+1,1) = std(y(1:lag+1));
% Loop over all datapoints y(lag+2),...,y(t)
for i=lag+2:length(y)
    % If new value is a specified number of deviations away
    if abs(y(i)-avgFilter(i-1)) > threshold*stdFilter(i-1)
        if y(i) > avgFilter(i-1)
            % Positive signal
            signals(i) = 1;
        else
            % Negative signal
            signals(i) = -1;
        end
        % Adjust the filters
        filteredY(i) = influence*y(i)+(1-influence)*filteredY(i-1);
        avgFilter(i) = mean(filteredY(i-lag:i));
        stdFilter(i) = std(filteredY(i-lag:i));
    else
        % No signal
        signals(i) = 0;
        % Adjust the filters
        filteredY(i) = y(i);
        avgFilter(i) = mean(filteredY(i-lag:i));
        stdFilter(i) = std(filteredY(i-lag:i));
    end
end
% Done, now return results
end