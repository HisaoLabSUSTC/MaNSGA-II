function [selData,selind] = DisSel(data,idealpoint,selNum)
   
   % get the size of the input solutions
   [numdata, M] = size(data);

    % if the number of solutions is smaller than the selected num, return the
    % solutions
    if numdata <= selNum
        selData = data;
        selind = 1:numdata;
    else    
        % the first solution is the solution closest to the estimated ideal point
        maxdata = max(data);  
        norm_data = (data - idealpoint)./(maxdata - idealpoint); % normalization
        idealpoint = zeros(1, M);  % after normalization, the estimated ideal point becomes a zero vector
        dist_ideal = pdist2(norm_data, idealpoint);  % calculate the distance between each normalized point and the estimated ideal point 
        [~, firstIndex] = min(dist_ideal); % get the index of the solution that is closest to the estimated ideal point   
    
        % projection to the hyperplane
        anorm = norm_data./sum(norm_data, 2);

        % initialize the selData as an empty set 
        selData = []; 
        selind = []; 

        % initialize v and d 
        v = zeros(numdata, 1);     % zeroes represent not selected
        d = ones(numdata, 1)*inf;   % initialize distance to a very large number, i.e., inf

        % select the first solution
        selData = [selData; norm_data(firstIndex, :)]; 
        v(firstIndex) = 1;  % 1 represent selected 

        for k = 1:numdata
            if v(k) == 0   % if it has not been selected
                d(k) = min(sqrt(sum((anorm(k, :) - anorm(firstIndex, :)).^2)), d(k));  % calculate the distance 
            end
        end

        sz_selData = size(selData, 1); % number of selected solutions
        selind = [selind; firstIndex]; % update the index of the selected solutions
        
        
        while sz_selData < selNum
            % Choose a solution from the candidate subset that is the furthest distance away from the selected subset
            z = find(v == 0); 
            max_distance = max(d(z)); 
            j = find(d == max_distance);
            j = j(1);
            selData = [selData; data(j, :)];
            v(j) = 1; 
            d(j) = inf;

            for k=1:numdata
                if v(k) == 0
                 d(k) = min(sqrt(sum((anorm(k, :) - anorm(j, :)).^2)), d(k));
                end
            end

        sz_selData = size(selData, 1); 
        selind = [selind; j];
        end   
    end
  
end
