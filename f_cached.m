function ret = f_cached(x, y, z, val, MODE)
	% Specify cache look up precision
	PRECISION = 6;
	
	persistent cache;
    
    if isempty(cache)
        cache = containers.Map;
    end
	
	% Generate lookup key in cache
	% Automatic rounding is built-in with num2str function
	cache_key = num2str([x y z], PRECISION);
	
	% Mode: 0 -> CALCULATE MODE
	% Used when first calculating all points
	if MODE == 0
		% Check if key exists in cache
        % Update in V0.4: isfield -> isKey
        %     cache.(cache_key) -> cache(cache_key)
        if cache.isKey(cache_key) == 0
			% Calculate and store if key doesn't exist
			fval = f(x, y, z);
			cache(cache_key) = fval;
			% fprintf('->Cache not used\n'); % DEBUG
        else
            fval = cache(cache_key);
			%fprintf('Something is wrong');
        end
        
		ret = fval;
        
	% Mode: 1 -> RETRIEVE MODE
	% Used within checker function
	elseif MODE == 1	
		% Check if key exists in cache
        if cache.isKey(cache_key) == 0
			% Variation in V0.2 to return NaN if the value is not stored
			% This is to ensure we don't do extra calculation unless instructed
			fval = NaN;
		else 
			% Retrieve if key exists
			fval = cache(cache_key);
        end
        ret = fval;
    
    % Mode: 2 -> DUMP MODE
    % Use to retrieve the cache in which all points and data are stored
    elseif MODE == 2
        ret = cache;
	
	% MODE: 3 -> SET MODE
	elseif MODE == 3
		fval = val;
		cache(cache_key) = fval;
        ret = 0;
	end
end

