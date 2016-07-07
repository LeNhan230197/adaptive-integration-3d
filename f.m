% Target function
function ret = f(x, y, z)
	% Detailed implementation of the target function goes here
    % pause(0.1);
    
    if x^2 + y^2 + z^2 > 1
        ret = 1;
    else
        ret = 0;
    end 
end
