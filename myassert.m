function myassert(cond, msg, bWarn)

if ~cond
    if nargin == 3 && bWarn
        warning(msg);
    else
        error(msg);
    end
end
    
end

