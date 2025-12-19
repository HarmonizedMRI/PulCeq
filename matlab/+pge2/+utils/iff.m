function out = iff(cond, valTrue, valFalse)
% iff - Conditional operator (like ternary operator in C)
%   out = iff(cond, valTrue, valFalse) returns valTrue if cond is true,
%   otherwise returns valFalse.

if cond
    out = valTrue;
else
    out = valFalse;
end

