function p = ensuretrailingslash(p)
% ensuretrailingslash - Adds slash/backslash (platform dependent) to string p
%
% function p = ensuretrailingslash(p)

    if ~endsWith(p, filesep)
        p = [p filesep];
    end

end

