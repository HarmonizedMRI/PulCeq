function p = ensuretrailingslash(p)
    if ~endsWith(p, filesep)
        p = [p filesep];
    end
end

