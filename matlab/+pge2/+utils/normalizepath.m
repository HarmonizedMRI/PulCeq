function pth = normalizepath(pth)
% normalizepath - Reformat path for pc/windows/mac, and add trailing slash/backslash 
%
% function p = normalizepath(p)

if ispc
    pth = strrep(pth, '/', '\');
else
    pth = strrep(pth, '\', '/');
end

if ~endsWith(pth, filesep)
    pth = [pth filesep];
end

