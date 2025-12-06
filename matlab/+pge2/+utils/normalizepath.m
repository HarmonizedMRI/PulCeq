function pth = normalizepath(pth)
% function pth = normalizepath(path)
%
% Reformat path for pc/windows/mac, and add trailing slash/backslash 

if ispc
    pth = strrep(pth, '/', '\');
else
    pth = strrep(pth, '\', '/');
end

if ~endsWith(pth, filesep)
    pth = [pth filesep];
end

