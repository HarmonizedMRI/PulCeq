function removefiles(s)
% removefiles - Calls rm/del (unix/Win) to remove files
%
% function removefiles(s)

if isunix
    rmCmd = 'rm -f ';
end
if ispc
    rmCmd = 'del ';
end

system(sprintf('%s %s', rmCmd, s));
