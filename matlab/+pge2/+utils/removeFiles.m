function removeFiles(s)

if isunix
    rmCmd = 'rm -f ';
end
if ispc
    rmCmd = 'del ';
end

system(sprintf('%s %s', rmCmd, s));
