function result = ispc_wsl

%%% S. Pasquet - V22.05.04
% Test if wsl is installed
% result = ispc_wsl

if ispc
    [s,~] = unix('wsl -l -v');
    if s==0
        result = 1;
    else
        result = 0;
    end
else
    result = 0;
end