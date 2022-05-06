function [status,result] = unix_cmd(cmd,wsl)

%%% S. Pasquet - V22.05.04
% Use appropriate command to run unix shells
% [status,result] = unix_cmd(cmd)

if nargin < 2
    wsl = ispc_wsl;
end

if wsl
    [status,result] = unix_wsl(cmd); % Use WSL if installed on Windows
else
    [status,result] = unix(cmd); % Regular command for Unix systems
end