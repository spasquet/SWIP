function [status,result] = unix_wsl(cmd)

%%% S. Pasquet - V22.05.04
% Unix command with WSL on Windows 11
% [status,result] = unix_wsl(cmd)

cmd_wsl = sprintf('wsl -e bash -li -c "%s"',cmd);
[status,result] = unix(cmd_wsl);
