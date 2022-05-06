function path_unix = unix_wsl_path(path,wsl)

% S. Pasquet - V22.05.04
% Create correct path for WSL in Windows 11
% new_path = unix_wsl_path(path)

if nargin < 2
    wsl = ispc_wsl;
end

if wsl
    path_unix = replace(path,'\','/');
    cmd1 = sprintf('wslpath %s',path_unix);
    [~,path_unix] = unix_cmd(cmd1);
    path_unix = strrep(path_unix, newline, '');
else
    path_unix = path;
end