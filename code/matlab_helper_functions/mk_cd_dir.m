function [status, msg] = mk_cd_dir(dirname, cd_dir)

% Check if a directory exists and if not, create it
if exist(dirname, 'dir') ~= 7
	[status,msg] = mkdir(dirname);
	if status ~= 1
		disp(msg)
	end
end

% Change directory if cd_dir is true
if cd_dir
    cd(dirname)
end

