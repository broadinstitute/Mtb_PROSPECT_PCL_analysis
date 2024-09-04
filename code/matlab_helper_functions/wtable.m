function wtable(tbl, out, varargin)

config = struct(...
    'name', {'--delimiter','--write_vn','--write_rn','--quiet'},...
    'default', {'\t',true,false,false},...
    'type',{'char','logical','logical','logical'},...
    'help', {'Field delimiter (default \t)',...
        'Write variable names',...
        'Write row names',...
        'Turn off verbosity'});
    
opt = struct('prog', mfilename, 'desc', 'Write a tsv file',...
    'undef_action', 'ignore');

[args, help_flag] = mortar.common.ArgParse.getArgs(config, opt, varargin{:});

if help_flag
   return
end


assert(istable(tbl),'%s> Structure "%s" in the input is not a table',...
    mfilename, inputname(1))

if ischar(out)
    p = fileparts(out);
    if ~isempty(p)
        if exist(p, 'dir')~=7
			sprintf('%s> Directory: %s does not exist', mfilename, p)
			mk_cd_dir(p,false);
		end
    end
    
    try
        writetable(tbl, out,'Delimiter',args.delimiter,'FileType','text',...
            'WriteRowNames',args.write_rn,...
            'WriteVariableNames',args.write_vn);
    catch ME
        disp(ME)
    end
else
    error('%s> Output filename should be a string',mfilename)
end