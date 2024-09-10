function saveas_png(fh, outdir, fname)

% SAVEASPNG(FH, OUTDIR, FNAME)
% 
% Function to save a figure as png
% FH - Handle to a figure
% OUTDIR - Output directory, will be created if does not exist
% FNAME - File name (with extension)

% Make sure that fh is a handle to a figure
assert(isa(fh, 'matlab.ui.Figure'), 'Handle to a figure is not valid')

% Create output directory if doesn't exist
if exist(outdir,'dir') ~= 7
	mkdir(outdir)
end

% Save the figure
saveas(fh, fullfile(outdir, fname), 'png')
