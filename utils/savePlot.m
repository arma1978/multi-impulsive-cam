function savePlot(name, path, type, res)
% savePlot  Save the current figure as both PNG and .fig.
%
%   savePlot(name, path, type, res)
%
%   Inputs:
%     name - filename without extension (string)
%     path - output directory (string)
%     type - print driver string passed to print(), e.g. '-dpng'
%     res  - resolution string, e.g. '-r300'
%   Output:
%     none (writes figure files to disk)

ax = gca; 
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';

print(strcat(fullfile(path,name),'.png'),type,res);
savefig(strcat(fullfile(path,name),'.fig'));
end

