function fig_style(ax, ttl, xl, yl)
% Applies journal formatting: Arial Bold, title/axis 30 pt, ticks/legend 24 pt.
% Edit FS_TITLE / FS_TICK below to change sizes globally.
FS_TITLE = 30;   % titles and axis labels
FS_TICK  = 24;   % tick labels, legends, in-figure text
set(ax,'FontName','Arial','FontWeight','bold','FontSize',FS_TICK, ...
        'LineWidth',1.5,'Box','off','TickDir','out');
if nargin>1 && ~isempty(ttl)
    title(ax,ttl,'FontName','Arial','FontWeight','bold','FontSize',FS_TITLE);
end
if nargin>2 && ~isempty(xl)
    xlabel(ax,xl,'FontName','Arial','FontWeight','bold','FontSize',FS_TITLE);
end
if nargin>3 && ~isempty(yl)
    ylabel(ax,yl,'FontName','Arial','FontWeight','bold','FontSize',FS_TITLE);
end
lg = findobj(ax.Parent,'Type','Legend');
if ~isempty(lg)
    set(lg,'FontName','Arial','FontWeight','bold','FontSize',FS_TICK);
end
end
