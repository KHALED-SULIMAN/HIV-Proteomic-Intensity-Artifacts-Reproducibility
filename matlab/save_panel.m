function save_panel(f, name, outdir)
% Robust save: editable .fig + 1200-DPI TIFF + PNG.
% Each format is attempted independently so one failure cannot abort the run.
if ~exist(outdir,'dir'); mkdir(outdir); end

ok_fig=false; ok_tif=false; ok_png=false;

% keep figure a sane physical size so 1200 dpi is not astronomically large
try
    set(f,'Units','inches');
    p=get(f,'Position');
    if p(3)>9 || p(4)>7
        scale=min(9/p(3),7/p(4));
        set(f,'Position',[p(1) p(2) p(3)*scale p(4)*scale]);
    end
    set(f,'PaperPositionMode','auto');
catch
end

% 1) editable .fig
try
    savefig(f, fullfile(outdir,[name '.fig'])); ok_fig=true;
catch ME
    fprintf('    [fig ] failed: %s\n', ME.message);
end

% 2) TIFF at 1200 dpi
try
    if isgraphics(f)
        exportgraphics(f, fullfile(outdir,[name '.tiff']), 'Resolution', 1200);
        ok_tif=true;
    end
catch
    try
        print(f, fullfile(outdir,[name '.tiff']), '-dtiff', '-r1200'); ok_tif=true;
    catch ME2
        fprintf('    [tiff] failed: %s\n', ME2.message);
    end
end

% 3) PNG at 1200 dpi (fall back to print, then to 600 dpi)
try
    if isgraphics(f)
        exportgraphics(f, fullfile(outdir,[name '.png']), 'Resolution', 1200);
        ok_png=true;
    end
catch
    try
        print(f, fullfile(outdir,[name '.png']), '-dpng', '-r1200'); ok_png=true;
    catch
        try
            print(f, fullfile(outdir,[name '.png']), '-dpng', '-r600'); ok_png=true;
            fprintf('    [png ] saved at 600 dpi (1200 exceeded limits)\n');
        catch ME3
            fprintf('    [png ] failed: %s\n', ME3.message);
        end
    end
end

fprintf('  %s  fig:%d tiff:%d png:%d\n', name, ok_fig, ok_tif, ok_png);
end
