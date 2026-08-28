%% =====================================================================
%  make_figure1a.m  —  Study-design schematic (Figure 1a)
%  Box-and-arrow schematic drawn natively in MATLAB.
%  Output: F1a_study_design .fig / .tiff / .png at 1200 dpi.
%  Fonts: Arial; headings 20 bold; body 13; small notes 11.
% =====================================================================
clear; close all; clc;

THISDIR = fileparts(mfilename('fullpath'));
PROJ = fileparts(THISDIR);
OUTDIR = fullfile(PROJ,'results','figures_matlab');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

FS_H = 20; FS_B = 13; FS_S = 11; FN = 'Arial';
cRNA=[0.20 0.35 0.65]; cPROT=[0.80 0.25 0.20]; cMETA=[0.20 0.60 0.35];
cCLIN=[0.55 0.35 0.70]; cGREY=[0.35 0.35 0.35];

f = figure('Color','w','Units','pixels','Position',[60 60 1750 1000]);
axes('Position',[0 0 1 1]); axis([0 175 0 100]); axis off; hold on;

% ---- TITLE ----
text(87.5,96,'Study design: assembly, artifact control, and analysis of a public four-omics HIV cohort', ...
    'HorizontalAlignment','center','FontName',FN,'FontSize',FS_H,'FontWeight','bold');

% ---- BLOCK 1: DATA SOURCES ----
text(20,89,'1.  Public data sources','HorizontalAlignment','center','FontName',FN,'FontSize',FS_H,'FontWeight','bold');
bx=4; bw=32; bh=8; gy=2.2; y0=78;
box_(bx,y0,        bw,bh,sprintf('RNA-seq  (SRA: PRJNA983231)\n95 libraries'),cRNA ,FS_B,'bold',FN);
box_(bx,y0-(bh+gy),bw,bh,sprintf('Olink affinity proteomics\n2,923 proteins'),cPROT,FS_B,'bold',FN);
box_(bx,y0-2*(bh+gy),bw,bh,sprintf('Untargeted metabolomics\n877 metabolites'),cMETA,FS_B,'bold',FN);
box_(bx,y0-3*(bh+gy),bw,bh,sprintf('Clinical phenotype\n(COCOMO cohort)'),cCLIN,FS_B,'bold',FN);

% ---- BLOCK 2: LINKAGE CHAIN ----
text(64,89,'2.  Identifier linkage (reconstructed)','HorizontalAlignment','center','FontName',FN,'FontSize',FS_H,'FontWeight','bold');
lx=46; lw=36; lh=6.5; ly=80;
box_(lx,ly,        lw,lh,'SRA run accession',cGREY,FS_B,'normal',FN);
box_(lx,ly-1*(lh+2),lw,lh,'library code  P20109_nnn',cGREY,FS_B,'normal',FN);
box_(lx,ly-2*(lh+2),lw,lh,sprintf('clinical "User" field\n(remove _R suffix)'),cGREY,FS_B,'normal',FN);
box_(lx,ly-3*(lh+2),lw,lh,'COCOMO_ID',cGREY,FS_B,'bold',FN);
for i=0:2
    darrow(lx+lw/2, ly-i*(lh+2), lx+lw/2, ly-(i+1)*(lh+2)+lh, cGREY);
end
text(64,ly-4*(lh+2)+3.5,'Linkage undocumented in public deposits; reconstructed in this study', ...
    'HorizontalAlignment','center','FontName',FN,'FontSize',FS_S,'FontAngle','italic','Color',cGREY);
for yy=[82 71.8 61.6 51.4]
    darrow(bx+bw+0.5, yy, lx-0.5, 70, [0.6 0.6 0.6]);
end

% ---- BLOCK 3: COHORT CONVERGENCE ----
text(122,89,'3.  Integrated cohort','HorizontalAlignment','center','FontName',FN,'FontSize',FS_H,'FontWeight','bold');
box_(104,74,36,9,'96 runs  \rightarrow  95 quantified  \rightarrow  89 linked',cGREY,FS_B,'normal',FN);
rectangle('Position',[110,58,24,12],'Curvature',0.2,'FaceColor',[0.93 0.95 1.0],'EdgeColor',[0.15 0.25 0.5],'LineWidth',3);
text(122,64,sprintf('n = 89\nall four layers matched'),'HorizontalAlignment','center','VerticalAlignment','middle','FontName',FN,'FontSize',FS_B+1,'FontWeight','bold');
darrow(lx+lw+0.5,64,103.5,78.5,[0.6 0.6 0.6]);
darrow(122,74,122,70.3,[0.15 0.25 0.5]);

% ---- BLOCK 4: ANALYSIS WORKFLOW ----
text(87.5,44,'4.  Analysis workflow','HorizontalAlignment','center','FontName',FN,'FontSize',FS_H,'FontWeight','bold');
steps={'Integration','Artifact diagnostics','Batch correction','Endotype derivation','Pathway enrichment','Drug mapping'};
subs ={'(deep vs linear)','(4 metadata-free tests)','(4 strategies)','(k = 2, held-out METS)','(lipid / platelet)','(ART-compatible)'};
scol ={cRNA,cPROT,cGREY,cCLIN,cMETA,[0.85 0.55 0.10]};
nw=26; nh=11; nx0=4; ny=26; gap=(175-2*nx0-6*nw)/5;
for i=1:6
    x=nx0+(i-1)*(nw+gap);
    rectangle('Position',[x ny nw nh],'Curvature',0.12,'FaceColor',[0.96 0.96 0.98],'EdgeColor',scol{i},'LineWidth',2.4);
    text(x+nw/2,ny+nh*0.62,steps{i},'HorizontalAlignment','center','FontName',FN,'FontSize',FS_B,'FontWeight','bold');
    text(x+nw/2,ny+nh*0.26,subs{i},'HorizontalAlignment','center','FontName',FN,'FontSize',FS_S,'FontAngle','italic','Color',[0.3 0.3 0.3]);
    if i<6
        darrow(x+nw+0.4, ny+nh/2, x+nw+gap-0.4, ny+nh/2, [0.4 0.4 0.4]);
    end
end
darrow(122,58,87.5,ny+nh+1.0,[0.15 0.25 0.5]);

% ---- SAVE ----
name='F1a_study_design';
% fix physical size so 1200 dpi stays within image limits (~11 x 6.3 in)
set(f,'Units','inches','PaperUnits','inches');
set(f,'PaperPositionMode','manual','PaperPosition',[0 0 11 6.3],'PaperSize',[11 6.3]);

savefig(f,fullfile(OUTDIR,[name '.fig']));

ok_t=false; ok_p=false;
try, exportgraphics(f,fullfile(OUTDIR,[name '.tiff']),'Resolution',1200); ok_t=true; catch ME1
    fprintf('tiff exportgraphics failed: %s\n',ME1.message); end
try, exportgraphics(f,fullfile(OUTDIR,[name '.png']),'Resolution',1200);  ok_p=true; catch ME2
    fprintf('png exportgraphics failed: %s\n',ME2.message); end
if ~ok_t
    try, print(f,fullfile(OUTDIR,[name '.tiff']),'-dtiff','-r600'); ok_t=true;
        fprintf('tiff saved at 600 dpi via print\n'); catch, end
end
if ~ok_p
    try, print(f,fullfile(OUTDIR,[name '.png']),'-dpng','-r600'); ok_p=true;
        fprintf('png saved at 600 dpi via print\n'); catch, end
end
fprintf('Saved %s  fig:1 tiff:%d png:%d  to %s\n', name, ok_t, ok_p, OUTDIR);

%% ================= local helper functions =================
function box_(x,y,w,h,txt,edgecol,fs,fw,fn)
    rectangle('Position',[x y w h],'Curvature',0.12,'FaceColor',[0.96 0.96 0.98], ...
        'EdgeColor',edgecol,'LineWidth',2.2);
    text(x+w/2,y+h/2,txt,'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName',fn,'FontSize',fs,'FontWeight',fw);
end

function darrow(x1,y1,x2,y2,col)
    plot([x1 x2],[y1 y2],'-','Color',col,'LineWidth',2.6);
    L=2.0; ang=atan2(y2-y1,x2-x1);
    p1=[x2-L*cos(ang-pi/7), y2-L*sin(ang-pi/7)];
    p2=[x2-L*cos(ang+pi/7), y2-L*sin(ang+pi/7)];
    patch([x2 p1(1) p2(1)],[y2 p1(2) p2(2)],col,'EdgeColor',col);
end
