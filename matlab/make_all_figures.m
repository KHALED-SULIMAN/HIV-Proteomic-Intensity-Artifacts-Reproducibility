%% =====================================================================
%  make_all_figures.m
%  Generates every main-figure panel as a separate editable .fig plus
%  1200-DPI TIFF/PNG.  Edit DATA and OUTDIR below, then run.
%  Fonts: titles/axis labels Arial Bold 30, ticks/legends Arial Bold 24
%         (change inside fig_style.m to restyle everything at once)
% =====================================================================
clear; close all; clc;

THISDIR = fileparts(mfilename('fullpath'));
PROJ = fileparts(THISDIR);
DATA = fullfile(PROJ,'results','figure_data');
OUTDIR = fullfile(PROJ,'results','figures_matlab');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end
rd = @(f) readtable(fullfile(DATA,f),'VariableNamingRule','preserve');
NW = @() figure('Color','w','Units','pixels','Position',[100 100 1100 850]);
C1=[0.20 0.35 0.65]; C2=[0.80 0.25 0.20]; C3=[0.20 0.60 0.35]; C4=[0.55 0.35 0.70];

%% ---------------- FIGURE 1b : cohort flow ----------------
T = rd('F1b_cohort_flow.csv');
f=NW(); ax=axes(f);
b=bar(ax,T.n,'FaceColor',C1,'EdgeColor','k','LineWidth',1.5);
xticks(ax,1:height(T)); xticklabels(ax,T.stage); xtickangle(ax,30);
text(ax,1:height(T),T.n+2,string(T.n),'HorizontalAlignment','center', ...
     'FontName','Arial','FontWeight','bold','FontSize',24);
ylim(ax,[0 110]); fig_style(ax,'Cohort assembly','','Number of samples');
save_panel(f,'F1b_cohort_flow',OUTDIR);

%% ---------------- FIGURE 2a : method comparison ----------------
T = rd('F2a_method_auroc.csv');
f=NW(); ax=axes(f); hold(ax,'on');
bar(ax,T.AUROC_METS,'FaceColor',C1,'EdgeColor','k','LineWidth',1.5);
errorbar(ax,1:height(T),T.AUROC_METS,T.AUROC_METS-T.ci_low,T.ci_high-T.AUROC_METS, ...
         'k','LineStyle','none','LineWidth',2.5,'CapSize',18);
yline(ax,0.5,'--','LineWidth',2,'Color',[.4 .4 .4]);
xticks(ax,1:height(T)); xticklabels(ax,T.method); ylim(ax,[0.4 1]);
fig_style(ax,'Held-out clinical stratification','','AUROC (metabolic syndrome)');
save_panel(f,'F2a_method_auroc',OUTDIR);

%% ---------------- FIGURE 2b : lambda sweep ----------------
T = rd('F2b_lambda_sweep.csv');
f=NW(); ax=axes(f); hold(ax,'on');
plot(ax,T.lambda,T.AUROC_cluster,'-o','Color',C1,'LineWidth',3,'MarkerSize',12,'MarkerFaceColor',C1);
plot(ax,T.lambda,T.AUROC_METS   ,'-s','Color',C2,'LineWidth',3,'MarkerSize',12,'MarkerFaceColor',C2);
plot(ax,T.lambda,T.AUROC_highTgl,'-^','Color',C3,'LineWidth',3,'MarkerSize',12,'MarkerFaceColor',C3);
legend(ax,{'cluster (anchor)','METS (held out)','high triglycerides (held out)'},'Location','best');
ylim(ax,[0.5 1.05]);
fig_style(ax,'Anchor-weight sweep','Anchor weight \lambda','AUROC');
save_panel(f,'F2b_lambda_sweep',OUTDIR);

%% ---------------- FIGURE 2c : silhouette vs k ----------------
T = rd('F2c_silhouette_vs_k.csv');
f=NW(); ax=axes(f);
plot(ax,T.k,T.silhouette,'-o','Color',C4,'LineWidth',3,'MarkerSize',14,'MarkerFaceColor',C4);
fig_style(ax,'Endotype number selection','Number of clusters k','Silhouette score');
save_panel(f,'F2c_silhouette_vs_k',OUTDIR);

%% ---------------- FIGURE 3a : per-sample medians ----------------
T = rd('F3a_sample_medians.csv');
f=NW(); ax=axes(f); hold(ax,'on');
g=T.endotype; raw=T.median_raw; cor=T.median_corrected;
grp=[repmat({'Raw E0'},sum(g==0),1);repmat({'Raw E1'},sum(g==1),1); ...
     repmat({'Corr E0'},sum(g==0),1);repmat({'Corr E1'},sum(g==1),1)];
vals=[raw(g==0);raw(g==1);cor(g==0);cor(g==1)];
boxplot(ax,vals,grp,'Colors','k','Widths',0.6,'Symbol','o');
set(findobj(ax,'Type','Line'),'LineWidth',2);
fig_style(ax,'Global per-sample shift','','Median protein level');
save_panel(f,'F3a_sample_medians',OUTDIR);

%% ---------------- FIGURE 3b/3c : volcano before & after ----------------
for v = ["raw","corrected"]
    T = rd(sprintf('F3%s_volcano_%s.csv', ternary(v=="raw",'b','c'), v));
    f=NW(); ax=axes(f); hold(ax,'on');
    sig = T.q < 0.05;
    scatter(ax,T.diff(~sig),-log10(T.q(~sig)+1e-300),40,[.75 .75 .75],'filled');
    scatter(ax,T.diff(sig) ,-log10(T.q(sig) +1e-300),50,C2,'filled');
    yline(ax,-log10(0.05),'--k','LineWidth',2);
    xline(ax,0,'-','Color',[.5 .5 .5],'LineWidth',1.5);
    nUp=sum(sig & T.diff>0); nDn=sum(sig & T.diff<0);
    text(ax,0.98,0.95,sprintf('%d up / %d down',nUp,nDn),'Units','normalized', ...
        'HorizontalAlignment','right','FontName','Arial','FontWeight','bold','FontSize',24);
    fig_style(ax,sprintf('Differential proteins (%s)',v), ...
        'Mean difference (high - low risk)','-log_{10} FDR');
    save_panel(f,sprintf('F3%s_volcano_%s',ternary(v=="raw",'b','c'),v),OUTDIR);
end

%% ---------------- FIGURE 3d : ARI ----------------
T = rd('F3d_ARI.csv');
f=NW(); ax=axes(f);
b=bar(ax,[T.ARI_raw T.ARI_corrected],'EdgeColor','k','LineWidth',1.5);
b(1).FaceColor=C2; b(2).FaceColor=C3;
xticks(ax,1:height(T)); xticklabels(ax,strrep(T.comparison,'_',' '));
legend(ax,{'raw','corrected'},'Location','best');
fig_style(ax,'Which layer defines the endotype','','Adjusted Rand Index');
save_panel(f,'F3d_ARI',OUTDIR);

%% ---------------- FIGURE 3e : artifact flags ----------------
T = rd('F3e_flags.csv');
f=NW(); ax=axes(f);
b=bar(ax,T.flags,'FaceColor',C2,'EdgeColor','k','LineWidth',1.5); ylim(ax,[0 4.5]);
text(ax,1:height(T),T.flags+0.15,string(T.flags)+"/4",'HorizontalAlignment','center', ...
     'FontName','Arial','FontWeight','bold','FontSize',24);
xticks(ax,1:height(T)); xticklabels(ax,T.stage);
fig_style(ax,'Artifact diagnostic flags','','Flags raised (of 4)');
save_panel(f,'F3e_flags',OUTDIR);

%% ---------------- FIGURE 4a : embedding ----------------
T = rd('F4a_embedding.csv');
f=NW(); ax=axes(f); hold(ax,'on');
e0=T.endotype==0; e1=T.endotype==1;
scatter(ax,T.PC1(e0),T.PC2(e0),140,C2,'filled','MarkerEdgeColor','k','LineWidth',1.2);
scatter(ax,T.PC1(e1),T.PC2(e1),140,C1,'filled','MarkerEdgeColor','k','LineWidth',1.2);
legend(ax,{sprintf('E0 high-risk (n=%d)',sum(e0)),sprintf('E1 low-risk (n=%d)',sum(e1))}, ...
       'Location','best');
fig_style(ax,'Integrated multi-omics endotypes','PC1','PC2');
save_panel(f,'F4a_embedding',OUTDIR);

%% ---------------- FIGURE 4b : clinical profile ----------------
T = rd('F4b_clinical_profile.csv');
f=NW(); ax=axes(f); hold(ax,'on');
M=[T.high_risk_mean T.low_risk_mean];
Mn=M./max(abs(M),[],2);                       % scale each variable to compare
b=bar(ax,Mn,'EdgeColor','k','LineWidth',1.5); b(1).FaceColor=C2; b(2).FaceColor=C1;
for i=1:height(T)
    if ~isnan(T.p(i)) && T.p(i)<0.05
        text(ax,i,1.06,'*','HorizontalAlignment','center', ...
             'FontName','Arial','FontWeight','bold','FontSize',34);
    end
end
xticks(ax,1:height(T)); xticklabels(ax,T.variable); xtickangle(ax,30); ylim(ax,[0 1.2]);
legend(ax,{'high risk (E0)','low risk (E1)'},'Location','best');
fig_style(ax,'Clinical profile by endotype','','Scaled mean (* p<0.05)');
save_panel(f,'F4b_clinical_profile',OUTDIR);

%% ---------------- FIGURE 4c : ROC ----------------
T = rd('F4c_roc.csv');
f=NW(); ax=axes(f); hold(ax,'on');
plot(ax,T.fpr,T.tpr,'-','Color',C1,'LineWidth',4);
plot(ax,[0 1],[0 1],'--','Color',[.5 .5 .5],'LineWidth',2);
axis(ax,[0 1 0 1]);
fig_style(ax,'Metabolic syndrome prediction','False positive rate','True positive rate');
save_panel(f,'F4c_roc',OUTDIR);

%% ---------------- FIGURE 4d : heatmap ----------------
T = rd('F4d_heatmap_top40.csv');
prot = T{:,3:end}; names = T.Properties.VariableNames(3:end);
[~,ord]=sort(T.endotype); Z=zscore(prot(ord,:));
f=figure('Color','w','Position',[100 100 1500 950]); ax=axes(f);
imagesc(ax,Z'); colormap(ax,redbluecmap_local()); caxis(ax,[-2.5 2.5]);
cb=colorbar(ax); cb.Label.String='z-score';
set(cb,'FontName','Arial','FontWeight','bold','FontSize',20);
yticks(ax,1:numel(names)); yticklabels(ax,names);
set(ax,'FontName','Arial','FontWeight','bold','FontSize',14);
xlabel(ax,'Patients (sorted by endotype)','FontName','Arial','FontWeight','bold','FontSize',30);
title(ax,'Top 40 differential proteins','FontName','Arial','FontWeight','bold','FontSize',30);
save_panel(f,'F4d_heatmap_top40',OUTDIR);

%% ---------------- FIGURE 5a/5b : pathways ----------------
pw = {'F5a_pathways_UP_in_high_risk.csv','UP in high-risk','F5a_pathways_UP'; ...
      'F5b_pathways_DOWN_in_high_risk.csv','DOWN in high-risk','F5b_pathways_DOWN'};
for i=1:2
    T = rd(pw{i,1}); T = T(1:min(12,height(T)),:);
    q = T.("Adjusted P-value"); [~,o]=sort(q,'descend'); T=T(o,:); q=q(o);
    f=figure('Color','w','Position',[100 100 1500 900]); ax=axes(f);
    barh(ax,-log10(q),'FaceColor',ternary(i==1,C2,C1),'EdgeColor','k','LineWidth',1.5);
    yticks(ax,1:height(T)); yticklabels(ax,cellfun(@(s) s(1:min(45,numel(s))), ...
        T.Term,'UniformOutput',false));
    set(ax,'FontName','Arial','FontWeight','bold','FontSize',18);
    xlabel(ax,'-log_{10} FDR','FontName','Arial','FontWeight','bold','FontSize',30);
    title(ax,sprintf('Enriched pathways (%s)',pw{i,2}), ...
        'FontName','Arial','FontWeight','bold','FontSize',30);
    save_panel(f,pw{i,3},OUTDIR);
end

%% ---------------- FIGURE 5c/5d : drug candidates ----------------
dg = {'F5c_drugs_setA.csv','Set A: lipid / platelet','F5c_drugs_setA'; ...
      'F5d_drugs_setB.csv','Set B: cytokine / immune','F5d_drugs_setB'};
for i=1:2
    T = rd(dg{i,1}); T = T(1:min(15,height(T)),:);
    [~,o]=sort(T.n_targets); T=T(o,:);
    f=figure('Color','w','Position',[100 100 1400 900]); ax=axes(f);
    barh(ax,T.n_targets,'FaceColor',C3,'EdgeColor','k','LineWidth',1.5);
    yticks(ax,1:height(T)); yticklabels(ax,T.drug);
    set(ax,'FontName','Arial','FontWeight','bold','FontSize',20);
    xlabel(ax,'Number of pathway targets','FontName','Arial','FontWeight','bold','FontSize',30);
    title(ax,dg{i,2},'FontName','Arial','FontWeight','bold','FontSize',30);
    save_panel(f,dg{i,3},OUTDIR);
end

%% ---------------- FIGURE 6 : missing modality ----------------
T = rd('F6_missing_modality.csv');
f=NW(); ax=axes(f); hold(ax,'on');
x=T.missing_frac*100;
plot(ax,x,T.MOSAIC       ,'-o','Color',C2,'LineWidth',3.5,'MarkerSize',13,'MarkerFaceColor',C2);
plot(ax,x,T.PCA_mean     ,'--s','Color',C1,'LineWidth',3.5,'MarkerSize',13,'MarkerFaceColor',C1);
plot(ax,x,T.PCA_knn      ,'--^','Color',C3,'LineWidth',3.5,'MarkerSize',13,'MarkerFaceColor',C3);
plot(ax,x,T.complete_case,':d','Color',C4,'LineWidth',3.5,'MarkerSize',13,'MarkerFaceColor',C4);
legend(ax,{'MOSAIC (native masking)','PCA + mean imputation', ...
           'PCA + KNN imputation','complete-case only'},'Location','southwest');
fig_style(ax,'Robustness to missing omics blocks', ...
    '% patients missing one whole layer','AUROC (metabolic syndrome)');
save_panel(f,'F6_missing_modality',OUTDIR);

%% ---------------- FIGURE 7 : batch correction (3 panels) ----------------
T = rd('F7_batch_correction.csv');
lbl = strrep(T.method,'_',' ');
specs = {'artifact_flags','Artifact flags (of 4)','Artifact removal',C2,'F7a_flags'; ...
         'METS_AUROC','AUROC (metabolic syndrome)','Clinical signal retained',C1,'F7b_auroc'; ...
         'stability_ARI','Bootstrap ARI','Endotype stability',C3,'F7c_stability'};
for i=1:3
    f=NW(); ax=axes(f);
    bar(ax,T.(specs{i,1}),'FaceColor',specs{i,4},'EdgeColor','k','LineWidth',1.5);
    xticks(ax,1:height(T)); xticklabels(ax,lbl); xtickangle(ax,25);
    if i==2; ylim(ax,[0.5 1]); end
    fig_style(ax,specs{i,3},'',specs{i,2});
    save_panel(f,specs{i,5},OUTDIR);
end

%% ---------------- SUPPLEMENTARY S1a : mapping rates ----------------
T = rd('S1a_mapping_rates.csv');
f=NW(); ax=axes(f);
histogram(ax,T.mapping_rate,15,'FaceColor',C1,'EdgeColor','k','LineWidth',1.5);
fig_style(ax,'RNA-seq mapping rates','Mapping rate (%)','Number of samples');
save_panel(f,'S1a_mapping_rates',OUTDIR);

fprintf('\nAll figures written to:\n  %s\n', OUTDIR);

%% ---------------- local helpers ----------------
function out = ternary(cond,a,b)
if cond; out=a; else; out=b; end
end
function m = redbluecmap_local()
n=256; h=floor(n/2);
r=[linspace(0.15,1,h) linspace(1,0.75,n-h)]';
g=[linspace(0.35,1,h) linspace(1,0.15,n-h)]';
b=[linspace(0.70,1,h) linspace(1,0.15,n-h)]';
m=[r g b];
end
