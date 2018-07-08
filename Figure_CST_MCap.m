tic
clc
clear all
close all

cf = pwd;

load('pc_prob_cst.mat');
load('pc_prob_mca_fr10.mat');
load('pc_prob_mcp_fr10.mat');

%% Plot
prob_fiber = {1,2,3,4,5,6,7,8,9,10;};
pc_cst = results_cst.pc_prob_ct;
pc_mc = results_mca.pc_prob_ct;

cst_col = [0.6 0.6 0.6]; % color cst
mca_col = [0.3 0.3 0.3]; % color MC a
mcp_col = [0.1 0.1 0.1]; % color MC p

bw = 1; % Barwidth
lw = 0.9; % Line Width
lws = 1; % Line Width subplot
fss = 12; % Font Size small
fsu = 8; % Font Size subplot
fsg = 9; % Font Size big
fsgg = 14; % Font Size bigger
f_eb = 1; % Size error bar
tl = 0.04; % Tick length
y_val = 80;

% Plot CST, MC a/p
data_cst = results_cst.pc_affected;
data_mca = results_mca.pc_affected;
data_mcp = results_mcp.pc_affected;

label_analysis = {'Img: A - Behav: A' 'Img: A - Behav: A [rec excl]' 'Img: A - Behav: C' 'Img: A - Behav: C [rec excl]' 'Img: C - Behav: C' 'Img: C - Behav: C [rec excl]'};
pc_ct = [30.4380, 18.1377, 12.8338, 9.0801, 7.7756, 7.3835, 5.7268, 4.1848, 2.6407, 1.7990]; % = pc_rois from Figure_CST_MCap2.m

for i = 1:size(data_cst,1)
    
    cst_min_mca(i,:) = data_cst(i,:) - data_mca(i,:);
    cst_min_mcp(i,:) = data_cst(i,:) - data_mcp(i,:);
    mcp_min_mac(i,:) = data_mcp(i,:) - data_mca(i,:);
    
end
    
figure(1)
for i = 1:length(label_analysis)
    subplot(3,2,i)
    data_plot = data_cst(i,:);
    bar(data_plot,bw,'FaceColor',cst_col,'Edgecolor','k');    
    hold on
    data_plot = data_mcp(i,:);
    bar(data_plot,bw,'FaceColor',mcp_col,'Edgecolor','k');    
    hold on
    data_plot = data_mca(i,:);
    bar(data_plot,bw,'FaceColor',mca_col,'Edgecolor','k');    
    bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
    if i == length(label_analysis)
        hold on
        data_plot = data_cst(i,:);
        bar(data_plot,bw,'FaceColor',cst_col,'Edgecolor','k');
        hold on
        data_plot = data_mca(i,:);
        bar(data_plot,bw,'FaceColor',mca_col,'Edgecolor','k');    
    end
    
    if i == 1 || i == 3 || i == 5
        ylabel('overlap (%)','Fontsize',fsg)
    end
        
    step_acc = 20;
    yticks_acc = 0:step_acc:y_val;
    set(gca,'YTick',yticks_acc,'Fontsize',fsu,'ylim', [yticks_acc(1),yticks_acc(end)])
    
    xticks = prob_fiber;
    set(gca,'XTickLabel',xticks,'Fontsize',fsu);
    zu = 0.7;
    set(gca,'xlim', [1-zu,10+zu]);
    title(label_analysis(i));
    line([2.5 2.5],[2.5 y_val],'Linewidth',lw,'Color','k','LineStyle','--');
    box off
    
    set(gca,'ticklength',[0.02 0.02]);
end

text(-15.8, 312, 'A)', 'Fontsize',fsgg);
text(-2.2, 312, 'B)', 'Fontsize',fsgg);

f = gcf;
set(f,'Units','centimeters','Position',[10 5 15 14.5]);

% Legende
% legende_text_x = -5.5;
% legende_text_y = -25;
% fac = 3;
% text(legende_text_x, legende_text_y, 'CST', 'Fontsize',fsgg,'Color',cst_col-0.06);
% text(legende_text_x+fac, legende_text_y, 'MC a', 'Fontsize',fsgg,'Color',mca_col-0.06)
% text(legende_text_x+fac*2, legende_text_y, 'MC p', 'Fontsize',fsgg,'Color',mcp_col-0.06)

% Grafik speichern
name = 'Figure_3_Acute_Chronic_CST_MCap';
path = pwd;
rez = 600; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
close all

toc