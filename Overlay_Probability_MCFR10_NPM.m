% Produces overlay probability maps of VLSM z maps with the Jülich
% Corticospinal tract
tic
clc
clear all
close all

cf = pwd;

fiber_ct = load_nii('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Motor_4p_r_MNI_FR10.nii');
%fiber_ct = load_nii('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Motor_4a_r_MNI_FR10.nii');

binary_roi = fiber_ct.img > 0;
size_fiber_ct = sum(sum(sum(binary_roi)));

ct = fiber_ct.img;
ct_el = [unique(ct(:)), histc(ct(:),unique(ct(:)))];
ct_el = ct_el(:,1);
ct_el(ct_el==0) = [];
prob_fiber = 10:10:length(ct_el)*10;

% Get Z Maps
cd E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Acute_Chron_rec_excl
data_folders = dir('Img*.nii');
data_folders_neu = cell(length(data_folders),1);
for w=1:length(data_folders)
    data_folders_neu{w} = data_folders(w).name;
end    
data_folders = data_folders_neu;
clear data_folders_neu w

results_overall = zeros(numel(data_folders),length(ct_el));
names_ges = cell(length(data_folders),1);
prop_roi_ges = zeros(length(data_folders),1);
prop_res_ges = zeros(length(data_folders),1);

thresh = [3.80 3.71 3.89 3.61 3.61 3.43]; % Statistical thersholds

% Calculate overlay
for i = 1:numel(data_folders)
        
    file = data_folders{i}; % Daten aktuell
    res = load_nii(file); 

    bin_res = res.img >= thresh(i);
    size_bin_res = sum(sum(sum(bin_res)));
    %size_bin_res_ges(j) = size_bin_res;
    
    overlap_roi = binary_roi.*bin_res;
    overlap_roi_size = sum(sum(sum(overlap_roi)));
    
    prop_roi = overlap_roi_size/size_fiber_ct*100; % Percent of whole ROI affected
    prop_res = overlap_roi_size/size_bin_res*100; % Percent of results on whole ROI
    
    prop_roi_ges(i) = prop_roi;
    prop_res_ges(i) = prop_res;
    
    %size_prob_ges = zeros(length(ct_el),1);
    pc_prob_ct_ges = zeros(length(ct_el),1);
    overlap_pc_res_prob_ges = zeros(length(ct_el),1);
    abs_size_prob = zeros(length(ct_el),1);
    
    for j = 1:length(ct_el)
    
        ct_prob = fiber_ct.img == ct_el(j);
        size_prob = sum(sum(sum(ct_prob)));
        abs_size_prob(j) = size_prob;
        %size_prob_ges(j) = size_prob;
        pc_prob_ct_ges(j) = size_prob/size_fiber_ct*100;

        overlap = bin_res.*ct_prob;
        size_overlap = sum(sum(sum(overlap)));

        overlap_pc_res_prob = size_overlap/size_prob*100; % Prozent des CSTs von Result betroffen
        overlap_pc_res_prob_ges(j) = overlap_pc_res_prob;

    end
    
    [path,name,ext] = fileparts(data_folders{i});
    pos_img = strfind(name, 'Img');
    pos_r = strfind(name, 'Arm')-2;
    %pos_rec = strfind(name, 'Rec');
    name = name(pos_img:pos_r);
    
    %results.(name) = overlap_pc_res_prob_ges;
    results_overall(i,:) = overlap_pc_res_prob_ges;
    names_ges{i} = name;
end

prob_roi_res = [prop_roi_ges prop_res_ges];

results_cst.pc_prob_ct = pc_prob_ct_ges;
results_cst.prob_roi_res = prob_roi_res;
results_cst.pc_affected = results_overall;

results_mcp = results_cst;
save('pc_prob_mcp_fr10.mat','results_mcp');


%fiber_ct.img = ct_prob;
%save_nii(fiber_ct,'Motor_4a_r_MNI_FR10_1.nii');


%% Plot
prob_fiber = num2cell(prob_fiber/10);
pc_ct = results_cst.pc_prob_ct;

bw = 1; % Barwidth
barmap = [0.6 0.6 0.6]; % [0.7 0.7 0.7] is grey
lw = 0.9; % Line Width
lws = 1; % Line Width subplot
fss = 10; % Font Size small
fsu = 8; % Font Size subplot
fsg = 9; % Font Size big
fsgg = 14; % Font Size bigger
f_eb = 1; % Size error bar
tl = 0.04; % Tick length
y_val = 80;

% Plot Poster
data_ges = results_overall;
label_analysis = {'Img: A - Behav: A' 'Img: A - Behav: A [rec excl]' 'Img: A - Behav: C' 'Img: A - Behav: C [rec excl]' 'Img: C - Behav: C' 'Img: C - Behav: C [rec excl]'};

figure(1)
for i = 1:length(label_analysis)
    subplot(3,2,i)
    data_plot = data_ges(i,:);
    bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
    hold on
    bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
    
    if i == 1 || i == 3 || i == 5
        ylabel('overlap [%]','Fontsize',fsg)
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
end

text(-15.8, 312, 'A)', 'Fontsize',fsgg);
text(-2.2, 312, 'B)', 'Fontsize',fsgg)

f = gcf;
set(f,'Units','centimeters','Position',[10 10 16 15]);

% Grafik speichern
name = 'Figure_3_Acute_Chronic_MCp';
path = pwd;
rez = 300; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dbmp',['-r',num2str(rez)],'-opengl') % save file 
close all

%% Plot PC affected
% y_val = 20;
% bw = 0.9; % Barwidth
% 
% label_analysis = {'Img: A - Behav: A' 'rec excl' 'Img: A - Behav: C' 'rec excl' 'Img: C - Behav: C' 'rec excl'};
% 
% figure(2)
% data_plot =  results_cst.prob_roi_res(:,1);
% bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
% 
% ylabel('overlap [%]','Fontsize',fsg)
% 
% step_acc = 5;
% yticks_acc = 0:step_acc:y_val;
% set(gca,'YTick',yticks_acc,'Fontsize',fsu,'ylim', [yticks_acc(1),yticks_acc(end)])
% 
% xticks = label_analysis;
% set(gca,'XTickLabel',xticks,'Fontsize',fsu);
% zu = 0.7;
% set(gca,'xlim', [1-zu,length(data_plot)+zu]);
% title('Overlay CST','Fontsize',fsu);
% 
% box off
% 
% f = gcf;
% set(f,'Units','centimeters','Position',[10 10 16 15]);
% 
% % Grafik speichern
% name = 'Figure_X_Overlay_Acute_Chronic_CST_rec_excl_bin';
% path = pwd;
% rez = 300; % resolution (dpi) of final graphic
% f = gcf; % f is the handle of the figure you want to export
% figpos = getpixelposition(f);
% resolution = get(0,'ScreenPixelsPerInch');
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
% close all

% %% T1 Analyses
% data_ges = [results_overall(8,:); results_overall(2,:); results_overall(6,:); results_overall(10,:)];
% label_analysis = {'Img: T1 C - Behav: C' 'Img: Clin A [adj] - Behav: A' 'Img: Clin A [adj] - Behav: C' 'Img: Clin C [adj] - Behav: C'};
% 
% figure(2)
% for i = 1:length(label_analysis)
%     subplot(length(label_analysis),1,i)
%     data_plot = data_ges(i,:);
%     bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
%     hold on
%     bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
%     
%     ylabel('overlap [%]','Fontsize',fsg)
%         
%     step_acc = 20;
%     yticks_acc = 0:step_acc:y_val;
%     set(gca,'YTick',yticks_acc,'Fontsize',fsu,'ylim', [yticks_acc(1),yticks_acc(end)])
%     
%     xticks = prob_fiber;
%     set(gca,'XTickLabel',xticks,'Fontsize',fsu);
%     zu = 0.7;
%     set(gca,'xlim', [1-zu,10+zu]);
%     title(label_analysis(i));
%     line([2.5 2.5],[2.5 y_val],'Linewidth',lw,'Color','k','LineStyle','--');
%     box off
% end
% 
% f = gcf;
% set(f,'Units','centimeters','Position',[10 8 8 15]);
% 
% % Grafik speichern
% name = 'prob_dist_t1';
% rez = 300; % resolution (dpi) of final graphic
% f = gcf; % f is the handle of the figure you want to export
% figpos = getpixelposition(f);
% resolution = get(0,'ScreenPixelsPerInch');
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
% close all

display('done');
toc
