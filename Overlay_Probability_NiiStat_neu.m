% Produces overlay probability maps of VLSM z maps with the J�lich
% Corticospinal tract
tic
clc
clear all
close all

cf = pwd;

fiber_ct = load_nii('F:\ACLM\Matlab_Skripte\Fiber_ct_r.nii');
binary = fiber_ct.img > 0;
size_fiber_ct = sum(sum(sum(binary)));

ct = fiber_ct.img;
ct_el = [unique(ct(:)), histc(ct(:),unique(ct(:)))];
ct_el = ct_el(:,1);
ct_el(ct_el==0) = [];
prob_fiber = 10:10:length(ct_el)*10;

% Get Z Maps
data_folders = dir('threshZlesion*.nii');
data_folders_neu = cell(length(data_folders),1);
for w=1:length(data_folders)
    data_folders_neu{w} = data_folders(w).name;
end    
data_folders = data_folders_neu;
clear data_folders_neu w

results_overall = zeros(numel(data_folders),length(ct_el));
names_ges = cell(length(data_folders),1);

% Calculate overlay
for i = 1:numel(data_folders)
        
    file = data_folders{i}; % Daten aktuell
    res = load_nii(file); 

    bin_res = res.img > 0;
    size_bin_res = sum(sum(sum(bin_res)));
    %size_bin_res_ges(j) = size_bin_res;
    
    %size_prob_ges = zeros(length(ct_el),1);
    pc_prob_ct_ges = zeros(length(ct_el),1);
    overlap_pc_res_prob_ges = zeros(length(ct_el),1);
    
    for j = 1:length(ct_el)
    
        ct_prob = fiber_ct.img == ct_el(j);
        size_prob = sum(sum(sum(ct_prob)));
        %size_prob_ges(j) = size_prob;
        pc_prob_ct_ges(j) = size_prob/size_fiber_ct*100;

        overlap = bin_res.*ct_prob;
        size_overlap = sum(sum(sum(overlap)));

        overlap_pc_res_prob = size_overlap/size_prob*100; % Prozent des motorischen Systems von Result betroffen
        overlap_pc_res_prob_ges(j) = overlap_pc_res_prob;

    end
    
    [~,name,~] = fileparts(data_folders{i});
    pos_img = strfind(name, 'Img');
    pos_hp = strfind(name, 'HP')-1;
    name = name(pos_img:pos_hp);
    
    results.(name) = overlap_pc_res_prob_ges;
    results_overall(i,:) = overlap_pc_res_prob_ges;
    names_ges{i} = name;
end

results.pc_prob_ct = pc_prob_ct_ges;
save('pc_prob_ct.mat','results');

%% Plot
prob_fiber = num2cell(prob_fiber);
pc_ct = results.pc_prob_ct;

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
data_ges = [results_overall(1,:); results_overall(3,:); results_overall(5,:); results_overall(7,:); results_overall(9,:); results_overall(11,:)];
label_analysis = {'Img: A - Sc: A' 'Img: A - Sc: A [rec excl]' 'Img: A - Sc: C' 'Img: A - Sc: C [rec excl]' 'Img: C - Sc: C' 'Img: C - Sc: C [rec excl]'};

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

f = gcf;
set(f,'Units','centimeters','Position',[10 10 16 15]);

% Grafik speichern
name = 'prob_dist_poster';
path = pwd;
rez = 300; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
close all

%% T1 Analyses
data_ges = [results_overall(8,:); results_overall(2,:); results_overall(6,:); results_overall(10,:)];
label_analysis = {'Img: T1 C - Sc: C' 'Img: Clin A [adj] - Sc: A' 'Img: Clin A [adj] - Sc: C' 'Img: Clin C [adj] - Sc: C'};

figure(2)
for i = 1:length(label_analysis)
    subplot(length(label_analysis),1,i)
    data_plot = data_ges(i,:);
    bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
    hold on
    bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
    
    ylabel('overlap [%]','Fontsize',fsg)
        
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

f = gcf;
set(f,'Units','centimeters','Position',[10 8 8 15]);

% Grafik speichern
name = 'prob_dist_t1';
rez = 300; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
close all

%% Recovered only
data_ges = results_overall(4,:);
label_analysis = {'Img: Clin A [rec only] - Sc: A'};

figure(3)
data_plot = data_ges;
bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
hold on
bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');

ylabel('overlap [%]','Fontsize',fsg)

step_acc = 20;
yticks_acc = 0:step_acc:y_val;
set(gca,'YTick',yticks_acc,'Fontsize',fsu,'ylim', [yticks_acc(1),yticks_acc(end)])

xticks = prob_fiber;
set(gca,'XTickLabel',xticks,'Fontsize',fsu);
zu = 0.7;
set(gca,'xlim', [1-zu,10+zu]);
title(label_analysis);
line([2.5 2.5],[2.5 y_val],'Linewidth',lw,'Color','k','LineStyle','--');
box off

f = gcf;
set(f,'Units','centimeters','Position',[10 10 10 7]);

% Grafik speichern
name = 'prob_dist_rec_only';
rez = 300; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
close all

display('done');
toc
