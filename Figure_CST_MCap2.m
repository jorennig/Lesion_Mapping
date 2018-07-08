% Die absoluten Werte werden einfach zusammenaddiert; also z.B.:
% Überlapp Stufe 1 Jülich::   CTS: 8 + MC: 4 + MC p: 5 = 17
% 
% Jetzt werden für die Überlapp Stufe 1 Jülich die Atlasdaten der maximal möglichen ebenfalls addiert; also z.B.:
% Überlapp Stufe 1 Jülich::   CTS: 20 + MC: 10 + MC p: 10 = 40
% 
% Jetzt wir das ganze prozentual ausgedrückt; also z.B.:
% (17*100)/40 = 42.5%    Diese Wert wird dann wieder (wie zuvor bei der alten Figur 3) dunkelgrau in die Abbildung eingetragen und zwar vor die gestrichelte Linie (wie zuvor bei der alten Figur 3); letztere jetzt jedoch mit den addierten Werten

% Produces overlay probability maps of VLSM z maps with the Jülich
% Corticospinal tract
tic
clc
clear all
close all
cd ('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM');

img.mca = load_nii('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Motor_4a_r_MNI_FR10.nii');
img.mcp = load_nii('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Motor_4p_r_MNI_FR10.nii');
img.cst = load_nii('E:\ACLM\Matlab_Skripte\VLSM_Results_NPM\Fiber_CT_R.nii');
imgs_motor = fieldnames(img);

rois.mca.bin = img.mca.img > 0;
rois.mca.size = sum(sum(sum(rois.mca.bin)));
rois.mcp.bin = img.mcp.img > 0;
rois.mcp.size = sum(sum(sum(rois.mcp.bin)));
rois.cst.bin = img.cst.img > 0;
rois.cst.size = sum(sum(sum(rois.cst.bin)));
rois_motor = fieldnames(rois);

ct = img.cst.img;
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

% Pre-allocate arrays
results_overall = zeros(numel(data_folders),length(ct_el));
names_ges = cell(length(data_folders),1); % List of result names (statistical maps)
prop_roi_ges = zeros(length(data_folders),length(rois_motor));
prop_res_ges = zeros(length(data_folders),length(rois_motor));

sum_size_prob_ges_rois = zeros(length(ct_el),numel(data_folders));
sum_size_overlap_ges_rois = zeros(length(ct_el),numel(data_folders));

thresh = [3.80 3.71 3.89 3.61 3.61 3.43]; % Statistical thersholds

% Calculate overlay
for i = 1:numel(data_folders)
        
    file = data_folders{i}; % Select VLSM results
    res = load_nii(file); % Load VLSM results

    bin_res = res.img >= thresh(i); % Apply statistical theshold & binarize thresholded map
    size_bin_res = sum(sum(sum(bin_res))); % Size of thesholded map
    
    % Pre-allocate arrays
    pc_prob_ct_ges = zeros(length(ct_el),length(rois_motor)); % 
    overlap_pc_res_prob_ges = zeros(length(ct_el),length(rois_motor));
    abs_size_prob = zeros(length(ct_el),length(rois_motor));

    size_prob_ges = zeros(length(ct_el),length(rois_motor)); % Size frequencies per ROI
    size_overlap_ges = zeros(length(ct_el),length(rois_motor)); % Size frequencies per ROI

    for j = 1:numel(rois_motor)
        
        binary_roi = rois.(rois_motor{j}).bin; % Select ROI
        size_roi = rois.(rois_motor{j}).size; % Get ROI size
%         
%         overlap_roi = binary_roi.*bin_res; % Overlap of ROI with thersholded map
%         overlap_roi_size = sum(sum(sum(overlap_roi))); % Get size of overlap
% 
%         prop_roi = overlap_roi_size/size_roi*100; % Percent of whole ROI affected
%         prop_res = overlap_roi_size/size_bin_res*100; % Percent of results on whole ROI
% 
%         prop_roi_ges(i,j) = prop_roi;
%         prop_res_ges(i,j) = prop_res;

        for k = 1:length(ct_el)

            roi_prob = img.(imgs_motor{j}).img == ct_el(k); % Select k.th probabilistic frequency
            
            size_prob = sum(sum(sum(roi_prob))); % Get size of probabilistic frequency ROI
            size_prob_ges(k,j) = size_prob;
            
            pc_prob_ct_ges(k,j) = size_prob/size_roi*100;

            overlap = bin_res.*roi_prob; % Overlay of probabilistic frequency ROI with
            size_overlap = sum(sum(sum(overlap))); % Size of overlay
            size_overlap_ges(k,j) = size_overlap;
            
            overlap_pc_res_prob = size_overlap/size_prob*100; % Prozent des CSTs von Result betroffen
            overlap_pc_res_prob_ges(k,j) = overlap_pc_res_prob;

        end
                
    end
    
    sum_size_prob_ges = sum(size_prob_ges,2);
    sum_size_overlap_ges = sum(size_overlap_ges,2);

    sum_size_prob_ges_rois(:,i) = sum_size_prob_ges;
    sum_size_overlap_ges_rois(:,i) = sum_size_overlap_ges;
   
end

prop_ges = (sum_size_overlap_ges_rois./sum_size_prob_ges_rois)*100;

%% Plot
prob_fiber = {1,2,3,4,5,6,7,8,9,10};
sum_rois = sum(sum_size_prob_ges_rois(:,1));
pc_rois = (sum_size_prob_ges_rois(:,1)/sum_rois)*100;

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
data_ges = prop_ges';
label_analysis = {'Img: A - Behav: A' 'Img: A - Behav: A [rec excl]' 'Img: A - Behav: C' 'Img: A - Behav: C [rec excl]' 'Img: C - Behav: C' 'Img: C - Behav: C [rec excl]'};

figure(1)
for i = 1:length(label_analysis)
    subplot(3,2,i)
    data_plot = data_ges(i,:);
    bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
    hold on
    bar(pc_rois,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
    
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
name = 'Figure_3_Acute_Chronic_CST_MCap2';
path = pwd;
rez = 300; % resolution (dpi) of final graphic
f = gcf; % f is the handle of the figure you want to export
figpos = getpixelposition(f);
resolution = get(0,'ScreenPixelsPerInch');
set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
print(f,fullfile(path,name),'-dbmp',['-r',num2str(rez)],'-opengl') % save file 
close all

display('done');
toc
