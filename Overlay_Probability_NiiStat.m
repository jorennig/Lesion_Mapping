% Produces overlay probability maps of VLSM z maps with the Jülich
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

size_prob_ges = zeros(length(ct_el),1);
pc_prob_ct_ges = zeros(length(ct_el),1);

for i = 1:length(ct_el)
    
    ct_prob = fiber_ct.img == ct_el(i);
    size_prob = sum(sum(sum(ct_prob)));
    size_prob_ges(i) = size_prob;
    pc_prob_ct_ges(i) = size_prob/size_fiber_ct*100;
    
    data_folders = dir('Zlesion*.nii');
    data_folders_neu = cell(length(data_folders),1);
    for w=1:length(data_folders)
        data_folders_neu{w} = data_folders(w).name;
    end    
    data_folders = data_folders_neu;
    clear data_folders_neu w
    
    % get index of different analyses/results
    idxclin = strfind(data_folders,'_clinHP');
    idxclin = find(~cellfun(@isempty,idxclin));

    idxt1 = strfind(data_folders,'_T1');
    idxt1 = find(~cellfun(@isempty,idxt1));

    idxrecex = strfind(data_folders,'rec_excl');
    idxrecex = find(~cellfun(@isempty,idxrecex));
    
    idxrecin = strfind(data_folders,'rec_only');
    idxrecin = find(~cellfun(@isempty,idxrecin));
    
    idxt1adj = strfind(data_folders,'adj');
    idxt1adj = find(~cellfun(@isempty,idxt1adj));
        
    analyses.an_clin = data_folders(idxclin);
    analyses.an_recex = data_folders(idxrecex);
    analyses.an_recin = data_folders(idxrecin);
    analyses.an_t1adj = data_folders(idxt1adj);
    analyses.an_t1 = data_folders(idxt1);

    an_names = fieldnames(analyses);
    
    for j = 1:numel(an_names)
        
        an_curr = analyses.(an_names{j});
        
        size_bin_res_ges = zeros(numel(an_curr),1);
        size_overlap = zeros(numel(an_curr),1);
        overlap_pc_res_prob_ges = zeros(numel(an_curr),1);

        for p = 1:numel(an_curr)

            file = an_curr{p}; % Daten aktuell
            res = load_nii(file); 

            bin_res = res.img > 0;
            size_bin_res = sum(sum(sum(bin_res)));
            %size_bin_res_ges(j) = size_bin_res;

            overlap = bin_res.*ct_prob;
            size_overlap = sum(sum(sum(overlap)));

            overlap_pc_res_prob = size_overlap/size_prob*100; % Prozent des motorischen Systems von Result betroffen
            overlap_pc_res_prob_ges(p) = overlap_pc_res_prob;
        
            name_prob = ['prob_' num2str(prob_fiber(i))];
            results.(an_names{j}).(name_prob) = overlap_pc_res_prob_ges;

        end

    end
end

results.pc_prob_ct = pc_prob_ct_ges;
save('pc_prob_ct.mat','results');
%clearvars -except results prob_fiber

%% Plot
prob_fiber = num2cell(prob_fiber);
pc_ct = results.pc_prob_ct;

data_ges = zeros(numel(an_names)*3,length(pc_ct));
data_sub = zeros(3,length(pc_ct));
idxc = 1;

% for i = 1:numel(an_names)
% 
%     for j = 1:length(prob_fiber) % Daten auslesen/umformen
% 
%         data = results.(an_names{i}).(['prob_' num2str((prob_fiber{j}))]);
%         data_sub(:,j) = data;
% 
%     end
%     
%     data_ges(idxc:idxc+size(data_sub,1)-1,:) = data_sub;
%     idxc = idxc+3;
%     
% end
% 
% bw = 1; % Barwidth
% barmap = [0.6 0.6 0.6]; % [0.7 0.7 0.7] is grey
% lw = 0.9; % Line Width
% lws = 1; % Line Width subplot
% fss = 10; % Font Size small
% fsu = 8; % Font Size subplot
% fsg = 9; % Font Size big
% fsgg = 14; % Font Size bigger
% f_eb = 1; % Size error bar
% tl = 0.02; % Tick length
% 
% img_an = {'Clin' 'T1' 'rec excl' 'rec only' 'Clin T1 adj'};
% idxc = 1;
% 
% figure(1)
% for i = 1:size(data_ges,1)
%     subplot(numel(an_names),3,i)
%     data_plot = data_ges(i,:);
%     bar(data_plot,bw,'FaceColor',barmap,'Edgecolor','k');    
%     hold on
%     bar(pc_ct,bw,'FaceColor','none','Edgecolor','k','LineStyle','--');
% 
%     if i == 1 || 4 || 7 || 10 || 13
%         ylabel([img_an{idxc} 'overlap [%]'],'Fontsize',fsg)
%         %text(img_an{idxc});
%         %text(50,offset1,img_an{idxc},'Fontsize',fsg,'HorizontalAlignment','Center');
%         idxc = idxc+1;
%     end
% 
%     if i == 1
%         title('Sc A - Img A');
%     end
%     if i == 2
%         title('Sc C - Img A');
%     end
%     if i == 3
%         title('Sc C - Img C');
%     end
% 
%     step_acc = 20;
%     yticks_acc = 0:step_acc:100;
%     set(gca,'YTick',yticks_acc,'Fontsize',fsu,'ylim', [yticks_acc(1),yticks_acc(end)])
% 
%     xticks = prob_fiber;
%     set(gca,'XTickLabel',xticks,'Fontsize',fsu);
%     zu = 0.7;
%     set(gca,'xlim', [1-zu,10+zu]);
%     line([2.5 2.5],[2.5 80],'Linewidth',lw,'Color','k','LineStyle','--');
%     box off
% end
% 
% f = gcf;
% set(f,'Units','centimeters','Position',[10 10 22 15]);
% 
% % Grafik speichern
% name = 'prob_dist';
% path = pwd;
% rez = 300; % resolution (dpi) of final graphic
% f = gcf; % f is the handle of the figure you want to export
% figpos = getpixelposition(f);
% resolution = get(0,'ScreenPixelsPerInch');
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% print(f,fullfile(path,name),'-dtiff',['-r',num2str(rez)],'-opengl') % save file 
% close all
% 
% display('done');
% toc
