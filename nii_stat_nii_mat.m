% Wandelt nifit Dateien in mat Files für NiiStat um

clear all
close all
clc

cd F:\ACLM
currentFolder = pwd;

% Img Data
tp = {'akut' 'chron'}; % Zeitpunkte
mod = {'clin' 'T1'};

for i = 1:length(tp)
    
    for j = 1:length(mod)
    
        tp_a = tp{i}; % Zeitpunkt aktuell
        mod_a = mod{j}; % Modality

        cf = [currentFolder '\Maps_norm_' tp_a '\NII_' mod_a '_k'];
        cz = [currentFolder '\Maps_norm_' tp_a '\NII_' mod_a '_k_mat'];
        mkdir(cz);
        cd (cf);

        data_folders = dir('kACLM*.nii');
        data_folders_neu = cell(length(data_folders),1);
        for w=1:length(data_folders)    
            data_folders_neu{w} = data_folders(w).name;    
        end
        data_folders = data_folders_neu;
        clear data_folders_neu    
        
        for p = 1:length(data_folders)
            
            curr_data = data_folders{p};
            nii_nii2mat(curr_data,1);
            
            [pathstr,img_name,ext] = fileparts(curr_data);            
            img_mat = [img_name '.mat'];
            
            movefile(img_mat,cz);
            
        end        
        
    end
    
end
