% Conducts VLSM Analyses with NiiStat
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

        cf = [currentFolder '\Maps_norm_' tp_a '\NII_' mod_a '_k_mat'];
        cd (cf);

        data_folders = dir('DataSheet*');
        data_folders_neu = cell(length(data_folders),1);
        for w=1:length(data_folders)    
            data_folders_neu{w} = data_folders(w).name;    
        end
        data_folders = data_folders_neu;
        clear data_folders_neu    
        
        for p = 1:length(data_folders)
            
            curr_name = data_folders{p};
            roiIndices = 0;
            modalityIndices = 1;
            numPermute = 2000;
            pThresh = 0.05;            
            [num,txt,raw] = xlsread(curr_name,'Data (2)');
            vpn = numel(num);
            minOverlap = round(vpn*0.15);
            
            if vpn < 30
                display('not enough participants');
            else                
                nii_stat(curr_name,roiIndices,modalityIndices,numPermute,pThresh,minOverlap);
                %nii_stat(xlsname, roiIndices, modalityIndices,numPermute, pThresh, minOverlap, regressBehav, maskName, GrayMatterConnectivityOnly, deSkew, doTFCE, doSVM)
            end
            
        end
        
    end
end

%% Results ordnen
cd F:\ACLM
currentFolder = pwd;
ta_f = 'F:\ACLM\Matlab_Skripte\VLSM_Results_NiiStat';

for i = 1:length(tp)
    
    for j = 1:length(mod)
    
        tp_a = tp{i}; % Zeitpunkt aktuell
        mod_a = mod{j}; % Modality

        cf = [currentFolder '\Maps_norm_' tp_a '\NII_' mod_a '_k_mat'];
        cd (cf);
        
        data_folders = regexp(genpath(cf),['[^;]*'],'match');
        data_folders = data_folders(2:end);
        
        for p = 1:length(data_folders)
            curr_f = data_folders{p};
            cd(curr_f);
            
            z_thresh_f = dir('threshZlesion*');
            z_thresh_f = z_thresh_f.name;
            
            copyfile(z_thresh_f,ta_f);
        end
        
    end
    
end

cd F:\ACLM

clc
display('done');
