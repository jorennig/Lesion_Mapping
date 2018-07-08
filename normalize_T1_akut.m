clear all
clc

cd F:\ACLM;
currentFolder = pwd;

data = xlsread('F:\ACLM\ACLM_gesamt.xlsx');
data = data(:,end);
data(isnan(data)) = [];

tp_a = 'akut'; % Zeitpunkt aktuell
cf = [pwd '\Laesionen_' tp_a];       
cd (cf);

data_folders = dir(cf);
data_folders_neu = cell(length(data_folders),1);
for w=1:length(data_folders)    
    data_folders_neu{w} = data_folders(w).name;    
end
data_folders_neu = data_folders_neu(3:end-1);
data_folders = data_folders_neu;
clear data_folders_neu w

for j=1:length(data_folders)

    if data(j) == 0

        display('Kein T1');

    else % Normalisierung mit T1 Lesion Map
        
        data_a = data_folders{j}; % Daten aktuell
        cf_d = [cf '\' data_a];
        cd (cf_d);

        fileName_l = ls('ACLM*T1*.nii');
        fileName_img = [fileName_l(1:end-3) 'img'];
        fileName_hdr = [fileName_l(1:end-3) 'hdr'];

        t1_dir = 'T1_norm';
        mkdir(t1_dir);
        t1_f = [pwd '\' t1_dir];

        copyfile(fileName_l,t1_f)
        copyfile(fileName_img,t1_f)
        copyfile(fileName_hdr,t1_f)

        cd (t1_f); 

        lesion = fileName_l;
        anatomical = fileName_img;
        pathological = anatomical;
        Modality = 1;
        vox = [1 1 1];
        bb = [-78 -112 -50 78 76 85];
        DeleteIntermediateImages = 1;
        UseTemplateMask = 0;
        
        clinical_mrnorm(anatomical, lesion, pathological, vox, bb, DeleteIntermediateImages, UseTemplateMask, Modality);
        
        target_folder = 'F:\ACLM\Maps_norm_akut\NII_T1_org';
        fileName_norm = ls('bwsrACLM*.nii');
        
        copyfile(fileName_norm,target_folder);
    end
end

cd F:\ACLM;
display(' ');
display('Normalization done!');