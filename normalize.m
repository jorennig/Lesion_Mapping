clear all
clc

cd F:\ACLM;
currentFolder = pwd;

tp = {'akut' 'chron'}; % Zeitpunkte

for i = 1:length(tp)
    
    cd F:\ACLM
    tp_a = tp{i}; % Zeitpunkt aktuell
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
        
        data_a = data_folders{j}; % Daten aktuell
        cf_d = [cf '\' data_a];
        cd (cf_d);
        
        data_files = dir(cf_d);        
        data_files_neu = cell(length(data_files),1);
        for w=1:length(data_files)    
            data_files_neu{w} = data_files(w).name;    
        end        
        data_files_neu = data_files_neu(3:end);
        data_files = data_files_neu;
        clear data_files_neu w
                
        % Check if already normalized, CT, T2FLAIR, T1
        norm_ges = zeros(length(data_files),1);
        ct_ges = zeros(length(data_files),1);
        t1_ges = zeros(length(data_files),1);
        t2_ges = zeros(length(data_files),1);
        dwi_ges = zeros(length(data_files),1);
        nii_ges = [];
        
        for p=1:length(data_files)            
            
            data_af = data_files{p}; % Daten aktuell Folder          
            
            norm_check = strfind(data_af,'bws');                        
            a = isempty(norm_check);            
            norm_ges(p) = a;            

            ct_check = strfind(data_af,'CT');                        
            a = isempty(ct_check);            
            ct_ges(p) = a;
            
            t1_check = strfind(data_af,'T1');                        
            a = isempty(t1_check);            
            t1_ges(p) = a;
            
            t2_check = strfind(data_af,'T2');                        
            a = isempty(t2_check);            
            t2_ges(p) = a;
            
            dwi_check = strfind(data_af,'DWI');                        
            a = isempty(dwi_check);            
            dwi_ges(p) = a;
            
            nii_check = strfind(data_af,'.nii');                        
            a = isempty(nii_check);            
            nii_ges(p) = a;            
            
        end
        
        check_norm = isempty(find(norm_ges==0,1));
        check_ct = isempty(find(ct_ges==0,1));
        check_t1 = isempty(find(t1_ges==0,1));
        check_t2 = isempty(find(t2_ges==0,1));
        check_dwi = isempty(find(dwi_ges==0,1));
        
        % Anzahl Lesion maps
        lm = length(data_files) - sum(nii_ges);                
       
        if check_norm == 0

            display('already done');

        else 
            if check_ct == 0

                fileName_l = ls('*.nii');
                ales = [cf_d '\' fileName_l];
                fileName_img = [fileName_l(1:end-3) 'img'];                    
                filename = [cf_d '\' fileName_img];

                vox = [1 1 1];
                bb = [-78 -112 -50 78 76 85];
                DelIntermediate = 1;
                mask = 1;

                clinical_ctnorm(filename,ales,vox,bb,DelIntermediate,mask);
                
            elseif check_t1 == 0 || check_t2 == 0 || check_dwi == 0
                
                vox = [1 1 1];
                bb = [-78 -112 -50 78 76 85];
                DeleteIntermediateImages = 1;
                UseTemplateMask = 0;                
                
                if lm == 1
                    
                    if check_t1 == 0
                        
                        fileName_l = ls('*.nii');
                        lesion = [cf_d '\' fileName_l];
                        fileName_path = [fileName_l(1:end-3) 'img'];                    
                        pathological = [cf_d '\' fileName_path];
                        anatomical = ls('*T1*.img');
                        Modality = 1;
                    
                    else
                        
                        fileName_l = ls('*.nii');
                        lesion = [cf_d '\' fileName_l];
                        pathological = ls('*T2*.img');
                        anatomical = ls('*T2*.img');                                             
                        Modality = 3;
                        
                    end
                    
                    clinical_mrnorm(anatomical, lesion, pathological, vox, bb, DeleteIntermediateImages, UseTemplateMask, Modality);
                    
                elseif lm == 2
                                        
                    for r=1:lm
                        
                        if r == 1 % Normalisierung mit T2FLAIR Lesion Map
                            fileName_l = ls('*T2*.nii');
                            lesion = [cf_d '\' fileName_l];
                            fileName_path = [fileName_l(1:end-3) 'img'];                    
                            pathological = [cf_d '\' fileName_path];
                            anatomical = ls('*T1*.img');
                            Modality = 1;

                            clinical_mrnorm(anatomical, lesion, pathological, vox, bb, DeleteIntermediateImages, UseTemplateMask, Modality);

                        elseif r == 2 % Normalisierung mit T1 Lesion Map
                            fileName_l = ls('ACLM*T1*.nii');
                            lesion = [cf_d '\' fileName_l];
                            fileName_path = [fileName_l(1:end-3) 'img'];                    
                            pathological = ls('ACLM*T1*.img');
                            anatomical = pathological;
                            Modality = 1;
                            
                            t1_dir = 'T1_norm';
                            mkdir(t1_dir);
                            t1_f = [pwd '\' t1_dir];                           
                            
                            anatomical_hdr = [anatomical(1:end-3) 'hdr'];
                            
                            copyfile(fileName_l,t1_f)
                            copyfile(anatomical,t1_f)
                            copyfile(anatomical_hdr,t1_f)
                            
                            cd (t1_f); 
                            
                            lesion = [t1_f '\' fileName_l];
                            anatomical = [t1_f '\' anatomical];
                            pathological = anatomical;
                            
                            clinical_mrnorm(anatomical, lesion, pathological, vox, bb, DeleteIntermediateImages, UseTemplateMask, Modality);
                            
                        end                        
                    end                    
                end
            end
        end
    end
end

cd F:\ACLM;
display(' ');
display('Normalization done!');