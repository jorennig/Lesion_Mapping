% Erstellt Excel Files für VLSM Analysen ACLM
clear all
close all
clc

cd F:\ACLM
currentFolder = pwd;

% Img Data (erstellt struct 'img_data' mit Dateinamen der image mat files)
tp = {'akut' 'chron'}; % Zeitpunkte
mod = {'clin' 'T1'};

for i = 1:length(tp)
    
    for j = 1:length(mod)
    
        tp_a = tp{i}; % Zeitpunkt aktuell
        mod_a = mod{j}; % Modality

        cf = [currentFolder '\Maps_norm_' tp_a '\NII_' mod_a '_k_mat'];
        cd (cf);

        data_folders = dir('kACLM*.mat');
        data_folders_neu = cell(length(data_folders),1);
        for w=1:length(data_folders)    
            data_folders_neu{w} = data_folders(w).name;    
        end
        data_folders = data_folders_neu;
        clear data_folders_neu    
        
        img_data.(tp_a).(mod_a) = data_folders;
    end
    
end

cd F:\ACLM
clearvars -except img_data

%% VLSM Analyses
vpn = 60; % number of participants
data = xlsread('F:\ACLM\ACLM_gesamt.xlsx');
data = data(1:vpn,:);
sc_akut_clin = (data(:,21)-5)*-1;
sc_chron_clin = (data(:,27)-5)*-1;
rec_pat = data(:,26);

t1_akut = data(:,83); % Which patient has acute T1
t1_chron = data(:,47); % Which patient has chronic T1

sc_akut_t1_img_akut = sc_akut_clin;
sc_chron_t1_img_akut = sc_chron_clin;
sc_chron_t1_img_chron = sc_chron_clin;

for q = 1:length(t1_chron)
    
    if isnan(t1_akut(q)) == 1
        sc_akut_t1_img_akut(q) = NaN;
        sc_chron_t1_img_akut(q) = NaN;
    end
    
    if isnan(t1_chron(q)) == 1
        sc_chron_t1_img_chron(q) = NaN;
    end

end

sc_akut_t1_img_akut(isnan(sc_akut_t1_img_akut)) = [];
sc_chron_t1_img_akut(isnan(sc_chron_t1_img_akut)) = [];
sc_chron_t1_img_chron(isnan(sc_chron_t1_img_chron)) = [];

clear data

% Analysis all Subjects
tp = {'akut' 'chron'}; % Zeitpunkte
mod = {'clin' 'T1'};

headline = [{'ID Number'} {'HP'}];

% DataSheets clinical imaging
sc_akut = sc_akut_clin;
sc_chron = sc_chron_clin;

for i = 1:length(tp)
            
    tp_a = tp{i}; % Zeitpunkt aktuell
    mod_a = mod{1}; % Modality

    img_dat = img_data.(tp_a).(mod_a);

    for p = 1:length(img_dat) % Remove data extension
        img_dat{p} = img_dat{p}(1:end-4);
    end

    if i == 1 % for acute imaging

        dat_gesamt_sc_a = [img_dat num2cell(sc_akut)];
        dat_gesamt_sc_c = [img_dat num2cell(sc_chron)];
        
        dat_gesamt_sc_a = [headline;dat_gesamt_sc_a];
        dat_gesamt_sc_c = [headline;dat_gesamt_sc_c];
        
        data_name_sc_a = ['DataSheet_Img_' tp_a '_Sc_akut_' mod_a '.xls'];
        data_name_sc_c = ['DataSheet_Img_' tp_a '_Sc_chron_' mod_a '.xls'];

        xlswrite(data_name_sc_a,dat_gesamt_sc_a,'Data (2)');
        xlswrite(data_name_sc_c,dat_gesamt_sc_c,'Data (2)');
        
        cz = ['F:\ACLM\Maps_norm_' tp_a '\NII_clin_k_mat'];
        copyfile(data_name_sc_c,cz);
        copyfile(data_name_sc_a,cz);
        
    elseif i == 2 % for chronic imaging

        dat_gesamt_sc_c = [img_dat num2cell(sc_chron)];
        dat_gesamt_sc_c = [headline;dat_gesamt_sc_c];

        data_name_sc_c = ['DataSheet_Img_' tp_a '_Sc_chron_' mod_a '.xls'];
        xlswrite(data_name_sc_c,dat_gesamt_sc_c,'Data (2)');
        
        cz = ['F:\ACLM\Maps_norm_' tp_a '\NII_clin_k_mat'];
        copyfile(data_name_sc_c,cz);

    end
    
end

clear sc_akut sc_chron

%% DataSheets T1
% T1 akut Score akut
tp_a = tp{1}; % Imaging akut
mod_a = mod{2}; % Modality: T1
score = sc_akut_t1_img_akut; % Score 

img_dat = img_data.(tp_a).(mod_a);

for p = 1:length(img_dat) % Remove data extension
    img_dat{p} = img_dat{p}(1:end-4);
end

dat_gesamt_sc_c = [img_dat num2cell(score)];

dat_gesamt_sc_c = [headline;dat_gesamt_sc_c];

data_name_sc_c = ['DataSheet_Img_' tp_a '_Sc_akut_' mod_a '.xls'];

xlswrite(data_name_sc_c,dat_gesamt_sc_c,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_T1_k_mat';
movefile(data_name_sc_c,cz);

% T1 akut Score chron
tp_a = tp{1}; % Imaging akut
mod_a = mod{2}; % Modality: T1
score = sc_chron_t1_img_akut; % Score 

img_dat = img_data.(tp_a).(mod_a);

for p = 1:length(img_dat) % Remove data extension
    img_dat{p} = img_dat{p}(1:end-4);
end

dat_gesamt_sc_c = [img_dat num2cell(score)];

dat_gesamt_sc_c = [headline;dat_gesamt_sc_c];

data_name_sc_c = ['DataSheet_Img_' tp_a '_Sc_chron_' mod_a '.xls'];

xlswrite(data_name_sc_c,dat_gesamt_sc_c,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_T1_k_mat';
movefile(data_name_sc_c,cz);

% T1 chron Score chron
tp_a = tp{2}; % Imaging akut
mod_a = mod{2}; % Modality: T1
score = sc_chron_t1_img_chron; % Score 

img_dat = img_data.(tp_a).(mod_a);

for p = 1:length(img_dat) % Remove data extension
    img_dat{p} = img_dat{p}(1:end-4);
end

dat_gesamt_sc_c = [img_dat num2cell(score)];

dat_gesamt_sc_c = [headline;dat_gesamt_sc_c];

data_name_sc_c = ['DataSheet_Img_' tp_a '_Sc_chron_' mod_a '.xls'];

xlswrite(data_name_sc_c,dat_gesamt_sc_c,'Data (2)');

cz = 'F:\ACLM\Maps_norm_chron\NII_T1_k_mat';
movefile(data_name_sc_c,cz);

%% Clin Img Adjusted to number of T1 chron
% Img akut Score akut
data_name = 'DataSheet_Img_akut_Sc_akut_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(t1_chron)   
        if isnan(t1_chron(i)) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_adjT1'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

% Img akut Score chron
data_name = 'DataSheet_Img_akut_Sc_chron_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(t1_chron)   
        if isnan(t1_chron(i)) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_adjT1'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

% Img chron Score chron
data_name = 'DataSheet_Img_chron_Sc_chron_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(t1_chron)   
        if isnan(t1_chron(i)) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_adjT1'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_chron\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

%% Exclude Recovered Clin
% Img akut Score akut
data_name = 'DataSheet_Img_akut_Sc_akut_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(rec_pat)   
        if rec_pat(i) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_rec_excl'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

% Img akut Score akut
data_name = 'DataSheet_Img_akut_Sc_chron_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(rec_pat)   
        if rec_pat(i) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_rec_excl'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

% Img chron Score akut
data_name = 'DataSheet_Img_chron_Sc_chron_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

for i=1:length(rec_pat)   
        if rec_pat(i) == 1
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_rec_excl'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_chron\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

%% Recovered only Clin
% Img akut Score akut
data_name = 'DataSheet_Img_akut_Sc_akut_clin';
[num,txt,raw] = xlsread(data_name,'Data (2)');
txt = txt(2:end,:);
txt = txt(:,1);

rec_pat(isnan(rec_pat)) = 3;

for i=1:length(rec_pat)   
        if rec_pat(i) == 0
            txt{i} = NaN;
            num(i) = NaN;
        end
end

fh = @(x) all(isnan(x(:)));
txt(cellfun(fh, txt)) = [];
num(isnan(num)) = [];

data_vlsm = [txt num2cell(num)];
data_vlsm = [headline; data_vlsm];

name_new = [data_name '_rec_only'];

xlswrite(name_new,data_vlsm,'Data (2)');

cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
movefile([name_new '.xls'],cz);

% % Img akut Score chron
% data_name = 'DataSheet_Img_akut_Sc_chron_clin';
% [num,txt,raw] = xlsread(data_name,'Data (2)');
% txt = txt(2:end,:);
% txt = txt(:,1);
% 
% for i=1:length(rec_pat)   
%         if rec_pat(i) == 0
%             txt{i} = NaN;
%             num(i) = NaN;
%         end
% end
% 
% fh = @(x) all(isnan(x(:)));
% txt(cellfun(fh, txt)) = [];
% num(isnan(num)) = [];
% 
% data_vlsm = [txt num2cell(num)];
% data_vlsm = [headline; data_vlsm];
% 
% name_new = [data_name '_rec_only'];
% 
% xlswrite(name_new,data_vlsm,'Data (2)');
% 
% cz = 'F:\ACLM\Maps_norm_akut\NII_clin_k_mat';
% movefile([name_new '.xls'],cz);
% 
% % Img akut Score akut
% data_name = 'DataSheet_Img_chron_Sc_chron_clin';
% [num,txt,raw] = xlsread(data_name,'Data (2)');
% txt = txt(2:end,:);
% txt = txt(:,1);
% 
% for i=1:length(rec_pat)   
%         if rec_pat(i) == 0
%             txt{i} = NaN;
%             num(i) = NaN;
%         end
% end
% 
% fh = @(x) all(isnan(x(:)));
% txt(cellfun(fh, txt)) = [];
% num(isnan(num)) = [];
% 
% data_vlsm = [txt num2cell(num)];
% data_vlsm = [headline; data_vlsm];
% 
% name_new = [data_name '_rec_only'];
% 
% xlswrite(name_new,data_vlsm,'Data (2)');
% 
% cz = 'F:\ACLM\Maps_norm_chron\NII_clin_k_mat';
% movefile([name_new '.xls'],cz);

%% Delete rest
cd F:\ACLM;
data_sheet = ls('DataSheet*');
for i=1:size(data_sheet,1)
    dd = data_sheet(i,:);
    delete(dd);
end

clc
display('done');
