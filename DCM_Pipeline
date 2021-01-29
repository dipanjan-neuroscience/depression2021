%initialize spm
clear all

SPM_PATH = 'D:\1_DEPRESSION_WORK\Code\spm12';
addpath(SPM_PATH)



spm('Defaults','fMRI');
spm_jobman('initcfg');

%% extract voi

%subNames= {'sub-01','sub-02', 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07',  'sub-09', 'sub-10','sub-11','sub-12', 'sub-13', 'sub-14','sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20','sub-21','sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-30', 'sub-31', 'sub-33', 'sub-34', 'sub-35', 'sub-36', 'sub-37', 'sub-38', 'sub-39', 'sub-40','sub-41','sub-42', 'sub-43', 'sub-44', 'sub-45', 'sub-47', 'sub-48', 'sub-49','sub-51','sub-52', 'sub-53', 'sub-54', 'sub-56', 'sub-57', 'sub-58', 'sub-59', 'sub-60','sub-62', 'sub-63', 'sub-64', 'sub-66', 'sub-67', 'sub-68', 'sub-69', 'sub-70','sub-71', 'sub-72'};
%subNames= {'sub-01', 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08', 'sub-09', 'sub-10','sub-11','sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20','sub-21','sub-22', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29'};
subNames= {'sub-01', 'sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08', 'sub-09', 'sub-10','sub-11','sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-19', 'sub-21','sub-22', 'sub-23', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28'};
%subNames= {'sub-01', 'sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08', 'sub-09','sub-11','sub-12', 'sub-13', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18', 'sub-19', 'sub-20','sub-21','sub-22', 'sub-24', 'sub-25', 'sub-26', 'sub-27', 'sub-28', 'sub-29'};

% maskNames={ 'lBroca.nii','lThal.nii','rBroca.nii','rThal.nii'};
% voiNames= {'lBroca','lThal','rBroca','rThal'};

% maskNames={'lsThal.nii','rsThal.nii'};
% voiNames= {'lsThal','rsThal'};

% maskNames={'lp24cd.nii','lp24rv.nii','lp32.nii','lsg32.nii','rp24cd.nii','rp24rv.nii','rp32.nii','rsg32.nii'};
% voiNames= {'lp24cd','lp24rv','lp32','lsg32','rp24cd','rp24rv','rp32','rsg32'};

% maskNames={ 'laIns.nii','lpIns.nii','raIns.nii','rpIns.nii'};
% voiNames= {'laIns','lpIns','raIns','rpIns'};


maskNames={'lFp1.nii','lSSC.nii','lV1.nii','lA1.nii','lSMA.nii','lMC.nii','laIns.nii','lpIns.nii','lThal.nii','lBroca.nii','rFp1.nii','rSSC.nii','rV1.nii','rA1.nii','rSMA.nii','rMC.nii','raIns.nii','rpIns.nii','rThal.nii','rBroca.nii'};
voiNames={'lFp1','lSSC','lV1','lA1','lSMA','lMC','laIns','lpIns','lThal','lBroca','rFp1','rSSC','rV1','rA1','rSMA','rMC','raIns','rpIns','rThal','rBroca'};


% subNames= {'sub-10'}; 
% maskNames={'sub10_lV1.nii','sub10_rV1.nii'};
% voiNames={'lV1','rV1'};
% 
% 
% subNames= {'sub-14'}; 
% maskNames={'sub14_lfp1.nii','sub14_rfp1.nii'};
% voiNames={'lFp1','rFp1'};



for sI = 1: length(subNames)

spm_dir = fullfile('D:\1_DEPRESSION_WORK\Programs\Input\DepRestingState\GLM2_POST', subNames{sI},'SPM.mat');

for vI = 1: length(voiNames)

mask_dir= fullfile ('D:\1_DEPRESSION_WORK\Code\Masks\PC', maskNames{vI});


clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat = cellstr(spm_dir); % directory to spm.mat file
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1; % Session index
matlabbatch{1}.spm.util.voi.name = voiNames{vI}; % name you want to give 
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(mask_dir); % directory to mask image
%matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask.mtype = 0; % inclusion
matlabbatch{1}.spm.util.voi.expression = 'i1';
spm_jobman('run',matlabbatch);



end


end


%% define DCMs

voiNamesL={'VOI_lSMA_1.mat', 'VOI_lMC_1.mat'};
voiNamesR={'VOI_rSMA_1.mat', 'VOI_rMC_1.mat'};   

% voiNamesLintero5={'VOI_lp24cd_1.mat', 'VOI_lp24rv_1.mat','VOI_lp32_1.mat', 'VOI_lsg32_1.mat','VOI_lICC_1.mat'};
% voiNamesRintero5={'VOI_rp24cd_1.mat', 'VOI_rp24rv_1.mat','VOI_rp32_1.mat', 'VOI_rsg32_1.mat','VOI_rICC_1.mat'}; 

voiNamesL={'VOI_laIns_1.mat','VOI_lpIns_1.mat'};
voiNamesR={'VOI_raIns_1.mat','VOI_rpIns_1.mat'}; 

voiNamesL={'VOI_lFp1_1.mat','VOI_lSSC_1.mat','VOI_lV1_1.mat','VOI_lA1_1.mat'};
voiNamesR={'VOI_rFp1_1.mat','VOI_rSSC_1.mat','VOI_rV1_1.mat','VOI_rA1_1.mat'};




for sI = 1: length(subNames)
%cd( fullfile('D:\1_DEPRESSION_WORK\Programs\Input\RestingDepHealthy\GLM2',subNames{sI}));
cd(fullfile('D:\1_DEPRESSION_WORK\Programs\Input\DepRestingState\GLM2_POST', subNames{sI}));
model_name = 'L_Mot';

xY         = voiNamesL;

SPM        = 'SPM.mat';

n   = 4;    % number of regions

nu  = 1;    % number of inputs. For DCM for CSD we have one input: null

TR  = 2.5;    % volume repetition time (seconds)

TE  = 0.035; % echo time (seconds)

 

% Connectivity matrices

a  = ones(n,n);

% a(1,4)=0;
% a(4,1)=0;

b  = zeros(n,n,nu);

c  = zeros(n,nu);

d  = zeros(n,n,0);

 

% Specify DCM

s = struct();

s.name       = model_name;

s.u          = [];

s.delays     = repmat(TR/2, 1, n)';

s.TE         = TE;

s.nonlinear  = false;

s.two_state  = false;

s.stochastic = false;

s.centre     = false;

s.induced    = 1;       % indicates DCM for CSD

s.a          = a;

s.b          = b;

s.c          = c;

s.d          = d;

 

DCM = spm_dcm_specify(SPM,xY,s);



end

clear DCM

for sI = 1: length(subNames)
%cd( fullfile('D:\1_DEPRESSION_WORK\Programs\Input\RestingDepHealthy\GLM2',subNames{sI}));
cd(fullfile('D:\1_DEPRESSION_WORK\Programs\Input\DepRestingState\GLM2_POST', subNames{sI}));
model_name = 'R_Mot';

xY         = voiNamesR;

SPM        = 'SPM.mat';

n   = 4;    % number of regions

nu  = 1;    % number of inputs. For DCM for CSD we have one input: null

TR  = 2.5;    % volume repetition time (seconds)

TE  = 0.035; % echo time (seconds)

 

% Connectivity matrices

a  = ones(n,n);
% a(1,4)=0;
% a(4,1)=0;

b  = zeros(n,n,nu);

c  = zeros(n,nu);

d  = zeros(n,n,0);

% Specify DCM

s = struct();

s.name       = model_name;

s.u          = [];

s.delays     = repmat(TR/2, 1, n)';

s.TE         = TE;

s.nonlinear  = false;

s.two_state  = false;

s.stochastic = false;

s.centre     = false;

s.induced    = 1;       % indicates DCM for CSD

s.a          = a;

s.b          = b;

s.c          = c;

s.d          = d;

 
DCM = spm_dcm_specify(SPM,xY,s);



end

clear DCM
%%estimate DCMs

 for h=1: length(subNames) 
   
 GCM_L_Mot(h,1) = {fullfile('D:\1_DEPRESSION_WORK\Programs\Input\DepRestingState\GLM2_POST', subNames{h},'DCM_L_Mot.mat')}; 
 
 end
  
 
 for h=1: length(subNames) 
   
 GCM_R_Mot(h,1) = {fullfile('D:\1_DEPRESSION_WORK\Programs\Input\DepRestingState\GLM2_POST', subNames{h},'DCM_R_Mot.mat')}; 
 
  end

use_parfor = true ;
GCM_L_Mot = spm_dcm_fit(GCM_L_Mot);
save('GCM_L_Mot.mat','GCM_L_Mot');


GCM_R_Mot = spm_dcm_fit(GCM_R_Mot);
save('GCM_R_Mot.mat','GCM_R_Mot');


%%PEB
% load GCM_L_Mot.mat 
% load GCM_R_Mot.mat 

load GCM_L_Mot.mat 
load GCM_R_Mot.mat 


load BDIAgeSexPost.mat

%DCM for fMRI diagnostics
spm_dcm_fmri_check (GCM_L_Mot)

BDIAgeSex(:,1)=BDIAgeSex(:,1)-mean(BDIAgeSex(:,1));

BDIAgeSex(:,2)=BDIAgeSex(:,2)-mean(BDIAgeSex(:,2));

M   = struct();
M.Q = 'all'; 

% Specify design matrix for N subjects. It should start with a constant column
M.X = horzcat(ones(25,1),BDIAgeSex);

% Choose field
field = {'A'};

% Estimate model

PEB_L_MotTrt    = spm_dcm_peb(GCM_L_Mot,M,field);

save('PEB_L_MotTrt.mat','PEB_L_MotTrt'); 

M   = struct();
M.Q = 'all'; 

% Specify design matrix for N subjects. It should start with a constant column
M.X = horzcat(ones(25,1),BDIAgeSex);

% Choose field
field = {'A'};

% Estimate model

PEB_R_MotTrt    = spm_dcm_peb(GCM_R_Mot,M,field);

save('PEB_R_MotTrt.mat','PEB_R_MotTrt'); 

%%BMR & BMA
clear
filenames={'PEB_L_ExtTrt.mat', 'PEB_R_ExtTrt.mat', 'GCM_L_Ext.mat', 'GCM_R_Ext.mat'};

  
for kk = 1:numel(filenames)
    load(filenames{kk})
end

BMA_L_ExtTrt=spm_dcm_peb_bmc(PEB_L_ExtTrt);
save('BMA_L_ExtTrt.mat','BMA_L_ExtTrt');
spm_dcm_peb_review(BMA_L_ExtTrt,GCM_L_Ext);

BMA_R_ExtTrt=spm_dcm_peb_bmc(PEB_R_ExtTrt);
save('BMA_R_ExtTrt.mat','BMA_R_ExtTrt');
spm_dcm_peb_review(BMA_R_ExtTrt,GCM_R_Ext);

%%leave one out cross validation 

%clear 

load BDIAgeSex.mat

BDIAgeSex(:,1)=BDIAgeSex(:,1)-mean(BDIAgeSex(:,1));

BDIAgeSex(:,2)=BDIAgeSex(:,2)-mean(BDIAgeSex(:,2));

M   = struct();
M.Q = 'all'; 

% Specify design matrix for N subjects. It should start with a constant column
M.X = horzcat(ones(63,1),BDIAgeSex);

[qE,qC,Q]=spm_dcm_loo(GCM_L_ExtS,M,{'A'});% (to,from)

