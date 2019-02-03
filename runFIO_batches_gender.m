
% load HmHfNmNf

imageFile = 'NmNf_80';

du = 1E-4;
dd = 2.4E-4;
dh = du/2;
nn = 5;
sdeve = 14;
% sdevi = 171;
ps = .037;
c = 2;
cs = .14;
js = 0;
bs = 10;
na = 16;

trials = 500000;

fxs = 250;
fys = 5:10:495;
sdevi = 160;
% a = a(:,:,61:100);
mask =1;

% priors_happy = [0,.107,.35, 1];
% priors_neutral = 1- priors_happy;
% priors = [priors_happy; priors_neutral]';

save_dir = 'SaveDir_FIO_batches_gender';

% sdevi = 121;
%happy/neutral task with all 80 faces
% FIO_batchesTask(fxs,fys,500,10,imageFile,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'loadTemps', 0, 'intNoiseMethod', 0, 'task', 'hn_all', 'saveDir', 'SavedTemps_hn_Mparams', 'calcPC', 0);
sdevi = 160;
% gender task with all 80 faces
[pc_g_n, ~] = FIO_batches(fxs,fys,500,10,imageFile,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'loadTemps', 0, 'intNoiseMethod', 0, 'saveDir', save_dir);
save(fullfile(save_dir, 'pc_g_n'), 'pc_g_n')

% sdevi = 134;
%gender task with only 40 neutral faces
% [pc_g_n, ~] = FIO_batchesTask(fxs,fys,500,10,imageFile,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'loadTemps', 0, 'intNoiseMethod', 0, 'task', 'g_n', 'saveDir', 'SavedTemps_hn_Mparams');
% save(fullfile('.','SavedTemps_hn_Mparams', 'pc_g_n'), 'pc_g_n')
% 
% %gender task with only 40 happy faces
% [pc_g_h, ~] = FIO_batchesTask(fxs,fys,500,10,imageFile,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'loadTemps', 0, 'intNoiseMethod', 0, 'task', 'g_h', 'saveDir', 'SavedTemps_hn_Mparams');
% save(fullfile('.','SavedTemps_hn_Mparams', 'pc_g_h'), 'pc_g_h')