
load human_pics

% load HmHfNmNf_160
du= 1E-4;
dd = 2.4E-4;
dh = du/2;
nn = 5;
sdeve = 14;
sdevi = 150;
ps = .037;
c = 10;
% c = 2;
cs = .14;
js = 0;
bs = 10;
na = 16;

trials = 100000;

fxs = 150:100:350;
fys = fxs;

[pc, ~] = FIO_smallSet(fxs,fys,a,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'intNoiseMethod', 0, 'saveDir', 'SaveDir_IdTask');
