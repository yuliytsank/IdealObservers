%fits parameters of FIO for a 500x500 image and degrees/px parameter of
%.037

%% Forced Fixation Data to Fit Model to
forced_fixations = [28, 120, 165, 200, 280, 301, 360, 438];
forced_fix_pc = [0.685052625, 0.75896775, 0.771905125, 0.794511875, 0.784707125, 0.797532625, 0.746166, 0.718614625];
forced_fix_std_error = [0.0371833097, 0.0289628982, 0.0362316964, 0.0348166069, 0.0453317072, 0.0336991118, 0.0380966326, 0.041915842];

%% Parameters of Model That Are set Based on Psychophysics 

load('human_pics.mat');%load file of stimuli to use

ps = .037;% pixel size in degrees of visual angle
cs = .14; %contrast
sdeve = 14;%standard deviation of external noise
c = 10;%number of classes
js = 0;
bs = 10;
na = 16;
trials = 100000; %number of simulated trials
save_dir = 'SavedTemps_ParamFit';
fxs = 250; % horizontal fixation position in center of 500x500 image
fys = forced_fixations; % vertical fixation positions of model to fit to human data at same points

%% Parameters Ranges of Model To Be Fit 
params.noise = 140:2:170;
params.du = [2E-6, 3E-6];
params.dd = [1.3E-6, 1.5E-6, 1.7E-6]; 
params.nn = [4.8:.1:5.5];

%find PCs

for du_num = 1: length(params.du)
    du = params.du(du_num);
    params.results.pc.(['du' num2str(du_num)]) = NaN([ length(forced_fixations), length(params.dd), length(params.nn), length(params.noise)]);
    params.results.rms.(['du' num2str(du_num)]) = NaN(1, length(params.dd), length(params.nn), length(params.noise));
    for dd_num = 1:length(params.dd)
        dd = params.dd(dd_num);
        dh = dd/2;
        
        for nn_num = 1:length(params.nn)
            nn = params.nn(nn_num);
            for noise_num = 1:length(params.noise)
                sdevi = params.noise(noise_num);
                [pc, sem] = FIO_smallSet(fxs,fys,a,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, 'intNoiseMethod', 0, 'saveDir', save_dir);
                params.results.pc.(['du' num2str(du_num)])(:, dd_num, nn_num, noise_num) = pc.first;
                params.results.rms.(['du' num2str(du_num)])(:, dd_num, nn_num, noise_num) = norm(forced_fix_pc-pc.first')/sqrt(length(forced_fix_pc));
                save('params', 'params');
            end
            delete(fullfile(save_dir, 'faces_processed.mat'))
        end
    end

end

%find AIC for PCs

for du_num = 1: length(params.du)
        du = params.du(du_num);
        params.results.AIC.(['du' num2str(du_num)]) = NaN(1, length(params.dd), length(params.nn), length(params.noise));
        for dd_num = 1:length(params.dd)
            dd = params.dd(dd_num);
            dh = dd/2;

            for nn_num = 1:length(params.nn)
                nn = params.nn(nn_num);
                for noise_num = 1:length(params.noise)
                    noise = params.noise(noise_num);
    %                
                    pc.first = params.results.pc.(['du' num2str(du_num)])(:, dd_num, nn_num, noise_num);
                    params.results.AIC.(['du' num2str(du_num)])(:, dd_num, nn_num, noise_num) = sum((forced_fix_pc-pc.first').^2/(forced_fix_std_error).^2);
                    save('params', 'params');
                end
                
            end
        end
    
end