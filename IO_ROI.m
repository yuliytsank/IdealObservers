function [pc,sem]=IO_ROI(a,mask,c,cs,sdeve,js,trials, varargin) 

%[pc,sem]=getFIO(fxs,fys,a,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials,ch)
%
%This function computes foveated ideal observer performance for a set of
%stimuli, a, at given foveation locations, (fxs,fys).
%
%OUTPUT VARIABLES
%pc:        Proportion correct.
%sem:       Standard error of the mean.
%
%INPUT VARIABLES
%a:         Matrix containing the stimuli of size [ny,nx,ns] where
%           ny: # pixels in y-dimension
%           nx: # pixels in x-dimension
%           ns: # number of distinct stimuli
%
%mask:      Cropping mask. If no mask, enter 0. Otherwise, the mask should
%           be a logical (0s and 1s) matrix of size [ny,nx].
%
%c:         # number of classes (e.g., a gender task would have 2)
%           An identification task with 10 faces would have c=10
%
%           *** This program assumes each class has an equal number of
%               stimuli. If this is not the case, you will have to tweak
%               this code. The structure of the image matrix, a, is such
%               that the third dimension, ns, has the first n elements
%               belonging to the first class, the next n elements to the
%               second class, and so on, such that ns=n*c ***
%
%cs:        Signal contrasts, values between 0 and 1. A vector implies
%           contrast uncertainty (e.g., [.1 .12 .14 .16]), while a
%           scalar implies a single contrast multiplier (e.g., .15).
%           This is usually the values you used for psychophysics.
%
%sdeve:     Standard deviation of external, additive white noise. This
%           is usually the external noise you used for psychophysics.
%js:        Vector of jitter amounts, in pixels. A separate template is
%           created for each jitter value in each of the four cardinal
%           directions. A value of zero implies no spatial uncertainty
%           (e.g., 0), while a vector implies spatial uncertainty
%           (multiple shifted versions of each template, e.g., [3 6 9])
%       
%trials:    Number of Monte Carlo simulation trials.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Optional Parameters%%%%%%%%%%%%%%%%%%%%%

try
    gpuArray(1);                                                           %check if this computer has a GPU
    canUseGPU= 1;

catch
    canUseGPU= 0;
end

obj = inputParser;
addParamValue(obj,'saveDir', 'SavedIO')
addParamValue(obj,'useGPU', 1)%use gpu if available by default
addParamValue(obj,'printText', '')%choice to print more text on each trial
addParamValue(obj,'maxTrials', 10000)%maximum number of trials to run at once

parse(obj, varargin{:})
save_dir = obj.Results.saveDir;
useGPU = obj.Results.useGPU;
printText = obj.Results.printText;
maxTrials = obj.Results.maxTrials;

if ~exist(fullfile( '.',save_dir), 'file')
    mkdir(fullfile('.', save_dir))
end


seed=sum(100*clock);
stream=RandStream('mt19937ar','seed',seed);
[ny,nx,ns]=size(a);
a=double(a);

if any(mask(:)~= 0)
    masks = repmat(mask, [1,1,ns]);
    a = a.*masks;
end

n=ns/c;                                                                     %Number of distinct stimuli within each class
nc=length(cs);

window_size = 30;
step_size = 10;
half_size = (window_size)/2;
map_height = ny-window_size;
map_width = nx-window_size;

x_positions = [0:step_size:map_width]+half_size;
y_positions = [0:step_size:map_height]+half_size;

% x_positions = 170;
% y_positions = 200;

pc = NaN([length(y_positions),length(x_positions)]);
padded_pc = ones([length(y_positions)+2, length(x_positions)+2])*(1/c);

% if canUseGPU && useGPU
%     a = gpuArray(a);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%% Jitter templates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if js == 0
    nj=1;
    aa=a;
else
    jits=length(js);
    nj=1+4*jits;                                                            %Number of jitter positions in each direction
    aa=zeros(ny,nx,ns*nj);                                                  %Initialize jittered templates matrix. This replaces images that fall outside the original image size with zeros.
    aa(:,:,1:nj:ns*nj)=a;                                                   %Non-jittered templates
    for i=1:jits
        aa(:,1:nx-js(i),1+0*jits+i:nj:ns*nj+i+0*jits)=a(:,1+js(i):nx,:);    %Left-jittered templates
        aa(1:ny-js(i),:,1+1*jits+i:nj:ns*nj+i+1*jits)=a(1+js(i):ny,:,:);    %Up-jittered templates
        aa(:,1+js(i):nx,1+2*jits+i:nj:ns*nj+i+2*jits)=a(:,1:nx-js(i),:);    %Right-jittered templates
        aa(1+js(i):ny,:,1+3*jits+i:nj:ns*nj+i+3*jits)=a(1:ny-js(i),:,:);    %Down-jittered templates
    end
end

[ny,nx,nsj]=size(aa);     
%%%%%%%%%%%%%%%%%%%%%%%%%%% Mask templates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nsj
%             if mask == 0
        aa(:,:,i)=aa(:,:,i)-mean2(aa(:,:,i));                               %With no mask, just set stimuli to zero-mean (if your psychophysics stimuli all had the same mean luminance)
% %             else
%                 temp=aa(:,:,i);
%                 aa(:,:,i)=mask.*(aa(:,:,i)-mean(temp(mask)));                       %With a mask, set the visible stimulus areas (where mask==1) to zero mean and then apply the mask.
%             end
end

experiment_time = tic;

trialSets = ceil(trials/maxTrials);%round up to run a multiple of maxTrials 
pc_trialSets = NaN([trialSets, length(y_positions),length(x_positions)]);%preallocate pc space for each trial set which will then be averaged 

noise = sdeve*randn([window_size*window_size maxTrials]);

for trialSet =1:trialSets

    for vert_ind = 1:length(y_positions)
        for horz_ind = 1:length(x_positions)

            vert_position = y_positions(vert_ind);
            horz_position = x_positions(horz_ind);
    %         disp(['Y: ',num2str(vert_position), ' X: ', num2str(horz_position)])

            small_aa = aa( (vert_position-half_size+1):(vert_position+half_size),(horz_position-half_size+1):(horz_position+half_size) ,:);
            [ny,nx,nsj]=size(small_aa);      

            %%%%%%%%%%%%%% Various operations to expedite computation %%%%%%%%%%%%%%%%%
                                                            %Reshape the stimulus matrix for future vectorized computation

            
            aa_array=reshape(small_aa,ny*nx,nsj);                                           %Linearize templates for speed

            class=ceil(c*rand(stream,maxTrials,1));                                %Randomly sample a class for each trial
            exem=ceil(n*rand(stream,maxTrials,1));                                 %Randomly sample an exemplar from within the class for each trial
            con=ceil(nc*rand(stream,maxTrials,1));
            jit=ceil(nj*rand(stream,maxTrials,1));
            shown=(class-1)*n*nc*nj+(exem-1)*nc*nj+(jit-1)*nc+con;              %Compute the proper linear index for the sampled template for each trial

            %%%%% Compute response distribution (MVN) parameters and sample %%%%
            nsjc=nsj*nc;                                                        %Total number of templates (class x exemplars/class x jitters x contrasts)
            cons=repmat(cs'*cs,nsj,nsj);                                       %Contrast uncertainty matrix, to be applied to the mean and covariance

            temps_w_contrast = (imresize(aa_array, [ny*nx, nsjc], 'nearest').* repmat(cons(1,:), [ny*nx, 1])+127.5);

            E= sum(temps_w_contrast.*temps_w_contrast);   
 
            l=zeros(nsjc,maxTrials);

            fprintf([printText, ' Set: ', num2str(trialSet), '/', num2str(trialSets) ,' Y: ',num2str(vert_position), '/', num2str(y_positions(end)), ' X: ', num2str(horz_position), '/', num2str(x_positions(end)), '\n'])
        %     r = (normrnd(temps_w_contrast(:,shown(trial)), sdeve));%trialsX1 array
            r = temps_w_contrast(:,shown)+noise;
            for i=1:nsjc
            %     temp1=r-repmat(mu(i,:),[trials 1]);                          %These three lines were my old school way of computing likelihoods
            %     temp2=temp1/covm;
            %     l(i,:)=exp(-.5*sum(temp2.*temp1,2))';
            %             l(i,:)=mvnpdf(r,mu(i,:),covm);                                  %Once again, Matlab made it too easy.
        %         l(i,trial)=exp(((2*r'*aa_array(:,i) - mu_s(i))/(2*(sdeve^2)))-(-10^3)); %nsjcXtrials matrix
                l(i,:)=(bsxfun(@minus, 2*sum(bsxfun(@times, r,(temps_w_contrast(:,i)))),  E(i)))/(2*(sdeve^2)); %nsjcXtrials matrix -exponent taken out 
            end
    %         fprintf([' Trial Done: ' num2str(trial), '\n']);   


            subs=(1:c*maxTrials)';                                                 %This groups likelihoods from each class for each trial
            subs=imresize(subs,[nsjc*maxTrials 1],'nearest');

            consts = max(l);%different constant for each trial
            l=exp(bsxfun(@minus, l,consts));                                                            %Need to linearize the likelihood for this accumulating to work
            ll = l(:);
            % suml=-(10^3) + (accumarray(subs,log(ll)));                                           %Each entry in sl is the sum of likelihoods from one class from one trial
            suml=(accumarray(subs,ll));                                           %Each entry in sl is the sum of likelihoods from one class from one trial
            suml = reshape(repmat(consts, [c 1]), [1 c*maxTrials])'+log(suml);
            sl=reshape(suml,c,maxTrials);                                          %Reshape so that each trial is a column (each class a row)
            [~,dec]=max(sl);                                                    %The maximum of the summed likelihoods is the decision for each trial
            correct=dec==class';                                                %When the decision equals the true class we sampled, the trial is correct
            pc_trialSets(trialSet,vert_ind, horz_ind)=mean(correct); 
            pc(vert_ind, horz_ind) = nanmean(pc_trialSets(:,vert_ind, horz_ind));                                      %Proportion correct
            sem=std(correct)/sqrt(trials);                               %Binary standard error of the mean                    %Binary standard error of the mean
            


            save(fullfile( '.',save_dir, 'pc') ,'pc');

        end
    end
end

padded_pc(2:(end-1), 2:(end-1)) = pc;
pc = padded_pc;

run_time=toc(experiment_time);
        
% save(fullfile('.', ['IntermediateVars_IO_' date]),'a' , 'mask', 'c' ,'cs', 'sdeve', 'js', ...
%     'trials', 'run_time', 'pc', 'sem', 'cons', 'mu', 'covm', 'r', 'l', 'll', 'suml', 'correct');
save(fullfile('.', ['IntermediateVars_IO_' date]),'a' , 'mask', 'c' ,'cs', 'sdeve', 'js', ...
    'trials', 'run_time', 'pc', 'sem', 'cons', 'r', 'l', 'll', 'suml', 'correct');

end