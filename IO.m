function [pc,sem]=IO(a,mask,c,cs,sdeve,js,trials, varargin) 

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

obj = inputParser;
addParamValue(obj,'saveDir', 'SavedIO')
addParamValue(obj,'printText', '')%choice to print more text on each trial 

parse(obj, varargin{:})
save_dir = obj.Results.saveDir;
printText = obj.Results.printText;

if ~exist(fullfile( '.',save_dir), 'file')
    mkdir(fullfile('.', save_dir))
end


seed=sum(100*clock);
stream=RandStream('mt19937ar','seed',seed);
[ny,nx,ns]=size(a);
a=double(a);
n=ns/c;                                                                     %Number of distinct stimuli within each class
nc=length(cs);

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
[ny,nx,nsj]=size(aa);                                                       %nsj is the number of total templates after jittering

%%%%%%%%%%%%%%%%%%%%%%%%%%% Mask templates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nsj
    if mask == 0
        aa(:,:,i)=aa(:,:,i)-mean2(aa(:,:,i));                               %With no mask, just set stimuli to zero-mean (if your psychophysics stimuli all had the same mean luminance)
    else
        temp=aa(:,:,i);
        aa(:,:,i)= mask.*(aa(:,:,i)-mean(temp(mask)));                       %With a mask, set the visible stimulus areas (where mask==1) to zero mean and then apply the mask.
    end
end
   
%%%%%%%%%%%%%% Various operations to expedite computation %%%%%%%%%%%%%%%%%
                                                %Reshape the stimulus matrix for future vectorized computation

experiment_time = tic;
aa_array=reshape(aa,ny*nx,nsj);                                           %Linearize templates for speed

class=ceil(c*rand(stream,trials,1));                                %Randomly sample a class for each trial
exem=ceil(n*rand(stream,trials,1));                                 %Randomly sample an exemplar from within the class for each trial
con=ceil(nc*rand(stream,trials,1));
jit=ceil(nj*rand(stream,trials,1));
shown=(class-1)*n*nc*nj+(exem-1)*nc*nj+(jit-1)*nc+con;              %Compute the proper linear index for the sampled template for each trial

%%%%% Compute responses and sample %%%%
nsjc=nsj*nc;                                                        %Total number of templates (class x exemplars/class x jitters x contrasts)
cons=repmat(cs,[1 nsj]);                                       %Contrast uncertainty matrix, to be applied to the mean and covariance

temps_w_contrast = (imresize(aa_array, [ny*nx, nsjc], 'nearest').* repmat(cons(1,:), [ny*nx, 1])+127.5);

E= sum(temps_w_contrast.*temps_w_contrast);                                                          %The mean response matrix (dot products of single-filtered templates)
l=zeros(nsjc,trials);
for trial = 1:trials
    
%     r = (normrnd(temps_w_contrast(:,shown(trial)), sdeve));%trialsX1 array
    r = temps_w_contrast(:,shown(trial))+sdeve*randn([ny*nx 1]);
    for i=1:nsjc

        l(i,trial)=((2*r'*(temps_w_contrast(:,i))) - E(i))/(2*(sdeve^2)); %nsjcXtrials matrix -exponent taken out 
    end
    fprintf([printText, ' Trial Done: ' num2str(trial), '\n']);   
end

subs=(1:c*trials)';                                                 %This groups likelihoods from each class for each trial
subs=imresize(subs,[nsjc*trials 1],'nearest');
%log-sum-exp trick
%find constant for trick
consts = max(l);%different constant for each trial
l=exp(bsxfun(@minus, l,consts));                                                            %Need to linearize the likelihood for this accumulating to work
ll = l(:);
% suml=-(10^3) + (accumarray(subs,log(ll)));                                           %Each entry in sl is the sum of likelihoods from one class from one trial
suml=(accumarray(subs,ll));                                           %Each entry in sl is the sum of likelihoods from one class from one trial
suml = reshape(repmat(consts, [c 1]), [1 c*trials])'+log(suml);
sl=reshape(suml,c,trials);                                          %Reshape so that each trial is a column (each class a row)
[~,dec]=max(sl);                                                    %The maximum of the summed likelihoods is the decision for each trial
correct=dec==class';                                                %When the decision equals the true class we sampled, the trial is correct
pc=mean(correct);                                            %Proportion correct
sem=std(correct)/sqrt(trials);                               %Binary standard error of the mean
run_time=toc(experiment_time);
save(fullfile( '.',save_dir, 'run_time') ,'run_time');
save(fullfile( '.',save_dir, 'pc') ,'pc');
      

end