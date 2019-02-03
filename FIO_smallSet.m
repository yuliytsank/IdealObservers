function [pc,sem]= FIO_smallSet(fxs,fys,a,mask,ps,c,cs,sdeve,sdevi,js,bs,na,du,dd,dh,nn,trials, varargin) 
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
%fxs:       Vector containing the x-coordinates, in pixels, for foveation.
%
%fys:       Vector containing the y-coordinates, in pixels, for foveation.
%
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
%
%sdevi:     Standard deviation of internal, additive white noise. This is
%           a free parameter that modulates overall performance.
%
%js:        Vector of jitter amounts, in pixels. A separate template is
%           created for each jitter value in each of the four cardinal
%           directions. A value of zero implies no spatial uncertainty
%           (e.g., 0), while a vector implies spatial uncertainty
%           (multiple shifted versions of each template, e.g., [3 6 9])
%       
%bs:        Size, in pixels, of each eccentricity bin.
%           (I've found that anything less than 10 works well)
%
%na:        Number of angle bins (I usually use 16).
%
%du:        Eccentricity constant in the up direction.
%
%dd:        Eccentricity constant in the down direction.
%['scale', num2str(sizes(scale_ratio)), '_']
%dh:        Eccentricity constant in the horizontal directions.
%
%nn:        Steep eccentricity roll-off factor.
%
%trials:    Number of Monte Carlo simulation trials.
%
%For faces, I've found the following values to be in the ballpark for
%fitting human forced fixation performance:%sdevi ~ 100
%du ~ .0001;
%dd ~ .0001;
%dh ~ .00005;
%nn ~ 5;
%Obviously, these values are going to differ, either slightly or
%substantially, depending on your stimuli and task.

try
    gpuArray(1);                                                           %check if this computer has a GPU
    canUseGPU= 1;

catch
    canUseGPU= 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Optional Parameters%%%%%%%%%%%%%%%%%%%%%

fixations =1;
fixation = 1;

obj = inputParser;
% addParamValue(obj,'loadTemps', 1)
addParamValue(obj,'saveDir', 'SavedTemps')
%addParamValue(obj,'interFixNoise', .05)
addParamValue(obj,'noiseConst', .01)
addParamValue(obj,'intNoiseMethod', 1)%if set to 1, internal noise is not filtered 
addParamValue(obj,'saveText', '')%choice to print more text on each trial 

parse(obj, varargin{:})

saveText = obj.Results.saveText;

addParamValue(obj,'printText', saveText)%choice to print more text on each trial 

parse(obj, varargin{:})

% loadTemps = obj.Results.loadTemps;
save_dir = obj.Results.saveDir;
%interFixNoise = obj.Results.interFixNoise;
noiseConst = obj.Results.noiseConst;
intNoiseMethod = obj.Results.intNoiseMethod;
printText = obj.Results.printText;

seed=sum(100*clock);
stream=RandStream('mt19937ar','seed',seed);

% if canUseGPU
%     a = gpuArray(a);
% end

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
        aa(:,:,i)=mask.*(aa(:,:,i)-mean(temp(mask)));                       %With a mask, set the visible stimulus areas (where mask==1) to zero mean and then apply the mask.
    end
end
   
%%%%%%%%%%%%%%%%%%%%% Fourier transform templates %%%%%%%%%%%%%%%%%%%%%%%%%
ffta=fftshift(fft2(aa));
ffta=circshift(ffta,[0 0 nsj/2]);                                           %FFT2 does strange things to the order of things when the matrix is 3D

%%%%%%%% Make the foundation for the CSF in the frequency domain %%%%%%%%%%
fmax=.5/ps;                                                                 %Maximum spatial frequency in cycles per degree
[f1,f2]=freqspace([ny nx],'meshgrid');                                      %Create a meshed frequency grid (this function normalizes the max frequency at 1)
f=sqrt(f1.^2+f2.^2);                                                        %Create a 2-D frequency grid to be used for specifying the CSFs
f=f*fmax;                                                                   %Convert the normalized frequency space to our current stimulus parameters

%%%%%% Specify the CSF parameters that set max contrast sensitivity, %%%%%%
%%%%%% frequency of max sensitivity, and bandwidth                   %%%%%%
a0=1.2;
b0=0.3;
c0=5/8;

%%%%%%%%% Set the angular and eccentricity-dependent parameters %%%%%%%%%%%
nas=(na+4)/4;                                                               %Number of angular bins in each quarter of the visual field (i.e., upper-right, lower-right, upper-left, lower-left)
duh=linspace(du,dh,nas);                                                    %Eccentricity parameter is linearly interpolated from the up direction value, du, to the horizontal value,
ddh=linspace(dh,dd,nas);                                                    %dh, and then to the down value, dd, at each angular bin
d0=[duh ddh(2:nas)];                                                        %Left and right sides of the visual field are assumed symmetric
cc=180/na;                                                                  %2*cc is the size, in degrees, of each angular bin
ang=[0 cc:2*cc:180-cc 180];                                                 %Specify the boundaries of the angular bins for a hemifield

%%%%%%%%%% Create grids for specifying eccentricities and angles %%%%%%%%%%
[x,y]=meshgrid(-nx:1:nx,-ny:1:ny);                                          %Make a grid twice as big as the stimuli in each direction so that all possible eccentricities are accounted for.
% if canUseGPU
%     x = gpuArray(x);
%     y = gpuArray(y);
% end

ecc=sqrt(x.^2+y.^2);                                                        %Eccentricity grid
angs=(180/pi)*atan2(y,x);                                                   %Angular grid
angs=angs+90;
angs(:,nx:-1:1)=angs(:,nx+2:1:end);                                         %Make the up direction 0, down direction 180, and symmetric across the vertical midline

%%%%%%%%%%%%%% Various operations to expedite computation %%%%%%%%%%%%%%%%%
a=reshape(aa,ny*nx,nsj);                                                    %Reshape the stimulus matrix for future vectorized computation
nfx=length(fxs);
nfy=length(fys);
time=zeros(nfy,nfx);

%%%% The actual simulation code. Simulates each fixation, one at a time %%%
% sa = cell(nfx, nfy);
% da = cell(nfx, nfy);
% csf_values = cell(nfx, nfy);
temp_params = zeros(1, 2);

%find amount of free RAM in kbs then convert to bytes
[r,w] = unix('free');
stats = str2double(regexp(w, '[0-9]*', 'match'));

if  numel(stats) == 0 
    free_mem = 0;
else
    free_mem = stats(1)*1000;
end

saved_combined_temps.name = [];
saved_single_temps.name = [];
saved_double_temps.name = [];

if ~exist(fullfile( '.',save_dir), 'dir')
    mkdir(save_dir)
end

if exist(fullfile( '.',save_dir, [saveText, 'sa_mu.mat']), 'file')
    load(fullfile( '.',save_dir, 'temp_params.mat'))
    saved_temp_num = length(temp_params)+1;
%     compiled_temp_obj = matfile(fullfile( '.',save_dir, [saveText, 'sa_mu.mat']), 'Writable',false);
%     saved_combined_temps = whos( '-FILE',  fullfile( '.',save_dir, [saveText, 'sa_mu.mat']));
%     saved_single_temps = whos( '-FILE',  fullfile( '.',save_dir, 'single.mat'));
%     saved_double_temps = whos( '-FILE',  fullfile( '.',save_dir, 'double.mat'));
    
%     used_mem = 0; 
    
%     if loadTemps
%         for temp = 1:length(saved_single_temps)
%             if free_mem*.5>=used_mem
%                 loaded = load(fullfile( '.',save_dir, 'single'),  saved_single_temps(temp).name);
%                 temp_single.(saved_single_temps(temp).name) = loaded.(saved_single_temps(temp).name);
%                 loaded = load(fullfile( '.',save_dir, 'double'),  saved_double_temps(temp).name);
%                 temp_double.(saved_double_temps(temp).name) = loaded.(saved_double_temps(temp).name);
%             else
% 
%                 break
%             end
%             infos = whos('temp_single');
%             used_mem = 2*infos.bytes;
%         end
%     else
        temp_single.first = [];
        temp_double.first = [];
%     end
            
                
                
else
    saved_temp_num = 1;
    reused_temp_num=0;
    first =[];
    temp_single.first = [];
    temp_double.first = [];
%     save(fullfile( '.',save_dir, 'single'), 'temp_single', '-v6');
%     save(fullfile( '.',save_dir, 'double'), 'temp_double', '-v6');
    save(fullfile( '.',save_dir, [saveText, 'sa_mu']), 'first', '-v7.3');
    save(fullfile( '.',save_dir, [saveText, 'da_temp']), 'first', '-v7.3');
end

compiled_temp_obj = matfile(fullfile( '.',save_dir, [saveText, 'sa_mu.mat']), 'Writable',false);

experiment_time = tic;
% 
% sa_array_all = NaN([ ny*nx, nsj*fixations]);
% sa_array_all_unchanged = NaN([ ny*nx, nsj*fixations]);
% da_array_all = NaN([ ny*nx, nsj*fixations]);
% 
% sa_array_dot = NaN([nsj, nsj*fixations]);
% da_array_dot = NaN([nsj, nsj*fixations]);


   
pc.(['fix' num2str(fixation)])=zeros(nfy,nfx);
sem.(['fix' num2str(fixation)])=zeros(nfy,nfx);

% sa_array_dot_prev = sa_array_dot(:,1:(nsj*(fixation-1)));
% sa_array_all_prev = sa_array_all(:, 1:(nsj*(fixation-1)));
% sa_array_all_prev_unchanged = sa_array_all_unchanged(:, 1:(nsj*(fixation-1)));
% 
% da_array_all_prev = da_array_all(:, 1:(nsj*(fixation-1)));

for fx=1:nfx
    for fy=1:nfy
        fprintf([printText,' fix' num2str(fixation) ' X: ' int2str(fxs(fx)) '/' int2str(fxs(end)) ', Y: ' int2str(fys(fy)) ...
            '/' int2str(fys(end))]);
        startTime = tic;
        
        if isprop(compiled_temp_obj, ['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]) %check if combined template is already saved
            loaded = load(fullfile( '.',save_dir, [saveText, 'sa_mu']), ['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            mu = loaded.(['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            loaded = load(fullfile( '.',save_dir, [saveText, 'da_temp']), ['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            temp = loaded.(['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            
            fprintf(' Temp Was Saved ');

        else

            %%% Select the area of the large eccentricity grid relevant to %%%%
            %%% the given foveation point                                  %%%%
            ecc0=ecc(ny+2-fys(fy):2*ny+1-fys(fy),nx+2-fxs(fx):2*nx+1-fxs(fx));  %Relevant eccentricities
%             ecc0=repmat(ecc0,[1 1 nsj]);                                        %Replicate the eccentricity map for quick computation
            angs0=angs(ny+2-fys(fy):2*ny+1-fys(fy),nx+2-fxs(fx):2*nx+1-...
                fxs(fx));                                                       %Relevant angles
%             angs0=repmat(angs0,[1 1 nsj]);                                      %Replicate the angular map for quick computation

            %%%%%%%%%% Create double and single-filtered templates %%%%%%%%%%%%
            sa =zeros(ny,nx,nsj);                                                %Single filtered
            da =zeros(ny,nx,nsj);
%             csf_values{fx, fy} = zeros(ny,nx,nsj);
            %Double filtered
            for i=1:na/2+1                                                      %Each combination of eccentricity and angle has a unique CSF
                emax=max(ecc0(angs0>=ang(i) & angs0<=ang(i+1)));                %Find the maximum eccentricity
                ebins=0:bs:emax-bs;                                             %Specify the boundaries of the eccentricity bins
                nb=length(ebins);                                               %Number of eccentricity bins depends on the max eccentricity and the bin size
                for j=1:nb
                    temp_name = ['temp' num2str(num2hex(d0(i))) '_' num2str(num2hex(ebins(j)))]; %name of template based on angle and bin number
                    %if template already exits then take it from template table
                    if isfield(temp_single, temp_name)
%                     if any(strcmp(temp_name, {saved_single_temps.name})) %check if combined template is already saved
                        
                                         
                        if i == length(ang)
                            ang2 = ang(i)+cc;
                        else
                            ang2 = ang(i+1);
                        end
                        if j == length(ebins)
                            ebins2 = ebins(j)+bs;
                        else
                            ebins2 = ebins(j+1);
                        end

                        pts_angs = [ang(i), ang(i), ang2, ang2];
                        pts_ebins = [ebins(j), ebins2, ebins(j), ebins2];

                        pts.x = round(sind(pts_angs).*(pts_ebins)+fxs(fx));
                        pts.y = round(fys(fy)-cosd(pts_angs).*(pts_ebins));

                        midx = round(sind(ang(i)+(ang2-ang(i))/2)*...              %Find midpoint of current bin
                            (ebins(j)+(ebins2-ebins(j))/2)+fxs(fx));
                        midy = round(fys(fy)-cosd(ang(i)+(ang2-ang(i))/2)*...
                            (ebins(j)+(ebins2-ebins(j))/2));

                        box_len_x = max(abs(midx-pts.x))+3;  %   
                        box_len_y = max(abs(midy-pts.y))+3;
                        
                        for side = 1:2%1 is right side, and 2 is left side
                    
                            if side == 2%if left side
                                midx = fxs(fx)-(midx-fxs(fx));
                            end

                            if (midx>length(angs0)) || (midx<1) || (midy>length(angs0)) || (midy<1)
                                continue
                            end

                            if midx+box_len_x>length(angs0)                              %Check if bin is on the edge of the image 
                                midx = length(angs0)-box_len_x;
                            elseif midx<=box_len_x
                                midx = box_len_x+1;
                            end
                            if midy+box_len_y>length(angs0)
                                midy = length(angs0)-box_len_y;
                            elseif midy<=box_len_y
                                midy = box_len_y+1;
                            end
                            
%                             if isfield(temp_single, temp_name)
                                temp1 = temp_single.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
                                temp2 = temp_double.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
%                             else
%                                 temp = load(fullfile( '.',save_dir, 'single'),  temp_name);
%                                 temp1 = temp.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
% %                                 tempObject = matfile(fullfile( '.',save_dir, 'single'));
% %                                 temp1 = tempObject.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
%                                 temp = load(fullfile( '.',save_dir, 'double'),  temp_name);
%                                 temp2 = temp.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
% %                                 tempObject = matfile(fullfile( '.',save_dir, 'double'));
% %                                 temp2 = tempObject.(temp_name)((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
% % %                                 if free_mem*.5>used_mem
% %                                     infos = whos('temp_single');
% %                                     used_mem = 2*infos.bytes;
% %                                     temp_single.(temp_name) = temp1;
% %                                     temp_double.(temp_name) = temp2;
% %                                 end
% %                                 clear temp
%                             end

                            small_angs0 = angs0((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
                            small_ecc0 = ecc0((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
            %                 small_ffta = ffta((midy-box_len+1):(midy+box_len), (midx-box_len+1):(midx+box_len), :);
%                             small_ffta = fftshift(fft2(aa((midy-box_len+1):(midy+box_len), (midx-box_len+1):(midx+box_len), :)));
%                             small_ffta=circshift(small_ffta,[0 0 nsj/2]);

                            small_sa = sa((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);
                            small_da = da((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :);

                            small_sa(bsxfun(@ge,small_angs0,repmat(ang(i), [1,1,nsj])) ...
                                & bsxfun(@le,small_angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                & bsxfun(@ge,small_ecc0,repmat(ebins(j), [1,1,nsj])))=...
                                temp1(bsxfun(@ge,small_angs0,repmat(ang(i), [1,1,nsj])) ...
                                & bsxfun(@le,small_angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                & bsxfun(@ge,small_ecc0,repmat(ebins(j), [1,1,nsj]))); %Select the current bin's pixels from temp and place in the final single-filtered templates' corresponding spatial region

                            small_da(bsxfun(@ge,small_angs0,repmat(ang(i), [1,1,nsj])) ...
                                & bsxfun(@le,small_angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                & bsxfun(@ge,small_ecc0,repmat(ebins(j), [1,1,nsj])))=...
                                temp2(bsxfun(@ge,small_angs0,repmat(ang(i), [1,1,nsj])) ...
                                & bsxfun(@le,small_angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                & bsxfun(@ge,small_ecc0,repmat(ebins(j), [1,1,nsj])));

                            sa((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :)=small_sa;
                            da((midy-box_len_y+1):(midy+box_len_y), (midx-box_len_x+1):(midx+box_len_x), :)=small_da;
                        end
                        
                    
                        reused_temp_num = reused_temp_num+1;

                    else
                        csf=(c0*(f.^a0).*exp(-b0*f-d0(i)*...
                        (ebins(j)*ps)^nn*f));
                        %compute singe and double filtered templates 
                        temp_single.(temp_name)=ifft2(ifftshift(bsxfun(@times,ffta,csf)),'symmetric');               %Filter all stimuli by the current bin's CSF and transform back to the spatial domain
                        temp_double.(temp_name)=ifft2(ifftshift(bsxfun(@times,ffta,csf.^2)),'symmetric');            %Filter all stimuli by the square of the current bin's CSF and transform back to the spatial domain
                        %save templates and csf value
%                         infos = whos('temp_single');
%                         used_mem = 2*infos.bytes;
%                         temp_single.(temp_name) = temp1;
%                         temp_double.(temp_name) = temp2;
                        temp_params(saved_temp_num, :) = [ang(i) ebins(j)];
%                         save(fullfile( '.',save_dir, 'single') ,'-append', '-Struct' ,'temp_single',temp_name, '-v6');
%                         save(fullfile( '.',save_dir, 'double') ,'-append', '-Struct' ,'temp_double',temp_name, '-v6');
%                         saved_temp_num = saved_temp_num+1;
%                         if free_mem*.5<=used_mem
%                             temp_single = rmfield(temp_single,temp_name);
%                             temp_double = rmfield(temp_double,temp_name);
%                         end

                        %update list of saved single filtered templates -
                        %double filtered should be the same list
%                         saved_single_temps = whos( '-FILE',  fullfile( '.',save_dir, 'single.mat'));
                    
                        %Specify the CSF (convert eccentricity to visual angle)

                        sa(bsxfun(@ge,angs0,repmat(ang(i), [1,1,nsj])) ...
                                    & bsxfun(@le,angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                    & bsxfun(@ge,ecc0,repmat(ebins(j), [1,1,nsj])))=...
                            temp_single.(temp_name)(bsxfun(@ge,angs0,repmat(ang(i), [1,1,nsj])) ...
                                    & bsxfun(@le,angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                    & bsxfun(@ge,ecc0,repmat(ebins(j), [1,1,nsj]))); %Select the current bin's pixels from temp and place in the final single-filtered templates' corresponding spatial region

                        da(bsxfun(@ge,angs0,repmat(ang(i), [1,1,nsj])) ...
                                    & bsxfun(@le,angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                    & bsxfun(@ge,ecc0,repmat(ebins(j), [1,1,nsj])))=...
                            temp_double.(temp_name)(bsxfun(@ge,angs0,repmat(ang(i), [1,1,nsj])) ...
                                    & bsxfun(@le,angs0,repmat(ang(i+1), [1,1,nsj])) ...
                                    & bsxfun(@ge,ecc0,repmat(ebins(j), [1,1,nsj])));
    %                 csf_values{fx, fy}(angs0>=ang(i) & angs0<=ang(i+1) & ecc0>=ebins(j)) = ...
    %                     csf;
%                         temp_single = rmfield(temp_single,temp_name);
%                         temp_double = rmfield(temp_double,temp_name);
                    end
%                     clear temp1 temp2
                    
                    save(fullfile( '.',save_dir, 'temp_params'), 'temp_params', 'reused_temp_num');
                end
            end
%             sa=circshift(sa,[0 0 nsj/2]);                                       %FFT2 does strange things to the order of things when the matrix is 3D
%             da=circshift(da,[0 0 nsj/2]);

            sa_array=reshape(circshift(sa,[0 0 nsj/2]),ny*nx,nsj);                                           %Linearize templates for speed
            da_array=reshape(circshift(da,[0 0 nsj/2]),ny*nx,nsj);
            

            %%%%% Compute response distribution (MVN) parameters and sample %%%%



            mu=sa_array'*sa_array;                                                          %The mean response matrix (dot products of single-filtered templates)

    %         mu=mu+abs(min(eig(mu)))*eye(size(mu));                              %Assure the mu matrix is positive-definite
            temp=da_array'*da_array;                                                        %Response covariance, external noise part, dot product of double-filtered templates

    %         temp=temp+abs(min(eig(temp)))*eye(size(temp));                      %Again, assure covariance is positive-definite

            sa_mu.(['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]) = mu;
            da_temp.(['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]) = temp;
            save(fullfile( '.',save_dir, [saveText, 'sa_mu']) ,'-append', '-Struct' ,'sa_mu',['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            save(fullfile( '.',save_dir, [saveText, 'da_temp']) ,'-append', '-Struct' ,'da_temp',['temp' num2str(fxs(fx)) '_' num2str(fys(fy))]);
            clear sa_mu da_temp
        end
        
        class=ceil(c*rand(stream,trials,1));                                %Randomly sample a class for each trial
        exem=ceil(n*rand(stream,trials,1));                                 %Randomly sample an exemplar from within the class for each trial
        con=ceil(nc*rand(stream,trials,1));
        jit=ceil(nj*rand(stream,trials,1));
        shown=(class-1)*n*nc*nj+(exem-1)*nc*nj+(jit-1)*nc+con;              %Compute the proper linear index for the sampled template for each trial
        
        nsjc=nsj*nc;                                                        %Total number of templates (class x exemplars/class x jitters x contrasts)
        cons=repmat(cs'*cs,fixation*nsj,fixation*nsj);                                        %Contrast uncertainty matrix, to be applied to the mean and covariance
        mu=imresize(mu,[nsjc fixation*nsjc],'nearest');                              %Expand the mean matrix to make room for contrast uncertainty
        mu=cons.*mu;                                                        %Apply contrast uncertainty
        temp=cons.*imresize(temp,[fixation*nsjc fixation*nsjc],'nearest');                    %Expand external noise covariance part and apply contrast uncertainty
        
        if ~intNoiseMethod   
            mu=mu+abs(min(eig(mu)))*eye(size(mu));
            temp=temp+abs(min(eig(temp)))*eye(size(temp));
            covm=sdeve^2*temp+sdevi^2*mu;                                       %The response covariance is the external noise part + internal noise part
        else
            covm=sdeve^2*temp+sdevi^2; 
        end


        V=sqrtm(covm);                                                    %These four lines constitute the old school nuts and bolts method of sampling MVN distributions
        z=randn(stream,size(covm,1),trials);
        eps=V*z;
        r=mu(shown,:)+eps';
%         r=mvnrnd(mu(shown,:),covm);                                         %Matlab made it too easy to sample MVNs with this function

        % Compute likelihoods, sum likelihoods if necessary, and decisions %
        l=zeros(nsjc,trials);
        for i=1:nsjc
            temp1=r-repmat(mu(i,:),[trials 1]);                          %These three lines were my old school way of computing likelihoods
            temp2=temp1/covm;
            l(i,:)=exp(-.5*sum(temp2.*temp1,2))';
%             l(i,:)=mvnpdf(r,mu(i,:),covm);                                  %Once again, Matlab made it too easy.
        end
        subs=(1:c*trials)';                                                 %This groups likelihoods from each class for each trial
        subs=imresize(subs,[nsjc*trials 1],'nearest');
        ll=l(:);                                                            %Need to linearize the likelihood for this accumulating to work
        suml=accumarray(subs,ll);                                           %Each entry in sl is the sum of likelihoods from one class from one trial
        sl=reshape(suml,c,trials);                                          %Reshape so that each trial is a column (each class a row)
        [~,dec]=max(sl);                                                    %The maximum of the summed likelihoods is the decision for each trial
        correct=dec==class';                                                %When the decision equals the true class we sampled, the trial is correct
        pc.(['fix' num2str(fixation)])(fy,fx)=mean(correct);                                            %Proportion correct
        sem.(['fix' num2str(fixation)])(fy,fx)=std(correct)/sqrt(trials);                               %Binary standard error of the mean
        time(fy,fx)=toc(startTime);
        fprintf([', Time: ' num2str(time(fy,fx)), '\n']);
        experiment_run_time = toc(experiment_time);
        save(fullfile( '.',save_dir, 'experiment_run_time') ,'experiment_run_time');
        save(fullfile( '.',save_dir, 'pc') ,'pc');
    end
end

 
save(fullfile('.', ['IntermediateVars_Multiple_' date]), 'fxs' , 'fys' , 'a' , 'mask' ,'ps' , 'c' ,'cs', 'sdeve', 'sdevi', 'js' , 'bs' ,'na' ...
    ,'du', 'dd', 'dh', 'nn', 'trials', 'time', 'pc', 'sem','reused_temp_num', 'saved_temp_num', ...
    'cons', 'mu', 'covm', 'r', 'l', 'll', 'suml', 'correct', 'noiseConst', 'intNoiseMethod');