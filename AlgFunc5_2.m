%% AlgFunc
% Calculates the iterative de-convolution algorithm.
%
% Syntax:
% [Yrf,Yrec] = AlgFunc4(X,psf_struct,state,ver,fil,sat,param_struct)
%
% Input:
% X - the original image to be de-convoluted.
% psf_struct - includes: "psf_mat" and "psf_size".
% psf_mat (inside "psf_struct") - a matrix of PSF (Point Spread Function) seperated by depth.
% psf_size (inside "psf_struct") - a matrix of PSF (Point Spread Function) dimensions.
% state - the initial state of Yrf ('thr' - soft threshold, 'zero')
% ver (default) - the version of calculation method ('RegCons').
% fil - the filtering process permission ('on','off').
% sat - the saturaion process permission ('on','off').
% param_struct - parameters structure, including: max_dif, max_loop,
% change_delay, type, nn.
%
% Output:
% Yrf - Reflectivity function.
% Yrec - Reconstructed image.

function [Yrf,Yrec] = AlgFunc5_2(X,psf_struct,state,fil,sat,param_struct)
%% Define parameters:
if exist('param_struct','var')
    max_dif = param_struct{1};
    max_loop = param_struct{2};
    change_delay = param_struct{3};
    type = param_struct{4};
    nn = param_struct{5};
    clear param_struct;
end
% Stop criterions:
if ~exist('max_dif','var')
    max_dif = 0.0000000001; % Default: 0.05.
end
if ~exist('max_loop','var')
    max_loop = 5000; % Default: 100.
end
% 3rd version:
% Adding delay (repetition) for each change:
if ~exist('change_delay','var')
    change_delay = 50; % Takes integer values, where "0" means no delay.
end
% 3rd version (4 or 8 neighbors):
if ~exist('type','var')
    type = 1; % "1" - only same pixel, "-1" - only surrounding neighbors, "0" - both.
end
if ~exist('nn','var')
    nn = 4; % 4 or 8 neighbors.
end

%% Extracting PSF elements:
psf_mat = psf_struct{1};
psf_size = psf_struct{2};

%% After transformation, and then the dividing to areas:
Y = X;

if strcmp(state,'thr')
    %% Soft thresholding:
    % Threshold:
    thr = 0.5;
    % Thresholded image (Assuming values range [0,1]):
    % Yrf is the Reflectivity function ("real" image).
    Yrf = (Y>thr).*(Y-thr);
else % (state == 'zero')
    %% Initial reflectivity function (zero):
    Yrf = zeros(size(Y)); % Reflectivity function.
end

%% Initial recreation:
L = size(psf_mat,1)/psf_size(1);
% Yrec is the Reconstructed image (recreated original image).
Yrec = zeros(size(Yrf));
Ytemp = zeros(size(Yrf));
Lpart = floor(size(Yrf,1)/L);
for i = 1:L
    if (i == L) && (i ~= 1)
        % Covering the final border:
        Ytemp(( (i-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:) = ...
            conv2(Yrf(( (i-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:), ...
            psf_mat(((i-1)*psf_size(1)+1):(i*psf_size(1)),:),'same');
        % Cutting the relevant section:
        Yrec(((i-1)*Lpart+1):end,:) = Ytemp(((i-1)*Lpart+1):end,:);
    elseif (i == 1) && (i ~= L)
        % Covering the first border:
        Ytemp(((i-1)*Lpart+1):( i*Lpart + (psf_size(1)-1)/2 ),:) = ...
            conv2(Yrf(((i-1)*Lpart+1):( i*Lpart + (psf_size(1)-1)/2 ),:),...
            psf_mat(((i-1)*psf_size(1)+1):(i*psf_size(1)),:),'same');
        % Cutting the relevant section:
        Yrec(((i-1)*Lpart+1):(i*Lpart),:) = Ytemp(((i-1)*Lpart+1):(i*Lpart),:);
    elseif (i == 1)&&(i == L)
        % Covering the only border:
        Ytemp(1:end,:) = ...
            conv2(Yrf(1:end,:),...
            psf_mat(1:psf_size(1),:),'same');
        % Cutting the relevant section:
        Yrec(1:end,:) = Ytemp(1:end,:);
    else
        % Covering the borders:
        Ytemp(( (i-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i*Lpart + (psf_size(1)-1)/2 ),:) = ...
            conv2(Yrf(( (i-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i*Lpart + (psf_size(1)-1)/2 ),:),...
            psf_mat(((i-1)*psf_size(1)+1):(i*psf_size(1)),:),'same');
        % Cutting the relevant section:
        Yrec(((i-1)*Lpart+1):(i*Lpart),:) = Ytemp(((i-1)*Lpart+1):(i*Lpart),:);
    end
end

% Saturation process (for values above 1):
if strcmp(sat,'on') % and not 'off'
    Yrec = (Yrec < 1).*Yrec + double(Yrec >= 1);
end


%% The correction loop:
i_loop = 0;
% Differences matrix:
dif_mat = abs(Y-Yrec);
% Maximal change vector:
% [value,first_index] = max(matrix_columns);
% [max_vec_val,max_vec_ind] = max(abs(dif_mat.'));

% Make a progress plot (part 1):
iter(i_loop+1) = i_loop;
MM(i_loop+1) = max(dif_mat(:));
mm(i_loop+1) = mean(dif_mat(:));

while (i_loop < max_loop)&&(max(dif_mat(:)) > max_dif)
    % Saving the last Yrf:
    Yrf_last = Yrf;
    
    % Defining change step size (all the versions):
    
    % 3rd version:
    % Adding delay (repetition) for each change:
    i_change = floor((i_loop + change_delay)/(change_delay + 1));
    
    % ver = 'RegCons';
    % regular: 0.5->0.4->...->0.1->0.09->...:
    if (i_change < 5)
        change = 0.5-0.1*i_change;
    else
        Q = fix((i_change - 4)/9); % Quotient
        R = rem((i_change - 4),9); % Remainder
        if (R == 0)
            R = 9;
            Q = Q - 1;
        end
        change = 0.1^(Q+1)-0.1^(Q+2)*R;
    end
    
    % The preparing matrix:
    check_mat = zeros(size(Yrf));
    % The theoretical recreated matrix:
    val_mat = zeros(size(Yrf));
    % 3rd version addition:
    val_mat_avg = zeros(size(Yrf));
    val_mat_avg2 = zeros(size(Yrf));
    % Preparing / Correcting the reflectivity function (added max_dif range):
    for m = 1:size(Yrec,1)
        
        % Setting the right filter (PSF):
        if strcmp(fil,'on') % and not 'off'
            i2 = floor((m-1)/Lpart)+1;
            if (i2 > L)
                i2 = L;
            end
            filter = psf_mat(((i2-1)*psf_size(1)+1):(i2*psf_size(1)),:);
        end
        
        if strcmp(fil,'on') % and not 'off'
            for n = 1:size(Yrec,2)
                
                if (Yrec(m,n) < Y(m,n) - max_dif)&&(Yrf(m,n) < 1 - change)
                    Yrftemp = Yrf;
                    Yrftemp(m,n) = Yrf(m,n) + change;
                    
                    % 3rd version addition:
                    % [original_mat,c] = mat_crop(original,m,n);
                    
                    % val = shortconv(Y(m,n),Yrftemp,filter,m,n,sat);
                    % 3rd version addition:
                    [val,val_avg,val_avg2] = shortconv2(Y,Yrftemp,filter,m,n,sat,dif_mat);
                    % (val_avg > -1) is without restriction, the values are: [0,1].
                    if ((val < dif_mat(m,n))&&(type == 1)) || ...
                            ((val_avg > -1)&&(type == -1)) || ...
                            ((val < dif_mat(m,n))&&(val_avg > -1)&&(type == 0)) || ...
                            ((val < dif_mat(m,n))&&(val_avg2 > -1)&&(type == 0))
                        check_mat(m,n) = 1;
                        % 2nd version addition:
                        val_mat(m,n) = val;
                        % 3rd version addition:
                        val_mat_avg(m,n) = val_avg;
                        val_mat_avg2(m,n) = val_avg2;
                    end
                elseif (Yrec(m,n) > Y(m,n) + max_dif)&&(Yrf(m,n) > change)
                    Yrftemp = Yrf;
                    Yrftemp(m,n) = Yrf(m,n) - change;
                    
                    % 3rd version addition:
                    % [original_mat,c] = mat_crop(original,m,n);
                    
                    % val = shortconv(Y(m,n),Yrftemp,filter,m,n,sat);
                    % 3rd version addition:
                    [val,val_avg,val_avg2] = shortconv2(Y,Yrftemp,filter,m,n,sat,dif_mat);
                    % (val_avg > -1) is without restriction, the values are: [0,1].
                    if ((val < dif_mat(m,n))&&(type == 1)) || ...
                            ((val_avg > -1)&&(type == -1)) || ...
                            ((val < dif_mat(m,n))&&(val_avg > -1)&&(type == 0)) || ...
                            ((val < dif_mat(m,n))&&(val_avg2 > -1)&&(type == 0))
                        check_mat(m,n) = -1;
                        % 2nd version addition:
                        val_mat(m,n) = val;
                        % 3rd version addition:
                        val_mat_avg(m,n) = val_avg;
                        val_mat_avg2(m,n) = val_avg2;
                    end
                end
                
            end
        else % "fil = 'off'".
            for n = 1:size(Yrec,2)
                if (Yrec(m,n) < Y(m,n) - max_dif)&&(Yrf(m,n) < 1 - change)
                    Yrf(m,n) = Yrf(m,n) + change;
                elseif (Yrec(m,n) > Y(m,n) + max_dif)&&(Yrf(m,n) > change)
                    Yrf(m,n) = Yrf(m,n) - change;
                else
                    Yrf(m,n) = Yrf(m,n);
                end
            end
        end
        
    end
    
    if strcmp(fil,'on') % and not 'off'
        % Correcting the reflectivity function (added max_dif range):
        for m2 = 1:size(Yrec,1)
            for n2 = 1:size(Yrec,2)
                if (check_mat(m2,n2) ~= 0)
                    % 2nd version (4 neighbors):
                    % log_val = nnvalcheck2(m,n,dif_mat,check_mat,val_mat);
                    % 3rd version (4 or 8 neighbors):
                    if (type == 0)
                        log_val = nnvalcheck3(m2,n2,dif_mat,check_mat,val_mat,val_mat_avg2,type,nn);
                    else
                        log_val = nnvalcheck3(m2,n2,dif_mat,check_mat,val_mat,val_mat_avg,type,nn);
                    end
                    
                    Yrf(m2,n2) = Yrf(m2,n2) + log_val*check_mat(m2,n2)*change;
                end
            end
        end
    end
    
    % Recreate the image:
    for i3 = 1:L
        if (i3 == L) && (i3 ~= 1)
            % Covering the final border:
            Ytemp(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:) = ...
                conv2(Yrf(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:), ...
                psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
            % Cutting the relevant section:
            Yrec(((i3-1)*Lpart+1):end,:) = Ytemp(((i3-1)*Lpart+1):end,:);
        elseif (i3 == 1) && (i3 ~= L)
            % Covering the first border:
            Ytemp(((i3-1)*Lpart+1):( i3*Lpart + (psf_size(1)-1)/2 ),:) = ...
                conv2(Yrf(((i3-1)*Lpart+1):( i3*Lpart + (psf_size(1)-1)/2 ),:),...
                psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
            % Cutting the relevant section:
            Yrec(((i3-1)*Lpart+1):(i3*Lpart),:) = Ytemp(((i3-1)*Lpart+1):(i3*Lpart),:);
        elseif (i3 == 1)&&(i3 == L)
        % Covering the only border:
        Ytemp(1:end,:) = ...
            conv2(Yrf(1:end,:),...
            psf_mat(1:psf_size(1),:),'same');
        % Cutting the relevant section:
        Yrec(1:end,:) = Ytemp(1:end,:);
        else
            % Covering the borders:
            Ytemp(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i3*Lpart + (psf_size(1)-1)/2 ),:) = ...
                conv2(Yrf(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i3*Lpart + (psf_size(1)-1)/2 ),:),...
                psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
            % Cutting the relevant section:
            Yrec(((i3-1)*Lpart+1):(i3*Lpart),:) = Ytemp(((i3-1)*Lpart+1):(i3*Lpart),:);
        end
    end
    
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        Yrec = (Yrec < 1).*Yrec + double(Yrec >= 1);
    end
    
    % Differences matrix update:
    dif_mat = abs(Y-Yrec);
    
    % Advancing the loop index
    i_loop = i_loop + 1;
    % Notes when maximal iteration reached:
    if (i_loop == max_loop)
        % dif_rms = rms(Y(:)-Yrec(:));
        % dif_med = median(abs(Y(:)-Yrec(:)));
        dif_avg = mean(abs(Y(:)-Yrec(:)));
        fprintf(['Maximal iteration ',int2str(max_loop),' reached on "The correction loop"\n']);
        fprintf(['The maximal absolute difference (MAD) is: ',num2str(max(dif_mat(:))),'\n']);
        % fprintf(['The root mean square (RMS) is: ',num2str(dif_rms),'\n']);
        % fprintf(['The median absolute difference is: ',num2str(dif_med),'\n']);
        fprintf(['The average absolute difference (AAD) is: ',num2str(dif_avg),'\n']);
    end
    
    % Make a progress plot (part 2):
    iter(i_loop+1) = i_loop;
    MM(i_loop+1) = max(dif_mat(:));
    mm(i_loop+1) = mean(dif_mat(:));
    
    % Stop condition / criterion for minima:
    if mm(i_loop+1) > mm(i_loop)
        % Last backed-up matrix:
        Yrf = Yrf_last;
        
        % Recreate the image:
        for i3 = 1:L
            if (i3 == L) && (i3 ~= 1)
                % Covering the final border:
                Ytemp(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:) = ...
                    conv2(Yrf(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):end,:), ...
                    psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
                % Cutting the relevant section:
                Yrec(((i3-1)*Lpart+1):end,:) = Ytemp(((i3-1)*Lpart+1):end,:);
            elseif (i3 == 1) && (i3 ~= L)
                % Covering the first border:
                Ytemp(((i3-1)*Lpart+1):( i3*Lpart + (psf_size(1)-1)/2 ),:) = ...
                    conv2(Yrf(((i3-1)*Lpart+1):( i3*Lpart + (psf_size(1)-1)/2 ),:),...
                    psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
                % Cutting the relevant section:
                Yrec(((i3-1)*Lpart+1):(i3*Lpart),:) = Ytemp(((i3-1)*Lpart+1):(i3*Lpart),:);
            elseif (i3 == 1)&&(i3 == L)
                % Covering the only border:
                Ytemp(1:end,:) = ...
                    conv2(Yrf(1:end,:),...
                    psf_mat(1:psf_size(1),:),'same');
                % Cutting the relevant section:
                Yrec(1:end,:) = Ytemp(1:end,:);
            else
                % Covering the borders:
                Ytemp(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i3*Lpart + (psf_size(1)-1)/2 ),:) = ...
                    conv2(Yrf(( (i3-1)*Lpart+1 - (psf_size(1)-1)/2 ):( i3*Lpart + (psf_size(1)-1)/2 ),:),...
                    psf_mat(((i3-1)*psf_size(1)+1):(i3*psf_size(1)),:),'same');
                % Cutting the relevant section:
                Yrec(((i3-1)*Lpart+1):(i3*Lpart),:) = Ytemp(((i3-1)*Lpart+1):(i3*Lpart),:);
            end
        end
        
        % Saturation process (for values above 1):
        if strcmp(sat,'on') % and not 'off'
            Yrec = (Yrec < 1).*Yrec + double(Yrec >= 1);
        end
        
        % Differences matrix update:
        dif_mat = abs(Y-Yrec);
        
        % Make a progress plot update (part 3):
        iter(i_loop+1) = iter(i_loop);
        MM(i_loop+1) = max(dif_mat(:));
        mm(i_loop+1) = mean(dif_mat(:));
    end
    
end
%

% Plots the progress plot:
figure
plot(iter,MM,'-o');
title('Maximum difference vs. iterations');
xlabel('Iterations');
ylabel('Difference');
%
figure
plot(iter,mm,'-o');
title('Average difference vs. iterations');
xlabel('Iterations');
ylabel('Difference');

end


%% Nearest 4 or 8 neighbours values difference check (version 3):
% Relevant Matrices:
% Differences matrix:
% dif_mat = abs(Y-Yrec);
% The preparing matrix:
% "1" - if adding the change reduces the "dif_mat" values in [m,n].
% "-1" - if removing the change reduces the "dif_mat" values in [m,n].
% "0" - else.
% check_mat = zeros(size(Yrf));
% The theoretical recreated matrix:
% "val_mat(m,n)" is "dif_mat(m,n)" of the changed values.
% val_mat = zeros(size(Yrf));
% "log_val = 1" Means it's the highest (">=") change in its surrounding, "0" otherwise.
% log_val
function log_val = nnvalcheck3(m,n,dif_mat,check_mat,val_mat,val_mat_avg,type,nn)
% Center value:
val = val_mat(m,n);
% Neighbors:
if nn == 4 % 4 Neighbors
   nn = 0;
elseif nn == 8 % 8 Neighbors
   nn = 1;
end

% "log_val = 1" Means it's the highest (">=") change in its surrounding, "0" otherwise.
log_val = 1;

if (val_mat_avg(m,n) == 1)&&(type < 1)
    log_val = 0;
end

% Compare the change to its surroundings:
if (m > 1)&&(n > 1)&&(check_mat(m-1,n-1) ~= 0)&&(nn)
    val_ref = val_mat(m-1,n-1);
    if (dif_mat(m,n) - val < dif_mat(m-1,n-1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m-1,n-1))&&(type < 1)
        log_val = 0;
    end
end
if (m > 1)&&(check_mat(m-1,n) ~= 0)
    val_ref = val_mat(m-1,n);
    if (dif_mat(m,n) - val < dif_mat(m-1,n) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m-1,n))&&(type < 1)
        log_val = 0;
    end
end
if (n > 1)&&(check_mat(m,n-1) ~= 0)
    val_ref = val_mat(m,n-1);
    if (dif_mat(m,n) - val < dif_mat(m,n-1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m,n-1))&&(type < 1)
        log_val = 0;
    end
end

if (m < size(val_mat,1))&&(n < size(val_mat,2))&&(check_mat(m+1,n+1) ~= 0)&&(nn)
    val_ref = val_mat(m+1,n+1);
    if (dif_mat(m,n) - val < dif_mat(m+1,n+1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m+1,n+1)&&(type < 1))
        log_val = 0;
    end
end
if (m < size(val_mat,1))&&(check_mat(m+1,n) ~= 0)
    val_ref = val_mat(m+1,n);
    if (dif_mat(m,n) - val < dif_mat(m+1,n) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m+1,n))&&(type < 1)
        log_val = 0;
    end
end
if (n < size(val_mat,2))&&(check_mat(m,n+1) ~= 0)
    val_ref = val_mat(m,n+1);
    if (dif_mat(m,n) - val < dif_mat(m,n+1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m,n+1))&&(type < 1)
        log_val = 0;
    end
end

if (m > 1)&&(n < size(val_mat,2))&&(check_mat(m-1,n+1) ~= 0)&&(nn)
    val_ref = val_mat(m-1,n+1);
    if (dif_mat(m,n) - val < dif_mat(m-1,n+1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m-1,n+1))&&(type < 1)
        log_val = 0;
    end
end
if (n > 1)&&(m < size(val_mat,1))&&(check_mat(m+1,n-1) ~= 0)&&(nn)
    val_ref = val_mat(m+1,n-1);
    if (dif_mat(m,n) - val < dif_mat(m+1,n-1) - val_ref)&&(type > -1)
        log_val = 0;
    end
    if (val_mat_avg(m,n) > val_mat_avg(m+1,n-1))&&(type < 1)
        log_val = 0;
    end
end

end


%% Filtering a specific band and cut it to a surrounding filter size (version 2):
% Smaller "val" is closer to reality, where "0" is the best, and "1" is the
% worst. Calculates the effect on the specific pixel only.
% Smaller "val_avg" is closer to reality, where "0" is the best, and "1" is
% the worst. Calculates the kernel-wide surrounding effect, but doesn't calculate the effect of the bigger surrounding within these borders.
% Smaller "val_avg2" is closer to reality, where "0" is the best, and "1"
% is the worst. Normalizes "val_avg" so the central pixel will have bigger
% significance on the decision.
function [val,val_avg,val_avg2] = shortconv2(original,reflectivity,filter_in,m,n,sat,dif_mat)
half_filsize_m = 0.5*(size(filter_in,1)+1);
half_filsize_n = 0.5*(size(filter_in,2)+1);
m_max = size(reflectivity,1);
n_max = size(reflectivity,2);
original_pix = original(m,n);

% Normalization matrix:
norm_mat = mat_avg(filter_in);

% Reducing the filter size into 3x3:
% filter = filter_in((half_filsize_m-1):(half_filsize_m+1),(half_filsize_n-1):(half_filsize_n+1));
% half_filsize_m = 0.5*(size(filter,1)+1);
% half_filsize_n = 0.5*(size(filter,2)+1);

% Not reducing the filter size:
filter = filter_in;

if (m >= half_filsize_m)&&(n >= half_filsize_n)&&(m <= m_max - half_filsize_m + 1)&&(n <= n_max - half_filsize_n + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat(1:end,1:end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m >= half_filsize_m)&&(n >= half_filsize_n)&&(m <= m_max - half_filsize_m + 1) && (n > n_max - half_filsize_n + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):end) - update );
    val_avg2 = norm_mat(1:end,1:size(update,2)).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):(m+half_filsize_m-1),(n-half_filsize_n+1):end);
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m >= half_filsize_m)&&(n >= half_filsize_n)&&(n <= n_max - half_filsize_n + 1) && (m > m_max - half_filsize_m + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):end,(n-half_filsize_n+1):(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):end,(n-half_filsize_n+1):(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat(1:size(update,1),1:end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):end,(n-half_filsize_n+1):(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m >= half_filsize_m)&&(m <= m_max - half_filsize_m + 1)&&(n <= n_max - half_filsize_n + 1) && (n < half_filsize_n)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):(m+half_filsize_m-1),1:(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,end-half_filsize_n+1));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):(m+half_filsize_m-1),1:(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat(1:end,(size(norm_mat,2)-size(update,2)+1):end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):(m+half_filsize_m-1),1:(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (n >= half_filsize_n)&&(m <= m_max - half_filsize_m + 1)&&(n <= n_max - half_filsize_n + 1) && (m < half_filsize_m)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(end-half_filsize_m+1,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original(1:(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat((size(norm_mat,1)-size(update,1)+1):end,1:end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat(1:(m+half_filsize_m-1),(n-half_filsize_n+1):(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m >= half_filsize_m)&&(n >= half_filsize_n) && (m > m_max - half_filsize_m + 1)&&(n > n_max - half_filsize_n + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):end,(n-half_filsize_n+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):end,(n-half_filsize_n+1):end) - update );
    val_avg2 = norm_mat(1:size(update,1),1:size(update,2)).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):end,(n-half_filsize_n+1):end);
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m <= m_max - half_filsize_m + 1)&&(n <= n_max - half_filsize_n + 1) && (m < half_filsize_m)&&(n < half_filsize_n)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize_m-1),1:(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(end-half_filsize_m+1,end-half_filsize_n+1));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original(1:(m+half_filsize_m-1),1:(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat((size(norm_mat,1)-size(update,1)+1):end,(size(norm_mat,2)-size(update,2)+1):end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat(1:(m+half_filsize_m-1),1:(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m >= half_filsize_m)&&(n <= n_max - half_filsize_n + 1) && (n < half_filsize_n)&&(m > m_max - half_filsize_m + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize_m+1):end,1:(n+half_filsize_n-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(half_filsize_m,end-half_filsize_n+1));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original((m-half_filsize_m+1):end,1:(n+half_filsize_n-1)) - update );
    val_avg2 = norm_mat(1:size(update,1),(size(norm_mat,2)-size(update,2)+1):end).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat((m-half_filsize_m+1):end,1:(n+half_filsize_n-1));
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
elseif (m <= m_max - half_filsize_m + 1)&&(n >= half_filsize_n) && (m < half_filsize_m)&&(n > n_max - half_filsize_n + 1)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize_m-1),(n-half_filsize_n+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check (central pixel difference):
    val = abs(original_pix-update(end-half_filsize_m+1,half_filsize_n));
    
    % Value check (average surrounding pixels differences):
    val_avg = abs( original(1:(m+half_filsize_m-1),(n-half_filsize_n+1):end) - update );
    val_avg2 = norm_mat((size(norm_mat,1)-size(update,1)+1):end,1:size(update,2)).*val_avg;
    val_avg2 = sum(val_avg2(:));
    val_avg = mean(val_avg(:));
    dif_temp = dif_mat(1:(m+half_filsize_m-1),(n-half_filsize_n+1):end);
    dif_temp = mean(dif_temp(:));
    if val_avg > dif_temp
        val_avg = 1;
    end
else
    error('The filter is too big');
end

end


%% Matrix of average:
function mat_out = mat_avg(mat_in)
half_size_m = 0.5*(size(mat_in,1)+1);
half_size_n = 0.5*(size(mat_in,2)+1);
[m,n] = size(mat_in);
mat_out = zeros(m,n);

% Circles: 1,8,16,...
if m*n > 9
    mat_out(:,:) = 0.25/(m*n-9);
    mat_out((half_size_m-1):(half_size_m+1),(half_size_n-1):(half_size_n+1)) = 0.25/8;
else
    mat_out((half_size_m-1):(half_size_m+1),(half_size_n-1):(half_size_n+1)) = 0.5/8;
end
mat_out(half_size_m,half_size_n) = 0.5;
end



%% Not in use (was previously):
%% Nearest 4 neighbours values difference check (version 2):
function log_val = nnvalcheck2(m,n,dif_mat,check_mat,val_mat)
% Center value:
val = val_mat(m,n);

% "log_val = 1" Means it's the highest (">=") change in its surrounding, "0" otherwise.
log_val = 1;

% Compare the change to its surroundings:
if (m > 1)&&(check_mat(m-1,n) ~= 0)
    val_ref = val_mat(m-1,n);
    if (dif_mat(m,n) - val < dif_mat(m-1,n) - val_ref)
        log_val = 0;
    end
end
if (n > 1)&&(check_mat(m,n-1) ~= 0)
    val_ref = val_mat(m,n-1);
    if (dif_mat(m,n) - val < dif_mat(m,n-1) - val_ref)
        log_val = 0;
    end
end
if (m < size(val_mat,1))&&(check_mat(m+1,n) ~= 0)
    val_ref = val_mat(m+1,n);
    if (dif_mat(m,n) - val < dif_mat(m+1,n) - val_ref)
        log_val = 0;
    end
end
if (n < size(val_mat,2))&&(check_mat(m,n+1) ~= 0)
    val_ref = val_mat(m,n+1);
    if (dif_mat(m,n) - val < dif_mat(m,n+1) - val_ref)
        log_val = 0;
    end
end

end

%% Not in use (was previously):
%% Filtering a specific band and cut it to a surrounding filter size:
function val = shortconv(original_pix,reflectivity,filter,m,n,sat)
half_filsize = 0.5*(size(filter,2)+1);
m_max = size(reflectivity,1);
n_max = size(reflectivity,2);

if (m >= half_filsize)&&(n >= half_filsize)&&(m <= m_max - half_filsize + 1)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):(m+half_filsize-1),(n-half_filsize+1):(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,half_filsize));
elseif (m >= half_filsize)&&(n >= half_filsize)&&(m <= m_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):(m+half_filsize-1),(n-half_filsize+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,half_filsize));    
elseif (m >= half_filsize)&&(n >= half_filsize)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):end,(n-half_filsize+1):(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,half_filsize));
elseif (m >= half_filsize)&&(m <= m_max - half_filsize + 1)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):(m+half_filsize-1),1:(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,end-half_filsize+1));    
elseif (n >= half_filsize)&&(m <= m_max - half_filsize + 1)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize-1),(n-half_filsize+1):(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(end-half_filsize+1,half_filsize));    
elseif (m >= half_filsize)&&(n >= half_filsize)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):end,(n-half_filsize+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,half_filsize));    
elseif (m <= m_max - half_filsize + 1)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize-1),1:(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(end-half_filsize+1,end-half_filsize+1));    
elseif (m >= half_filsize)&&(n <= n_max - half_filsize + 1)
    % Cutting the relevant area:
    shortimage = reflectivity((m-half_filsize+1):end,1:(n+half_filsize-1));
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(half_filsize,end-half_filsize+1));    
elseif (m <= m_max - half_filsize + 1)&&(n >= half_filsize)
    % Cutting the relevant area:
    shortimage = reflectivity(1:(m+half_filsize-1),(n-half_filsize+1):end);
    % Areal convolution:
    update = conv2(shortimage,filter,'same');
    % Saturation process (for values above 1):
    if strcmp(sat,'on') % and not 'off'
        update = (update < 1).*update + double(update >= 1);
    end
    % Value check:
    val = abs(original_pix-update(end-half_filsize+1,half_filsize));    
else
    error('The filter is too big');
end

end

%% Not in use:
%% Crops the nearest neighbors (3x3):
% original_mat(1,n_c),original_mat(m_c,1),original_mat(m_c,n_c),...
% original_mat(m_c,end),original_mat(end,n_c)
function [original_mat,c] = mat_crop(original,m,n)
[m_max,n_max] = size(original);
m_c = 2;
n_c = 2;

if (m>1)&&(m<m_max)&&(n>1)&&(n<n_max)
    original_mat = original((m-1):(m+1),(n-1):(n+1));
    
elseif (m>1)&&(n>1)&&(n<n_max) && (m>m_max)
    original_mat = original((m-1):m,(n-1):(n+1));
elseif (m>1)&&(m<m_max)&&(n>1) && (n>n_max)
    original_mat = original((m-1):(m+1),(n-1):n);
elseif (m<m_max)&&(n>1)&&(n<n_max) && (m<1)
    original_mat = original(m:(m+1),(n-1):(n+1));
    m_c = 1;
elseif (m>1)&&(m<m_max)&&(n<n_max) && (n<1)
    original_mat = original((m-1):(m+1),n:(n+1));
    n_c = 1;
    
elseif (m>1)&&(n<n_max) && (m>m_max)&&(n<1)
    original_mat = original((m-1):m,n:(n+1));
    n_c = 1;
elseif (m<m_max)&&(n>1) && (n>n_max)&&(m<1)
    original_mat = original(m:(m+1),(n-1):n);
    m_c = 1;
elseif (m>1)&&(n>1) && (m>m_max)&&(n>n_max)
    original_mat = original((m-1):m,(n-1):n);
elseif (m<m_max)&&(n<n_max) && (m<1)&&(n<1)
    original_mat = original(m:(m+1),n:(n+1));
    m_c = 1;
    n_c = 1;
    
else
    error('Original image is smaller than 2x2');
end

c = [m_c,n_c];
end

%{
%% Nearest neighbours values difference check (version 1):
function log_val = nnvalcheck(original,reflectivity,filter,m,n,dif_mat,check_mat,change,sat)
% Center value:
ref = reflectivity;
ref(m,n) = ref(m,n) + check_mat(m,n)*change;
val = shortconv(original(m,n),ref,filter,m,n,sat,sat);

% "log_val = 1" Means it's the highest (">=") change in its surrounding, "0" otherwise.
log_val = 1;

% Compare the change to its surroundings:
if (m > 1)&&(check_mat(m-1,n) ~= 0)
    ref = reflectivity;
    ref(m-1,n) = ref(m-1,n) + check_mat(m-1,n)*change;
    val_ref = shortconv(original(m-1,n),ref,filter,m-1,n,sat);
    if (dif_mat(m,n) - val < dif_mat(m-1,n) - val_ref)
        log_val = 0;
    end
end
if (n > 1)&&(check_mat(m,n-1) ~= 0)
    ref = reflectivity;
    ref(m,n-1) = ref(m,n-1) + check_mat(m,n-1)*change;
    val_ref = shortconv(original(m,n-1),ref,filter,m,n-1,sat);
    if (dif_mat(m,n) - val < dif_mat(m,n-1) - val_ref)
        log_val = 0;
    end
end
if (m < size(reflectivity,1))&&(check_mat(m+1,n) ~= 0)
    ref = reflectivity;
    ref(m+1,n) = ref(m+1,n) + check_mat(m+1,n)*change;
    val_ref = shortconv(original(m+1,n),ref,filter,m+1,n,sat);
    if (dif_mat(m,n) - val < dif_mat(m+1,n) - val_ref)
        log_val = 0;
    end
end
if (n < size(reflectivity,2))&&(check_mat(m,n+1) ~= 0)
    ref = reflectivity;
    ref(m,n+1) = ref(m,n+1) + check_mat(m,n+1)*change;
    val_ref = shortconv(original(m,n+1),ref,filter,m,n+1,sat);
    if (dif_mat(m,n) - val < dif_mat(m,n+1) - val_ref)
        log_val = 0;
    end
end

end
%}

