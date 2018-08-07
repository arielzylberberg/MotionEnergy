function [motion_energy,eright,eleft] = motionfilter(stim,sampfreq_hz,varargin)
% stim: time, irrelevant direction, relevant direction
% http://www.lifesci.sussex.ac.uk/home/George_Mather/EnergyModel.htm

sigmac     = 0.5; %was 0.35 in kiani&shadlen 2008
sigmag     = 0.05; %deg
t_const    = 60;
degperbin  = 0.0260; 
doplot     = 0;
speed_deg_per_sec_for_plot = 5; %only used for plotting

for i=1:length(varargin)
    if isequal(varargin{i},'sigmac')
        sigmac = varargin{i+1};
    elseif isequal(varargin{i},'sigmag')
        sigmag = varargin{i+1};
    elseif isequal(varargin{i},'t_const')
        t_const = varargin{i+1};
    elseif isequal(varargin{i},'degperbin')    
        degperbin = varargin{i+1};
    elseif isequal(varargin{i},'speed_deg_per_sec')    
        speed_deg_per_sec_for_plot = varargin{i+1}; % only for plotting
    elseif isequal(varargin{i},'doplot')
        doplot = varargin{i+1};
    end
end

dt = 1/sampfreq_hz; % time per sample
dx = degperbin;

%--------------------------------------------------------------------------
%           Create component spatiotemporal filters
%--------------------------------------------------------------------------

max_x       = 0.8;         %Half-width of filter (deg)
x_filt_aux  = 0:dx:max_x;
x_filt      = sort(unique([x_filt_aux -x_filt_aux]));% A row vector holding spatial sampling intervals
nx          = length(x_filt);

alfa        = atan(x_filt/sigmac);
even_x      = (cos(alfa).^4 .* cos(4*alfa)); %even
even_x      = reshape(even_x,[1 1 nx]);
odd_x       = (cos(alfa).^4 .* sin(4*alfa)); %odd
odd_x       = reshape(odd_x,[1 1 nx]);

% spatial filter in the orthogonal dimension
max_y       = 0.2;
y_filt_aux  = 0:dx:max_y;
% A row vector holding spatial sampling intervals
y_filt      = sort(unique([y_filt_aux -y_filt_aux]));
ny          = length(y_filt);
gauss_y     = exp(-y_filt.^2/(2*sigmag.^2));% despues permuto dimensiones
gauss_y     = reshape(gauss_y,[1 ny 1]);


%tiempo
max_t       = 0.3;      % Duration of impulse response (sec)
t_filt      = [0:dt:max_t]'; % A column vector holding temporal sampling intervals
nt          = length(t_filt);

slow_t      = (t_const*t_filt).^5.*exp(-t_const*t_filt).*(1/factorial(5)-(t_const*t_filt).^2/factorial(5+2));
fast_t      = (t_const*t_filt).^3.*exp(-t_const*t_filt).*(1/factorial(3)-(t_const*t_filt).^2/factorial(3+2));

%--------------------------------------------------------------------------
%           Convolve
%--------------------------------------------------------------------------

% aux calculations
stimy               = convn(gauss_y,stim);
oddx_stimy          = convn(odd_x,stimy);
evenx_stimy         = convn(even_x,stimy);
fast_oddx_stimy     = convn(fast_t,oddx_stimy);
slow_oddx_stimy     = convn(slow_t,oddx_stimy);
fast_evenx_stimy    = convn(fast_t,evenx_stimy);
slow_evenx_stimy    = convn(slow_t,evenx_stimy);

% go
resp_right_1    = - fast_oddx_stimy + slow_evenx_stimy;
resp_right_2    =   slow_oddx_stimy + fast_evenx_stimy;
resp_left_1     =   fast_oddx_stimy + slow_evenx_stimy;
resp_left_2     = - slow_oddx_stimy + fast_evenx_stimy;

%--------------------------------------------------------------------------
%         Square the filter output
%--------------------------------------------------------------------------
resp_left_1     = resp_left_1.^2;
resp_left_2     = resp_left_2.^2;
resp_right_1    = resp_right_1.^2;
resp_right_2    = resp_right_2.^2;

%--------------------------------------------------------------------------
%         Add
%--------------------------------------------------------------------------

energy_right    = resp_right_1 + resp_right_2;
energy_left     = resp_left_1 + resp_left_2;

eright = sum(sum(energy_right,2),3);
eleft  = sum(sum(energy_left,2),3);

%--------------------------------------------------------------------------
%         Calculate net energy as the R-L difference
%--------------------------------------------------------------------------
motion_energy = eright - eleft;


%% plot the filters
% doplot = 1;
if (doplot)
    
    % % Step 1c: combine space and time to create spatiotemporal filters
    e_slow = slow_t * even_x(:)';    %SE/TS
    e_fast = fast_t * even_x(:)';   %SE/TF
    o_slow = slow_t * odd_x(:)';   %SO/TS
    o_fast = fast_t * odd_x(:)';   % SO/TF
    % %
    % % %--------------------------------------------------------------------------
    % % %         STEP 2: Create spatiotemporally oriented filters
    % % %--------------------------------------------------------------------------
    left_1  = o_fast + e_slow;      % L1
    left_2  = -o_slow + e_fast;     % L2
    right_1 = -o_fast + e_slow;    % R1
    right_2 = o_slow + e_fast;     % R2

    
    figure(3)
%     xx = [0 10]*degperbin;
    
    xx = [0 1]*speed_deg_per_sec_for_plot*(3*dt);%para plotear 5 grados por sec
    tt = [0 3]*dt;
    imagesc(x_filt,t_filt,left_1)
    hold on; plot(xx,tt,-xx,tt,'b')
    xlabel('deg'); ylabel('time')
    
    figure(4)
    subplot(2,2,1);
    imagesc(x_filt,t_filt,left_1)
    hold on; plot(xx,tt,-xx,tt,'b')
    xlabel('deg'); ylabel('time')
    subplot(2,2,2); imagesc(x_filt,t_filt,left_2)
    hold on; plot(xx,tt,-xx,tt,'b')
    subplot(2,2,3); imagesc(x_filt,t_filt,right_1)
    hold on; plot(xx,tt,-xx,tt,'b')
    subplot(2,2,4); imagesc(x_filt,t_filt,right_2)
    hold on; plot(xx,tt,-xx,tt,'b')
    
    %fourier - in prep
    NFFTt = 100;
    NFFTx = 100;
    ff = abs(fftshift(fft2(right_1,NFFTt,NFFTx)));
    Fsx = 1/dx;
    fx  = Fsx/2*[-1:2/NFFTx:1-2/NFFTx];
    Fst = 1/dt;
    ft  = Fst/2*[-1:2/NFFTt:1-2/NFFTt];
    
    figure(5)
    imagesc(ft,fx,ff')
    xlabel('temporal freq (Hz)')
    ylabel('spatial freq (c/deg)')
    hold on
    xx = xlim;
    speed = speed_deg_per_sec_for_plot;%deg/sec
    yy = -1/speed * xx; 
    plot(xx,yy,'r')
    
%     pause

    
end
