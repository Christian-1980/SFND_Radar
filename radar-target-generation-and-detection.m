clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Speed of light
c = 3*10^8; % [m/s]

%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains contant

% A struct to have a container to hold all relevant data
radar = struct();
radar.frequency_operational = 77*10^9; % [GHz]
radar.range_max = 200; % [m]
radar.res_max = 1; % [m]

% A struct to have a container to hold all relevant data
target = struct();
target.v_max = 70; % [m/s]
target.v_resolution = 3; % [m/s]
target.init_range = 110;
target.init_v = -20;

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of iradar.Ts parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requiremenradar.Ts above.

% radar bandwidth
radar.bsweep = c / (2 * radar.res_max);

% fitting factor
radar.fitting_factor = 5.5; % practical factor for fitting

% radar.Ts
radar.Ts = radar.fitting_factor * 2 * radar.range_max / c; 

% Slope
radar.chirp_slope = radar.bsweep / radar.Ts;


%% Print ou the radar props
fprintf('-----------------\n');
fprintf('Radar Properties:\n');
fprintf('  Operational Frequency: %.2e Hz\n', radar.frequency_operational);
fprintf('  Maximum Range: %.2f m\n', radar.range_max);
fprintf('  Range Resolution: %.2f m\n', radar.res_max);
fprintf('  Bandwidth: %.2f m\n', radar.bsweep);
fprintf('  Fitting Factor: %.2f m\n', radar.fitting_factor);
fprintf('  Swing time: %.2f m\n', radar.Ts);
fprintf('-----------------\n');

                                                          
%The number of chirps in one sequence. Iradar.Ts ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*radar.Ts,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
t_delta=zeros(1,length(t));


%% Signal generation and Moving Target simulation

% Running the radar scenario over the time. 
for i=1:length(t)         

    %For each time stamp update the Range of the Target for constant velocity.
    range = target.init_range + target.init_v*t(i)

    %For each time sample we need update the transmitted and
    %received signal. 
    t_delta= 2*range / c;
    t_new  = t(i)-t_delta;

    Tx(i) = cos( 2*pi*( radar.frequency_operational*t(i) + (0.5*radar.chirp_slope*t(i)^2) ) );
    Rx(i) = cos( 2*pi*( radar.frequency_operational*t_new + (0.5*radar.chirp_slope*t_new^2) ) );

    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end


%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
new_Mix = reshape(Mix,[Nr,Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(new_Mix);

% Take the absolute value of FFT output
signal_fft = abs(signal_fft/max(max(signal_fft)));

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal = signal_fft(1:length(t)/2+1);   % Taking only half of output

% Also plotting the range
figure ('Name','Range based on First FFT')
f = length(t)*(0:length(t)/2)/length(t);
plot(f,signal) 
axis ([0 200 0 1]);


%% RANGE DOPPLER RESPONSE

% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift(sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);

% Plot the Range Doppler Map
figure, surf(doppler_axis, range_axis, RDM);
title('Range Doppler Map');
xlabel('Doppler Velocity (m/s)');
ylabel('Range (m)');
zlabel('Amplitude (dB)');
colorbar;


%% CFAR implementation
%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
training_cells_chirps = 10;
training_cells_samples = 10;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
guard_cells_chrips = 4;
guard_cells_samples = 4;

% offset the threshold by SNR value in dB
offset_threshold = 4.5;

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR

% Getting the dimensions
[m,n] = size(RDM);
% Vector to hold threshold values 
threshold_cfar = zeros(m,n);
%Vector to hold final signal after thresholding
signal_cfar = zeros(m,n);

% helpers
chirps_width = training_cells_chirps+guard_cells_chrips;
samples_width = training_cells_samples+guard_cells_samples;

% number of cells for averaging
num_cells = (2*chirps_width + 1)*(2*samples_width + 1) - (2*guard_cells_chrips + 1)*(2*guard_cells_samples + 1);

% loop for CUT
for i =  (chirps_width + 1):( m - 2*chirps_width)
    for j = (samples_width + 1):(n - 2*(samples_width))
        % loop to calculate cfar
        threshold_cfar(i,j) = sum(sum(db2pow(RDM(i-(chirps_width) : i+(chirps_width),j-(samples_width) : j+(samples_width))))); 
        threshold_cfar(i,j) = threshold_cfar(i,j) - sum(sum(db2pow(RDM((i-guard_cells_chrips):(i+guard_cells_chrips),(j-guard_cells_samples):(j+guard_cells_samples)))));
        
        threshold_cfar(i,j) = threshold_cfar(i,j)/num_cells;
        % calculate the threshold
        threshold_cfar(i,j) = offset_threshold + pow2db(threshold_cfar(i,j));
        
        if RDM(i,j) > threshold_cfar(i,j)
            signal_cfar(i,j) = 1;
        else
            signal_cfar(i,j) = 0;
        end
        
    end
end

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar);
title('CFAR Output');
xlabel('Doppler Velocity (m/s)');
ylabel('Range (m)');
zlabel('Detection');
colorbar;