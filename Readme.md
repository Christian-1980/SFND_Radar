
# Radar Target Generation and Detection

This project uses Matlab to introduce frequency modulated continuous-wave (FMCW) radar and related post-processing techniques. The topics covered include:
- Fast Fourier transforms (FFT) and 2 dimensional FFT
- Clutter v. Target
- Engineer chirp bandwith to meet system requirements for max resolution
- Phased array beam steering to determine angle of arrival (AoA)
- Noise suppression by applying Constant false alarm rate (CFAR)
- Signal-to-noise ratio (SNR) and dynamic thresholding

## Visualization of the results


## Installing Matlab
Instructions for installing the latest version of Matlab can be found at https://www.mathworks.com/

I used the onlin version that is provided by Udacity account.

## Project requirements and solution

In order to have a clean data storage a struct is used for both, the radar and the target.

### Radar Specicfication

| Specifications | Value |
:----:|:----:
|Frequency | 77 Ghz|
|Range Resolution | 1 m|
|Max Range | 200 m|
|Max velocity | 70 m/s|
|Velocity resolution| m m/s|

### Target Spoecification

| Specifications | Value |
:----:|:----:
| Velocity | -20 m/s | 
| Range | 110 m |

### FMCW Waveform Design

#### Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.
```
%% FMCW Waveform Generation

% Design the FMCW waveform by giving the specs of each of iradar.Ts parameters.
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
```

### Simulation Loop

#### Simulate Target movement and calculate the beat or mixed signal for every timestamp
For given system requirements the calculated slope should be around 2e13

### Range FFT (1st FFT)

#### Implement the Range FFT on the Beat or Mixed Signal and plot the result

A beat signal should be generated such that once range FFT implemented, it gives the correct range i.e the initial position of target assigned with an error margin of +/- 10 meters.

- Implement the 1D FFT on the Mixed Signal

```
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
```



- Reshape the vector into Nr*Nd array.
- Run the FFT on the beat signal along the range bins dimension (Nr)
- Normalize the FFT output with length, L = Bsweep * Tchirp.
- Take the absolute value of that output.
- Keep one half of the signal
- Plot the output
- There should be a peak at the initial position of the target. In our case at 110 m.

```
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
```

### 2D CFAR

#### Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map
The 2D CFAR processing should be able to suppress the noise and separate the target signal. The output should match the image shared in walkthrough.

#### Create a CFAR README File
	
In a README file, write brief explanations for the following:

- Implementation steps for the 2D CFAR process.
- Selection of Training, Guard cells and offset.
- Steps taken to suppress the non-thresholded cells at the edges.

