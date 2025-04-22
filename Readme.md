
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

### FMCW Waveform Design

#### Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.


### Simulation Loop

#### Simulate Target movement and calculate the beat or mixed signal for every timestamp
For given system requirements the calculated slope should be around 2e13

### Range FFT (1st FFT)

#### Implement the Range FFT on the Beat or Mixed Signal and plot the result
A correct implementation should generate a peak at the correct range, i.e the initial position of target assigned with an error margin of +/- 10 meters.

#### Implement the Range FFT on the Beat or Mixed Signal and plot the result
A beat signal should be generated such that once range FFT implemented, it gives the correct range i.e the initial position of target assigned with an error margin of +/- 10 meters.

### 2D CFAR

#### Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map
The 2D CFAR processing should be able to suppress the noise and separate the target signal. The output should match the image shared in walkthrough.

#### Create a CFAR README File
	
In a README file, write brief explanations for the following:

- Implementation steps for the 2D CFAR process.
- Selection of Training, Guard cells and offset.
- Steps taken to suppress the non-thresholded cells at the edges.

