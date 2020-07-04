# SFND_Radar_Target_Generation_and_Detection
Udacity SFND Radar lesson project 

## Reflection
### FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (T_chirp) and slope of the chirp.
1. Radar specifications
```
Frequency of operation = 77 GHz
Max Range = 200 m
Range Resolution = 1 m
Max Velocity = 100 m/s
speed of light = 3e8 m/s

target init range = 100 m
target init velocity = 15 m
```
2. FMCW waveform design 
```
B=c/(2*range_resolution);
T_chirp=5.5*(2*max_range)/c;
slope=B/T_chirp;
slope = 2.0455e13
```

### Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.
1. Variable setup
```
% Operating carrier frequency of Radar 
fc= 77e9;                       % carrier freq

% The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT for Doppler Estimation. 
Nd=128;                         % #of doppler cells OR #of sent periods % number of chirps

% The number of samples on each chirp. 
Nr=1024;                        %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each chirp
t=linspace(0,Nd*T_chirp,Nr*Nd); %total time for samples

% Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t));          %transmitted signal
Rx=zeros(1,length(t));          %received signal
Mix = zeros(1,length(t));       %beat signal

% Similar vectors for range_covered and time delay.
range_covered=zeros(1,length(t));
time_delay=zeros(1,length(t));
```
2. Loop
```
for i=1:length(t)         
    % For each time stamp update the Range of the Target for constant velocity. 
    range_covered(i)=target_init_range+target_init_velocity*t(i);
    time_delay=2*range_covered/c;
    
    % For each time sample we need update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-time_delay(i)) + (slope*(t(i)-time_delay(i))^2)/2));
    
    % Now by mixing the Transmit and Receive generate the beat signal
    % This is done by element wise matrix multiplication of Transmit and Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end
```

### Range FFT (1st FFT)
Implement the Range FFT on the Beat or Mixed Signal.
```
% Run the FFT on the beat signal along the range bins dimension (Nr) and normalize.
mix_fft=fft(Mix,Nr)./Nr;

% Take the absolute value of FFT output
mix_fft=abs(mix_fft);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
mix_fft=mix_fft(1:Nr/2);
```
![1st_FFT](https://github.com/allenhyp/SFND_Radar_Target_Generation_and_Detection/blob/master/result_image/1D_FFT.png?raw=true)

### 2D CFAR
1. Output of 2D FFT operation, i.e the Range Doppler Map.
![2D_FFT](https://github.com/allenhyp/SFND_Radar_Target_Generation_and_Detection/blob/master/result_image/2D_FFT.png?raw=true)

2. CFAR setup
```
% Select the number of Training Cells in both the dimensions.
Tr=20;
Td=16;

% Select the number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation
Gr=10;
Gd=8;

% Offset the threshold by SNR value in dB
offset=8;

% Create a vector to store noise_level for each iteration on training cells
window_size=(2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
num_training_cells=window_size-(2*Gr+1)*(2*Gd+1);
signals=zeros(size(RDM));
```
3. Slide window though the Range Doppler Map
```
[RDM_r, RDM_d]=size(RDM);
for i=1:RDM_d-2*(Td+Gd)
    for j=1:RDM_r-2*(Tr+Gr)
        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR
        training_window=db2pow(RDM(j:j+2*(Tr+Gr),i:i+2*(Td+Gd)));
        training_window(Tr+1:Tr+2*Gr+2,Td+1:Td+2*Gd+2)=0;
        noise_level=pow2db(sum(training_window)/num_training_cells);
        threshold=noise_level*offset;
        
        % The process above will generate a thresholded block, which is smaller
        % than the Range Doppler Map as the CUT cannot be located at the edges of
        % matrix. Hence,few cells will not be thresholded. To keep the map size same
        % set those values to 0.
        target_d=i+Td+Gd;
        target_r=j+Tr+Gr;
        if RDM(target_r,target_d) > threshold
            signals(target_r,target_d)=1;
        end
    end
end
```
4. Ouput of 2D CA-CFAR result
![2D_CA-CFAR](https://github.com/allenhyp/SFND_Radar_Target_Generation_and_Detection/blob/master/result_image/2D_CA-CFAR.png?raw=true)
As the image shown above, the target with range ~= 100 and velocity ~= 15 was detected.

