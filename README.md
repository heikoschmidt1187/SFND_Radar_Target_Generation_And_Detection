# Project Writeup
The following project layout has been implemented:

<img src="imgs/projectLayout.png" width="700" />

## FMCW Waveform Design
For the given system requirements, an FMCW waveform has to be designed. These
requirements are the following:

<img src="imgs/systemRequirements.png" width="700" />

The waveform is characterized by Bandwith (Bsweep), Sweep Time (Tchirp) and
the slope of the chirp. This is done in the MATLAB code the following way:

```MATLAB
%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B = c / (2 * dres);                 % bandwith Bsweep
Tchirp = 5.5 * 2 * (Rmax / c);      % chirp time
slope = B / Tchirp;                 % slope of chirp

display(slope);
```

The resulting chirps look in general like the following example:

<img src="imgs/chirp.png" width="700" />

The calculated slope value is:

```console
slope =
    2.0455e+13
```

## Simulation Loop
In the simulation loop, a moving target has to be simulated and based on the
movement the beat or mixed signal for every timestamp has to be calculated:

```MATLAB
%Operating carrier frequency of Radar
fc= 77e9;             %carrier freq


%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)         


    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
    r_t(i) = R + v * t(i);      % based on s = v * t, mapped to range
    td(i) = 2 * r_t(i) / c;     % time delay

    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.

    % update transmit sample
    Tx(i) = cos(2 * pi * (fc * t(i) + slope * t(i)^2 / 2));                      

    % update receive sample - shifted by time delay tau
    Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + (slope * (t(i) - td(i))^2) / 2));

    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i) * Rx(i);     % mixed signal

end
```

## Range FFT (1st FFT)
A 1D-FFT on the range beat signal needs to be done and plotted:

```MATLAB
%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
signal_fft = reshape(Mix, [Nr, Nd]);    % reshape the linspace to 2D matrix Nr rows, Nd cols

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(signal_fft, Nr) / Nr;  % 1D FFT in Nr direction (col), normalized

 % *%TODO* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft);   % we are interested in the absolute value

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:(Nr / 2));  

%plotting the range
figure ('Name','Range from First FFT')
subplot(2, 1, 1)

 % *%TODO* :
 % plot FFT output
plot(signal_fft);
axis ([0 200 0 1]);
xlabel('range measured');
```

The resulting graph is depicted below, the initial position of the object is
chosen to be at 100 meters.

<img src="imgs/1DFFT.jpg" width="700" />

## 2D CFAR
For the 2D CFAR Implementation, first a doppler FFT is done to get the range
doppler matrix:

```MATLAB
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
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
```

The plot looks like this:

<img src="imgs/2DFFT.jpg" width="700" />


After that, the 2D CFAR is done an printed as graph:

```MATLAB
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;    % training cells for range
Td = 4;     % training cells for doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 8;     % guard cells for range
Gd = 4;     % guard cells for doppler

no_of_tcells = (2 * Tr + 2 * Gr + 1) * (2 * Td + 2 * Gd + 1) - ...
    (2 * Gr + 1) * (2 * Gd + 1);

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 10;   % threshold offset

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1, 1);      % noise level in range and doppler direction


% *%TODO* :
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

% slide over range dimension - /2 because only one half of the signal is
% used, keep spacing at the edges
sig_CFAR = zeros(size(RDM));
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;    % training cells for range
Td = 4;     % training cells for doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 8;     % guard cells for range
Gd = 4;     % guard cells for doppler

no_of_tcells = (2 * Tr + 2 * Gr + 1) * (2 * Td + 2 * Gd + 1) - ...
    (2 * Gr + 1) * (2 * Gd + 1);

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 10;   % threshold offset

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1, 1);      % noise level in range and doppler direction


% *%TODO* :
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

% slide over range dimension - /2 because only one half of the signal is
% used, keep spacing at the edges
sig_CFAR = zeros(size(RDM));

for r = (Tr + Gr + 1):(Nr / 2 - (Gr + Tr))

% slide over doppler dimension - full signal range is used, keep
% spacing at the edges
for d = (Td + Gd + 1):(Nd - (Gd + Td))

        % TODO: sum the signal level within the training cells, convert db2pow

        % get the patch with all training and guard cell + CUT (converted)
        train_patch = db2pow(RDM(r - (Tr + Gr) : r + (Tr + Gr), d - (Td + Gd) : d + (Td + Gd)));

        % zero CUT and guard cells as we are not interested in these now
        train_patch(Tr + 1 : end - Tr, Td + 1 : end - Td) = 0;

        (sum(sum(train_patch)) / no_of_tcells);

        % average summed values for used training cells and convert back
        % pow2db
        noise_level = pow2db(sum(sum(train_patch)) / no_of_tcells);

        % add the offset to determine the threshold
        threshold = noise_level + offset;

        % compare signal under CUT with threshold --> 0 or one
        if RDM(r, d) < threshold
            sig_CFAR(r, d) = 0;
        else
            sig_CFAR(r, d) = 1;
        end
    end
end


% *%TODO* :
% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.

% the not handled values are a boarder of training + guard cell
sig_CFAR(1 : Nr / 2, 1 : Td + Gd) = 0;          % front cols
sig_CFAR(1 : Nr / 2, Nd - (Td + Gd) : Nd) = 0;  % back cols
sig_CFAR(1 : Tr + Gr, 1 : Nd) = 0;                  % top rows
sig_CFAR(Nr / 2 - (Tr + Gr) : Nr / 2, 1 : Nd) = 0;  % bottom rows


% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered RDM'),surf(doppler_axis, range_axis, sig_CFAR);
colorbar;
for r = (Tr + Gr + 1):(Nr / 2 - (Gr + Tr))

    % slide over doppler dimension - full signal range is used, keep
    % spacing at the edges
    for d = (Td + Gd + 1):(Nd - (Gd + Td))

        % TODO: sum the signal level within the training cells, convert db2pow

        % get the patch with all training and guard cell + CUT (converted)
        train_patch = db2pow(RDM(r - (Tr + Gr) : r + (Tr + Gr), d - (Td + Gd) : d + (Td + Gd)));

        % zero CUT and guard cells as we are not interested in these now
        train_patch(Tr + 1 : end - Tr, Td + 1 : end - Td) = 0;

        (sum(sum(train_patch)) / no_of_tcells);

        % average summed values for used training cells and convert back
        % pow2db
        noise_level = pow2db(sum(sum(train_patch)) / no_of_tcells);

        % add the offset to determine the threshold
        threshold = noise_level + offset;

        % compare signal under CUT with threshold --> 0 or one
        if RDM(r, d) < threshold
            sig_CFAR(r, d) = 0;
        else
            sig_CFAR(r, d) = 1;
        end
    end
end


% *%TODO* :
% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.

% the not handled values are a boarder of training + guard cell
sig_CFAR(1 : Nr / 2, 1 : Td + Gd) = 0;          % front cols
sig_CFAR(1 : Nr / 2, Nd - (Td + Gd) : Nd) = 0;  % back cols
sig_CFAR(1 : Tr + Gr, 1 : Nd) = 0;                  % top rows
sig_CFAR(Nr / 2 - (Tr + Gr) : Nr / 2, 1 : Nd) = 0;  % bottom rows


% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered RDM'),surf(doppler_axis, range_axis, sig_CFAR);
colorbar;
```

<img src="imgs/CACFAR.jpg" width="700" />

### Implementation Steps for 2D CFAR

For the CFAR, first the number of training and guard cells around the cell
under test (CUT) in both, range and doppler direction, have been chosen with a
trial and error approach. The same has been done for the offset value:

```MATLAB
% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;    % training cells for range
Td = 4;     % training cells for doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 8;     % guard cells for range
Gd = 4;     % guard cells for doppler

no_of_tcells = (2 * Tr + 2 * Gr + 1) * (2 * Td + 2 * Gd + 1) - ...
    (2 * Gr + 1) * (2 * Gd + 1);

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 10;   % threshold offset
```

The 2D CFAR itself has been implemented in a loop, where in the indices the
number of training and guard cells have been taken into account to not "shoot"
over the edges on looping:

```MATLAB
% slide over range dimension - /2 because only one half of the signal is
% used, keep spacing at the edges
sig_CFAR = zeros(size(RDM));

for r = (Tr + Gr + 1):(Nr / 2 - (Gr + Tr))

    % slide over doppler dimension - full signal range is used, keep
    % spacing at the edges
    for d = (Td + Gd + 1):(Nd - (Gd + Td))
        ...
    end
end
```

Inside the loop, first the patch around the CUT is extracted based on the
training and guard cells. As the values of the trainings cells have to be
summed up, the guard cells and the CUT are set to 0 for now to not influence the
sum and therefor the average of all training cells:

```MATLAB
% get the patch with all training and guard cell + CUT (converted)
train_patch = db2pow(RDM(r - (Tr + Gr) : r + (Tr + Gr), d - (Td + Gd) : d + (Td + Gd)));

% zero CUT and guard cells as we are not interested in these now
train_patch(Tr + 1 : end - Tr, Td + 1 : end - Td) = 0;

% average summed values for used training cells and convert back
% pow2db
noise_level = pow2db(sum(sum(train_patch)) / no_of_tcells);
```

The resulting noise level is then offset by the previously defined value. All
range doppler values above the threshold lead to a CFAR signal of 1, all others
to 0:

```MATLAB
% add the offset to determine the threshold
threshold = noise_level + offset;

% compare signal under CUT with threshold --> 0 or one
if RDM(r, d) < threshold
    sig_CFAR(r, d) = 0;
else
    sig_CFAR(r, d) = 1;
end
```

Because the loops are not covering all values on the given matrix, there may be
false positives at the edges. To suppress these, all the cells that have never
been considered as CUT are set to zero:

```MATLAB
% the not handled values are a boarder of training + guard cell
sig_CFAR(1 : Nr / 2, 1 : Td + Gd) = 0;          % front cols
sig_CFAR(1 : Nr / 2, Nd - (Td + Gd) : Nd) = 0;  % back cols
sig_CFAR(1 : Tr + Gr, 1 : Nd) = 0;                  % top rows
sig_CFAR(Nr / 2 - (Tr + Gr) : Nr / 2, 1 : Nd) = 0;  % bottom rows
```
