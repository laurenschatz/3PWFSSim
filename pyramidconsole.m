%Pyramid Wavefront Sensor Module
%Simulates a refractive PWFS

% Lauren Schatz
% October 19, 2017


npix= 512; %number of pixels across the pupil diameter
Npix= 2048; % number of pixels across the total roster. Represents the amount of zero padding.
lambda= 700.*10^-9; %wavelength in nanometers
error=1*10^-9; % nanometers of error you will be applying though the WFS
sampling=64; % 128; 32 %Numberof pixels across each pupil. chose 32, 64,or 128 le
rooferror=2; % in pixels. If you have a roof error, number of pixels in diameter creating the roof effect. 
tabletoperror=2; %in pixels. If you have a tabletop error, number of pixels in diameter creating the flattened tabletop of the pyramid tip. 

%% Set Up

%Simulator setup

tripyramid=false; %toggle true/false for 3PWFS simulation. If false, will run 4PWFS simulation
have_reconstructor=false; %toggle true/false for reconstructor matrix
have_pyramidmask=true;  %toggle true/fasle for pyramid focal plane maske.

roof=false; %toggle true/false for roof error
tabletop=false; % toggle true/false for tabletop error
    
MVM=true; %Toggle true/false. True: Will run Matrix-Vector-Multiply (intensity Centroids) for reconstructor matrix
                             %Flase: Will use the full frame (Olivier
                             %Guyon's method) for the reconstructor matrix

%WFS pupils: the type of phase error to be thrown into the PWFS sim.
    %zpupil: Feeds in zernike polynomial errors one by one. Does this
        %for Z -35:35;
    %phasescrn: Feeds in a single phase screen generated using zernike
        %polynomials in a Kolmogorov turbulence power spectrum.
    %modulatedzern: Same as zernpupil but with a 16 point modulation
        %pattern. (COMING SOON)
    %modulatedphasescrn: Same as phasescrn but with a 16 point modulation
        %pattern. (COMING SOON)

%% Run Simulation

%% Load in or generate the pyramid focal plane masks
if have_pyramidmask==false
    fprintf('generating pyramid mask')
    pyramidmask=maskgenerator(Npix, tripyramid, tabletoperror,tabletop, roof, rooferror);
    
end
    
if have_pyramidmask==true && tripyramid==true
    fprintf('loading pyramid mask')
    load 'tripyramidmask.mat';
   
end

if have_pyramidmask==true && tripyramid==false
    fprintf('loading pyramid mask')
    load 'quadpyramidmask.mat';
end

%% Errors on the pyramid

if tabletop==true;
    
    
end



% Load in or generate the reconstructor matrix
if have_reconstructor==false
    fprintf('Generating Reconstructor matrix')
    [rmatrix,success]=reconstructorgenerator(npix,Npix, pyramidmask, lambda, error, sampling, tripyramid, MVM);
end

if have_reconstructor==true
    if tripyramid==true && MVM ==true
        fprintf('3PWFS MVM')
        r=load('MVMtrireconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==true && MVM ==false
        fprintf('3PWFS full frame')
        r=load('trireconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==false && MVM ==true
        fprintf('4PWFS MVM')
        r=load('MVMquadreconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==false && MVM ==false
        fprintf('4PWFS full frame')
        r=load('quadreconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
end
        
%%   Run the WFS

%Returns variable WF, that is a data cube of the reconstructed wavefronts
%from PWFS.
rpupil=[];
%Calls on function pyramidsim
if zpupil==true
    [rpupil]=zernpupil(error,success,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pyramidmask);
end
%rpupil= pyramidsim(npix, Npix, sampling, rmatrix, lambda, sampling, tripyramid, MVM);

    

    

    
    
    