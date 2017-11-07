function rpupil=pyramidsim(error, npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pupil, pyramidmask)

%% Set up

% Centers of the pyramid pupils
% Sampling= # of pixels across the pyramid pupil 
if sampling==32
    cen1=[38,38];
    cen2=[92,38];
    cen3=[38,92];
    cen4=[92,92];
end

if sampling==64
    cen1=[75,75];
    cen2=[183,75];
    cen3=[75,183];
    cen4=[183,183];
end

if sampling==128
    cen1=[150,150];
    cen2=[366,150];
    cen3=[150,366];
    cen4=[366,366];
end
    

%% Fourier Transform to focal plane
FTpupil=fftshift(fft2(fftshift(pupil)))/(length(Npix).*length(Npix));
%figure, imagesc(abs(FTpupil).^2), title('PSF')

%% PWFS 
pyramid=FTpupil.*pyramidmask; %Apply the OPD mask to simulate the pyramid tip
% figure;imagesc(angle(pyramid));axis equal; title('Pyramid focal plane
% splitting')

%% Back to Pupil Plane and simulate detection
Pupilpyramid=abs(ifftshift(fft2(ifftshift(pyramid)))).^2;

%% Bin the pupil oixels 
[a,b]=size(Pupilpyramid);
p=npix/sampling;

Pupilpyramid=sum(reshape(Pupilpyramid, p, []),1);
Pupilpyramid=reshape(Pupilpyramid, a/p, []).';
Pupilpyramid=sum(reshape(Pupilpyramid, p, []),1);
Pupilpyramid=reshape(Pupilpyramid, b/p, []).';

%% Intensity Centroid Calculation MVM
if MVM==true
PupilOne=Pupilpyramid(cen1(1)-sampling/2-1:cen1(1)+sampling/2+1, cen1(2)-sampling/2-1:cen1(2)+sampling/2+1);
PupilTwo=Pupilpyramid(cen2(1)-sampling/2-1:cen2(1)+sampling/2+1, cen2(2)-sampling/2-1:cen2(2)+sampling/2+1);
PupilThree=Pupilpyramid(cen3(1)-sampling/2-1:cen3(1)+sampling/2+1, cen3(2)-sampling/2-1:cen3(2)+sampling/2+1);
PupilFour=Pupilpyramid(cen4(1)-sampling/2-1:cen4(1)+sampling/2+1, cen4(2)-sampling/2-1:cen4(2)+sampling/2+1);
end

if tripyramid == true
    PupilFour=PupilFour.*0;
end

%% Calculate Centroids 

sz=size(PupilOne);
Sx=zeros(sz(1),sz(2));
Sy=Sx;
if MVM==true

% 4PWFS quad-cell equation
if tripyramid==false
    Sx= (PupilOne+PupilTwo-PupilThree-PupilFour)/(PupilOne+PupilTwo+PupilThree+PupilFour);
    Sy= (PupilOne-PupilTwo-PupilThree+PupilFour)/(PupilOne+PupilTwo+PupilThree+PupilFour);
end

% 3PWFS tri-cell equation
if tripyramid==true
    Sx= (PupilTwo*(sqrt(3)/2)-PupilThree*(sqrt(3)/2))/(PupilOne+PupilTwo+PupilThree);
    Sy= (PupilOne-PupilTwo*(0.5)-PupilThree*(0.5))/(PupilOne+PupilTwo+PupilThree); 
end 
 

Centroids=[];
for k=1:sz(2)
        Centroids=[Centroids, Sx(:,k)];
        Centroids=[Centroids, Sy(:,k)];    
end
end

%% Full Frame Method
if MVM==false
    Centroids=Pupilpyramid;
end

Gcolumn=reshape(Centroids',[size(Centroids,1)*size(Centroids,2) 1]);
%apply reconstructor matrix to find Rwavefront. Rwavefront is a single
%vector containing the magnitude of each zernike detected.

Rwavefront= rmatrix*Gcolumn;
Rwavefront(abs(Rwavefront)<1e-6)=0;
count=1;
pupilsize=size(pupil);
el=zeros(npix,npix);

ncount=[];
mcount=[];
counter=0;
success=0;
for n=0:5
    for m=-5:5
         
        counter=counter+1;
        ma = abs(m);
    if n==0 & m == 0
        continue
    elseif mod(n-ma,2)~=0
        continue
    elseif n<ma
        continue
    else
        success=success+1;
        ncount(success)=n;
        mcount(success)=m;
        %figure; imagesc(zernike(n,m,npix).*Rwavefront(success,1));
        el=el+zernike(n,m,npix).*Rwavefront(success,1);
%         figure; imagesc(angle(rpupil)); title('Complex part of Reconstructed wavefront')
%         figure; imagesc(real(rpupil)); title('Real part of Reconstructed wavefront')
    end
    end
end
WF=zernike(0,0,npix).*exp(1i*(2*pi./lambda)*el*error);
rpupil = complex(zeros(Npix));
rpupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1)= WF;




end