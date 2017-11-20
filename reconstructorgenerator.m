function [rmatrix,success]= reconstructorgenerator(npix, Npix, pyramidmask, lambda, error, sampling, tripyramid, MVM)

%% Set up

% Centers of the pyramid pupils
% Sampling= # of pixels across the pyramid pupil 
if sampling==32
%     cen1=[38,38];
%     cen2=[92,38];
%     cen3=[38,92];
%     cen4=[92,92];

      cen1=[161,255];
      cen2=[255,161];
      cen3=[161,66];
      cen4=[66,161];
end

if sampling==64
%     cen1=[75,75];
%     cen2=[183,75];
%     cen3=[75,183];
%     cen4=[183,183];

      cen1=[321,510];
      cen2=[510,321];
      cen3=[321,132];
      cen4=[132,321];
end

if sampling==128
%     cen1=[150,150];
%     cen2=[366,150];
%     cen3=[150,366];
%     cen4=[366,366];
      cen1=[641,1019];
      cen2=[1019,641];
      cen3=[641,264];
      cen4=[264,641];
end


%% Zernike Generation
ncount=[];
mcount=[];
counter=1;
G_mat=[];
success=0;
%Scanning through Z -35 to Z 35
for n=0:5
    for m=-5:5
ncount(counter)=n;
mcount(counter)=m;
counter=counter+1;
ma = abs(m);
    if n==0 & m == 0
        continue
    elseif mod(n-ma,2)~=0
        continue
    elseif n<ma
        continue
    else  
% Generate the pupil with WF error        
    ef= zernike(0,0,npix).*exp(1i*((2*pi)/lambda)*error*zernike(n,m, npix));        
    pupil = complex(zeros(Npix));
    pupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1) =ef;
    success=success+1;
    
%% Fourier Transform to focal plane
FTpupil=fftshift(fft2(fftshift(pupil)))/(length(Npix).*length(Npix));
%figure, imagesc(abs(FTpupil).^2), title('PSF')

%% PWFS 
pyramid=FTpupil.*pyramidmask; %Apply the OPD mask to simulate the pyramid tip
% figure;imagesc(angle(pyramid));axis equal; title('Pyramid focal plane
% splitting')
%figure; imagesc(angle(pyramid));axis equal

%% Back to Pupil Plane and simulate detection
Pupilpyramid=abs(ifftshift(fft2(ifftshift(pyramid)))).^2;
figure; imagesc(Pupilpyramid); axis equal
%% Bin the pupil oixels 
[a,b]=size(Pupilpyramid);
p=npix/sampling;

Pupilpyramid=sum(reshape(Pupilpyramid, p, []),1);
Pupilpyramid=reshape(Pupilpyramid, a/p, []).';
Pupilpyramid=sum(reshape(Pupilpyramid, p, []),1);
Pupilpyramid=reshape(Pupilpyramid, b/p, []).';
%figure;imagesc(Pupilpyramid); axis equal

%% Intensity Centroid Calculation MVM
if MVM==true
PupilOne=Pupilpyramid(cen1(1)-sampling/2-1:cen1(1)+sampling/2+1, cen1(2)-sampling/2-1:cen1(2)+sampling/2+1);
PupilTwo=Pupilpyramid(cen2(1)-sampling/2-1:cen2(1)+sampling/2+1, cen2(2)-sampling/2-1:cen2(2)+sampling/2+1);
PupilThree=Pupilpyramid(cen3(1)-sampling/2-1:cen3(1)+sampling/2+1, cen3(2)-sampling/2-1:cen3(2)+sampling/2+1);
PupilFour=Pupilpyramid(cen4(1)-sampling/2-1:cen4(1)+sampling/2+1, cen4(2)-sampling/2-1:cen4(2)+sampling/2+1);

if tripyramid == true
    PupilFour=PupilFour.*0;
end

% figure
% subplot(2,2,1)
% imagesc(PupilOne);
% axis equal
% subplot(2,2,2)
% imagesc(PupilTwo);
% axis equal
% subplot(2,2,3)
% imagesc(PupilThree);
% axis equal
% subplot(2,2,4);
% imagesc(PupilFour);
% axis equal

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
end

%% Full Frame Method
if MVM==false
    Centroids=Pupilpyramid;
end

%% G-Matrix Generation
Gcolumn=reshape(Centroids',[size(Centroids,1)*size(Centroids,2) 1]);
G_mat(:,success)=Gcolumn;
G_mat(isnan(G_mat))=0;   

 
end
end
end

%% Build Reconstructor Matrix and save

[U Sig V] = svd(G_mat);
SigInT=1./Sig';
SigInT(isinf(SigInT))=0;
rmatrix=V*SigInT*U';

if tripyramid==true
    if MVM==true
        save('MVMtrireconstructormatrix.mat', 'rmatrix')
    end
    if MVM==false
        save('trireconstructormatrix.mat', 'rmatrix')
    end
end

if tripyramid==false
    if MVM==true
        save('MVMquadreconstructormatrix.mat', 'rmatrix')
    end
    if MVM==false
    save('quadreconstructormatrix.mat', 'rmatrix')
    end
end
   
    
end