function [rpupil]=zernpupil(error,success,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pyramidmask)

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


%call pyramidsim
rpupil(:,:,success)= pyramidsim(error,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pupil, pyramidmask);

    end
    end
end
end
    




%%end