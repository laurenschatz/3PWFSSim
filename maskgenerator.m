%Generates pyramid focal plane mask.

%Saves the pyramid mask generated in a .mat file
function [pyramidmask]=maskgenerator(Npix, tripyramid, tabletoperror ,tabletop, roof, rooferror)
%%
% Npix=2048;
% npix=512;
% tripyramid=true;
% pyramid=tripyramid;

    ef= zernike(0,0,npix).*exp(1i*((2*pi)/lambda).*zernike(0,0, npix));        
    pupil = complex(zeros(Npix));
    pupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1) =ef;
    
FTpupil=fftshift(fft2(fftshift(pupil)))/(length(Npix).*length(Npix));
%%
%coordinate system setup.
x= (-Npix/2:Npix/2-1);
[X,Y]=meshgrid(x,x);
[phi,r]=cart2pol(X,Y); % phi from -pi to pi

%zernike polynomials (tip/tilt)
tilt=.295;
zernikez={r.*cos(phi),r.*sin(phi)};
% W1=tilt*(cos(-3*pi/4).*zernike{1}+sin(-3*pi/4).*zernike{2});
%W2= tilt*(cos(-pi/4).*zernikez{1}+sin(-pi/4).*zernikez{2});
%W3= -tilt*(cos(pi/4).*zernikez{1}-sin(pi/4).*zernikez{2});
%W4= -tilt*(cos(3*pi/4).*zernikez{1}-sin(3*pi/4).*zernikez{2});

W1=tilt*(zernikez{2});
W2=tilt*zernikez{1};
W3=tilt*-zernikez{2};
W4=tilt*-zernikez{1};
%Masks;
%%

if tripyramid==true

%Mask Generation
    %first 1/3 facet
    part1= phi>=pi/4-pi;
    part2= phi<=2*pi/3+pi/4-pi;
    tri1=part1.*part2;
        %second 1/3 facet
        part3= phi>=2*pi/3+pi/4-pi;
        part4= phi<=4*pi/3+pi/4-pi;
        tri2=part3.*part4;
            %third 1/3 facet
            tri3=tri1+tri2<1;
            %figure; imagesc(tri3+2*tri2+3*tri1);
%Add in tilts

focalt1= exp(1i*2*pi*W1).*tri1;%Complex wavefront
focalt2= exp(1i*2*pi*W2).*tri2;%Complex wavefront
focalt3= exp(1i*2*pi*W3).*tri3;%Complex wavefront

pyramidmask=focalt1+focalt2+focalt3;
figure; imagesc(angle(pyramidmask));
save('tripyramidmask.mat','pyramidmask')

if tabletop==true
flat=pyramidmask(Npix/2, Npix/2);
pyramidmask(Npix/2-tabletoperror/2:Npix/2+tabletoperror/2,Npix/2-tabletoperror/2:Npix/2+tabletoperror/2)=flat;
figure; imagesc(angle(pyramidmask));
figure; imagesc(angle(pyramidmask));
save('quadpyramidmask.mat','pyramidmask')
end

end
    
%%
if tripyramid==false
%Mask Generation
    %first 1/4 facet
     part1= phi>=pi/4-pi;
     part2= phi<= pi/2+pi/4-pi;
     quad1= part1.*part2;
            %second 1/4 facet
            part3= phi>= pi/2+pi/4-pi;
            part4= phi<= pi+pi/4-pi;
            quad2=part3.*part4;
                %third 1/4 facet
                part5= phi>=pi+pi/4-pi;
                part6= phi<= 3*pi/2+pi/4-pi;
                quad3=part5.*part6;
                    %fourth 1/4 facet
                    quad4= quad1+quad2+quad3<1;
                    %figure; imagesc(quad1+2*quad2+3*quad3+4*quad4)
%%                    
if roof==true
%% roof error
quad2=imtranslate(quad2, [rooferror/2,0]);
    quad4=imtranslate(quad4,[-rooferror/2,0]);
        quad1=quad4+quad2<1;
        quad3=quad1;
            quad1(Npix/2+1:Npix,:)=0;
                quad3(1:Npix/2,:)=0;
            figure; imagesc(quad1+2*quad2+3*quad3+4*quad4)








end
                    
% Add In Tip/Tilt
focalq1= exp(1i*2*pi*W1).*quad1;%Complex wavefront
focalq2= exp(1i*2*pi*W2).*quad2;%Complex wavefront
focalq3= exp(1i*2*pi*W3).*quad3;%Complex wavefront
focalq4= exp(1i*2*pi*W4).*quad4;%Complex wavefront
pyramidmask= focalq1+focalq2+focalq3+focalq4;
figure; imagesc(angle(pyramidmask));
save('quadpyramidmask.mat','pyramidmask')

if tabletop==true 
flat=pyramidmask(Npix/2, Npix/2);
pyramidmask(Npix/2-tabletoperror/2:Npix/2+tabletoperror/2,Npix/2-tabletoperror/2:Npix/2+tabletoperror/2)=flat;
figure; imagesc(angle(pyramidmask));
figure; imagesc(angle(pyramidmask));
save('quadpyramidmask.mat','pyramidmask')
end
end

%% Continuing the debugging code
pyramid=FTpupil.*pyramidmask;
Pupilpyramid=abs(ifftshift(fft2(ifftshift(pyramid)))).^2;
figure; imagesc(Pupilpyramid); axis equal












end