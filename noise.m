function n = noise(im, nph, rdns)

% image_out = noise(image_in, nph, rdns)
% Add photon and read noise to an image
%  - Poisson noise appropriate to a total photon count of nph per frame
%  - Gaussian read noise with rdns rms counts per pixel
% Output counts are in units of photons.
% Setting nph = -1 treats counts in image_in as photons directly

n = zeros(size(im));

nframe = size(im,3);
do_norm = (nph~=-1);
ms = ones(size(im,1),size(im,2))*rdns^2;

if ~do_norm
    norm = ones(nframe,1);
else
    norm = nph./squeeze(sum(sum(im,1),2));
end

progressbar(0);
for k=1:nframe
   n(:,:,k) = randp(im(:,:,k)*norm(k))+round(randgauss(0,ms));
   progressbar(k/nframe);
end
