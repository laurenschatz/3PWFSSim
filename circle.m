function c = circle(s,d,varargin)

% c = circle(s, d)
% c = circle(s, d, [x y])
% Create an array of zeros of size s by v with a circle of diameter d filled with ones
% either in the center of the array or with center at [x y]

if nargin>2
   cen = varargin{1};
   if (size(cen,1)~=1)|(size(cen,2)~= 2)
      fprintf('Error in center definition\n')
      return
   end
else
   m = bitshift(s,-1)+1;
   cen = [m m];
end

c = zeros(s);
[x,y] = meshgrid(-cen(1)+1:s-cen(1),-cen(2)+1:s-cen(2));
r = sqrt(double(x.^2+y.^2));
c(r<=d/2) = 1;
