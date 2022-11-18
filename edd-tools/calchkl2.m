function [hkl,F,s2] = calchkl2(s2max,Struc)
% [hkl,S,s2] = calchkl2(s2max,Struc)
% where Struc = 'BCC', 'FCC', 'DIA' or 'HCP'
%           S = the geometric factor
%          s2 = h^2+k^2+l^2
%                              AC@Mar.2012

if nargin ~= 2
    fprintf('[hkl,S,s2] = calchkl2(s2max,Struc)\n');
    return;
end
hkl = []; Struc = upper(Struc); tp = 2*pi*sqrt(-1);

if strcmp(Struc,'HCP')
  for h = 0:floor(sqrt(s2max)),
   More = 1; k = 0;
   while More,
     if  k > h, More = 0;
     else
       if h==0&&k==0;l=1;else l = 0;end
       while h^2+h*k+k^2 <= s2max && l <= 10
         hkl = [hkl; h k l]; l = l+1;
       end;
     end;
     k = k+1;
   end;
  end;
else
  for h = 1:floor(sqrt(s2max)),
   More = 1; k = 0;
   while More,
     if  k > h, More = 0;
     else
       l = 0;
       while h^2+k^2+l^2 <= s2max && l <= k,
         hkl = [hkl; h k l]; l = l+1;
       end;
     end;
     k = k+1;
   end;
  end;
end


switch Struc
    case ('BCC')
        R = [0 0 0; 1 1 1]'/2; b = [0 0 0]';
    case('FCC')
        R = [0 0 0; 1 1 0; 1 0 1; 0 1 1]'/2; b = [0 0 0]';
    case('DIA')
        R = [0 0 0; 1 1 0; 1 0 1; 0 1 1]'/2; b = [0 0 0; 1 1 1]'/4;
    case('AL2O3')
        R = [0 0 0; 1 -sqrt(3) 0;1 sqrt(3) 0;0 0 5.45955]'/2;
        b1= [0 0 -4.57; 0 0 -1.92; 0 0 1.92; 0 0 4.57]';
        b2= [0.73 1.26 -3.25; -1.46 0 -3.25; 0.73 -1.26 -3.25; -0.73 -1.26 3.25; 1.46 0 3.25; -0.73 1.26 3.25]';    
    case('SC')
        R = [0 0 0]'; b = [0 0 0]';
    case('HCP')
        R = [0 0 0; 1/3 2/3 1/2]'; b = [0 0 0]';
end;
F = abs(sum(exp(tp*hkl*R),2).*sum(exp(tp*hkl*b),2));
%F =abs(sum(exp(tp*hkl*R),2).*sum(exp(tp*hkl*b1),2)+sum(exp(tp*hkl*R),2).*sum(exp(tp*hkl*b2),2));
s2= sum(hkl.^2,2);
i = find(1e-6 < F); hkl = hkl(i,:); F = F(i); s2 = s2(i);
if strcmp(Struc,'HCP')
    s2 = hkl(:,1).^2 + hkl(:,1).*hkl(:,2) + hkl(:,2).^2;
    %[~,i] = sort(hkl(:,1).^2+hkl(:,2).^2+hkl(:,1).*hkl(:,2)); hkl = hkl(i,:); F = F(i); s2 = s2(i);
    [~,i] = sort(sqrt(4/3/1^2*(hkl(:,1).^2+hkl(:,1).*hkl(:,2)+hkl(:,2).^2)+hkl(:,3).^2/1.633^2)); hkl = hkl(i,:); F = F(i); s2 = s2(i);
else
    [~,i] = sort(sum(hkl.^2,2)); hkl = hkl(i,:); F = F(i); s2 = s2(i);
end


