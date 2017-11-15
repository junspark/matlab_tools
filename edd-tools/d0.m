% [d-spacing] = d0('system',a,'# of peaks to show');          
% calculate d-spacing for a given lattice (cubic system only)
%
% system -> 'fcc','bcc','sc','dia','hcp'
% required sub-function: calchkl2.m
%                        d0calc.m
%                                            AC@Feb.2012

function [d_,lat_]=d0(system,a,range,output)
%hkl;
if ~exist('a','var')||isempty(a); a=3.6;end  % default a = 3.6 Angstron
if ~exist('range','var')||isempty(range); range=5;end  % output 5 peaks in default
if ~exist('output','var')||isempty(output); output=1;end % show d-spacing on screen(default: 'on')
c = 1.633*a;
sys='cubic';
switch upper(system)
    case 'SC';      [lat,F,s2]=calchkl2(300,'SC'); 
    case 'FCC';     [lat,F,s2]=calchkl2(300,'FCC'); 
    case 'BCC';     [lat,F,s2]=calchkl2(300,'BCC');
    case 'DIA';     [lat,F,s2]=calchkl2(300,'DIA');
    case 'HCP';     [lat,F,s2]=calchkl2(300,'HCP'); sys='hexag'; if size(a,2)==3;c=a(3);elseif size(a,2)==2; c=a(2);end
    case 'CERIA';   [lat,F,s2]=calchkl2(300,'FCC'); a=5.41165;
    case 'LAB6';    [lat,F,s2]=calchkl2(300,'SC');  a=4.15689;
    case 'ALUMINA'; [lat,F,s2]=calchkl2(300,'HCP'); a=4.7590914; c=12.991779; sys='hexag';
    case 'GRAPHITE'; [lat,F,s2]=calchkl2(300,'HCP'); a=2.47; c=6.79; sys='hexag';
end
d=zeros(range+30,1);
for i=1:range+30;d(i,1)=d0calc(sys,lat(i,1),lat(i,2),lat(i,3),a(1),a(1),c);end
[~,k]=sort(d,'descend');d=d(k);lat=lat(k,:);F=F(k);s2=s2(k); d_=d(1:range); lat_=lat(1:range,:);

if output
fprintf('Allowed diffraction for %s lattic: \n',system);
   switch upper(system)
       case {'HCP','ALUMINA','GRAPHITE'}
           fprintf(' a = %s %s; c = %s %s \n',num2str(a(1)),char(197),num2str(c),char(197));
       otherwise
           fprintf(' a = %s %s;\n',num2str(a(1)),char(197));
   end
fprintf(' (hkl) ; s2 ; F2 ; d-spac (%s)\n',char(197));
for i=1:range
    fprintf('(%s) : %s : %s : %s %s \n',num2str(lat(i,:)),num2str(s2(i)),num2str(F(i).^2),num2str(d(i)),char(197)); 
end
end