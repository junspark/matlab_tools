ri=6;       % mm
ro=15;      % mm
pi=400;     % MPa
as_sig=zeros(6,numnp);   
% Calculate the stresses at the each node
for j=1:1:numnp
   r2=x(j)^2+y(j)^2;
   th=atan2(y(j),x(j));
   s_r=pi*ri^2*(1-(ro^2/r2))/(ro^2-ri^2);
   s_h=pi*ri^2*(1+(ro^2/r2))/(ro^2-ri^2);
   S=[s_r 0  0;
       0 s_h 0;
       0  0  0];
   g=[  cos(th) -sin(th) 0;
        sin(th) cos(th) 0;
        0 0 1];
   % Transform the stress
   Str=g*S*g';
   % Order of stress components is important!!!
   as_sig(1:6,j)=[Str(1,1) Str(2,2) Str(3,3) Str(1,2) Str(2,3) Str(3,1)];
end

