function outvec = Decomp2DevSphDecomp(invec)
% usage:
%  mat = Decomp2Matrix(vec)
%  
%   vec is 6 x n
%   mat is 3 x 3 x n (symm) 
num_vec = size(vec, 2);
%
onebyroot2 = 1/sqrt(2);
onebyroot3 = 1/sqrt(3);
onebyroot6 = 1/sqrt(6);
T = [onebyroot2   -onebyroot2              0   0   0   0;
    -onebyroot6   -onebyroot6   2*onebyroot6   0   0   0;
              0             0              0   1   0   0;
              0             0              0   0   1   0;
              0             0              0   0   0   1;
     onebyroot3    onebyroot3     onebyroot3   0   0   0;
 ];
 
outvec = T*invec;