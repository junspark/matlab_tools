function T = VectorizedCOBMatrix(R, varargin)
% VECTORIZEDCOBMATRIX - Generate array of 6 x 6 basis transformation
%  matrices for the vectorized tensor representation in 3-D given by:
%         
%   [A] = [A_11, A_12, A_13;A_12, A_22, A_23;A_13, A_23, A_33]
%       |
%       |
%       V
%   {A} = [A_11, A_22, A_33, sqrt(2)*A_23, sqrt(2)*A_13, sqrt(2)*A_12]
%
%   or
%   
%   [A] = [A_11, A_12, A_13; A_12, A_22, A_23; A_13, A_23, A_33]
%       |
%       |
%       V
%   {A} = [A_11, A_22, A_33, sqrt(2)*A_12, sqrt(2)*A_13, sqrt(2)*A_23]
%           
%   where the operation R*A*R' (in tensor notation) is obtained by the
%   matrix-vector product [T]*{A}.
%
%   Usage:
%   
%       T = VectorizedCOBMatrix(R)
%
%   Inputs:
%
%       1. R: The 3 x 3 x n array of rotation matrices representing a
%           change of basis
%
%   Outputs:
%
%       1. T: The 6 x 6 x n array transformation matrices as described above
%

% default options
optcell = {...
    'Order', '11-22-33-23-13-12', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

num_pts = size(R, 3);

% Can't avoid loop easily considering the 3-d input
T = zeros(6, 6, num_pts);
for i = 1:num_pts
    % R(:, :, i) = expm(skew(:, :, i));
    if strcmp(opts.Order, '11-22-33-23-13-12')
        T(1, 1, i) = R(1, 1, i)^2;
        T(1, 2, i) = R(1, 2, i)^2;
        T(1, 3, i) = R(1, 3, i)^2;
        T(1, 4, i) = sqrt(2)*R(1, 2, i)*R(1, 3, i);
        T(1, 5, i) = sqrt(2)*R(1, 1, i)*R(1, 3, i);
        T(1, 6, i) = sqrt(2)*R(1, 1, i)*R(1, 2, i);
        
        T(2, 1, i) = R(2, 1, i)^2;
        T(2, 2, i) = R(2, 2, i)^2;
        T(2, 3, i) = R(2, 3, i)^2;
        T(2, 4, i) = sqrt(2)*R(2, 2, i)*R(2, 3, i);
        T(2, 5, i) = sqrt(2)*R(2, 1, i)*R(2, 3, i);
        T(2, 6, i) = sqrt(2)*R(2, 1, i)*R(2, 2, i);
        
        T(3, 1, i) = R(3, 1, i)^2;
        T(3, 2, i) = R(3, 2, i)^2;
        T(3, 3, i) = R(3, 3, i)^2;
        T(3, 4, i) = sqrt(2)*R(3, 2, i)*R(3, 3, i);
        T(3, 5, i) = sqrt(2)*R(3, 1, i)*R(3, 3, i);
        T(3, 6, i) = sqrt(2)*R(3, 1, i)*R(3, 2, i);
        
        T(4, 1, i) = sqrt(2)*R(2, 1, i)*R(3, 1, i);
        T(4, 2, i) = sqrt(2)*R(2, 2, i)*R(3, 2, i);
        T(4, 3, i) = sqrt(2)*R(2, 3, i)*R(3, 3, i);
        T(4, 4, i) = R(2, 3, i)*R(3, 2, i) + R(2, 2, i)*R(3, 3, i);
        T(4, 5, i) = R(2, 3, i)*R(3, 1, i) + R(2, 1, i)*R(3, 3, i);
        T(4, 6, i) = R(2, 2, i)*R(3, 1, i) + R(2, 1, i)*R(3, 2, i);
        
        T(5, 1, i) = sqrt(2)*R(1, 1, i)*R(3, 1, i);
        T(5, 2, i) = sqrt(2)*R(1, 2, i)*R(3, 2, i);
        T(5, 3, i) = sqrt(2)*R(1, 3, i)*R(3, 3, i);
        T(5, 4, i) = R(1, 3, i)*R(3, 2, i) + R(1, 2, i)*R(3, 3, i);
        T(5, 5, i) = R(1, 3, i)*R(3, 1, i) + R(1, 1, i)*R(3, 3, i);
        T(5, 6, i) = R(1, 2, i)*R(3, 1, i) + R(1, 1, i)*R(3, 2, i);
        
        T(6, 1, i) = sqrt(2)*R(1, 1, i)*R(2, 1, i);
        T(6, 2, i) = sqrt(2)*R(1, 2, i)*R(2, 2, i);
        T(6, 3, i) = sqrt(2)*R(1, 3, i)*R(2, 3, i);
        T(6, 4, i) = R(1, 3, i)*R(2, 2, i) + R(1, 2, i)*R(2, 3, i);
        T(6, 5, i) = R(1, 1, i)*R(2, 3, i) + R(1, 3, i)*R(2, 1, i);
        T(6, 6, i) = R(1, 2, i)*R(2, 1, i) + R(1, 1, i)*R(2, 2, i);
    elseif strcmp(opts.Order, '11-22-33-12-13-23')
        T(1, 1, i) = R(1, 1, i)^2;
        T(1, 2, i) = R(1, 2, i)^2;
        T(1, 3, i) = R(1, 3, i)^2;
        T(1, 4, i) = sqrt(2)*R(1, 1, i)*R(1, 2, i);
        T(1, 5, i) = sqrt(2)*R(1, 1, i)*R(1, 3, i);
        T(1, 6, i) = sqrt(2)*R(1, 2, i)*R(1, 3, i);
        
        T(2, 1, i) = R(2, 1, i)^2;
        T(2, 2, i) = R(2, 2, i)^2;
        T(2, 3, i) = R(2, 3, i)^2;
        T(2, 4, i) = sqrt(2)*R(2, 1, i)*R(2, 2, i);
        T(2, 5, i) = sqrt(2)*R(2, 1, i)*R(2, 3, i);
        T(2, 6, i) = sqrt(2)*R(2, 2, i)*R(2, 3, i);
        
        T(3, 1, i) = R(3, 1, i)^2;
        T(3, 2, i) = R(3, 2, i)^2;
        T(3, 3, i) = R(3, 3, i)^2;
        T(3, 4, i) = sqrt(2)*R(3, 1, i)*R(3, 2, i);
        T(3, 5, i) = sqrt(2)*R(3, 1, i)*R(3, 3, i);
        T(3, 6, i) = sqrt(2)*R(3, 2, i)*R(3, 3, i);
        
        T(4, 1, i) = sqrt(2)*R(2, 1, i)*R(1, 1, i);
        T(4, 2, i) = sqrt(2)*R(1, 2, i)*R(2, 2, i);
        T(4, 3, i) = sqrt(2)*R(1, 3, i)*R(2, 3, i);
        T(4, 4, i) = (R(1, 1, i)*R(2, 2, i) + R(1, 2, i)*R(2, 1, i));
        T(4, 5, i) = (R(1, 1, i)*R(2, 3, i) + R(1, 3, i)*R(2, 1, i));
        T(4, 6, i) = (R(1, 2, i)*R(2, 3, i) + R(1, 3, i)*R(2, 2, i));
        
        T(5, 1, i) = sqrt(2)*R(1, 1, i)*R(3, 1, i);
        T(5, 2, i) = sqrt(2)*R(1, 2, i)*R(3, 2, i);
        T(5, 3, i) = sqrt(2)*R(1, 3, i)*R(3, 3, i);
        T(5, 4, i) = (R(1, 1, i)*R(3, 2, i) + R(1, 2, i)*R(3, 1, i));
        T(5, 5, i) = (R(1, 1, i)*R(3, 3, i) + R(1, 3, i)*R(3, 1, i));
        T(5, 6, i) = (R(1, 2, i)*R(3, 3, i) + R(1, 3, i)*R(3, 2, i));
        
        T(6, 1, i) = sqrt(2)*R(2, 1, i)*R(3, 1, i);
        T(6, 2, i) = sqrt(2)*R(2, 2, i)*R(3, 2, i);
        T(6, 3, i) = sqrt(2)*R(2, 3, i)*R(3, 3, i);
        T(6, 4, i) = (R(2, 1, i)*R(3, 2, i) + R(2, 2, i)*R(3, 1, i));
        T(6, 5, i) = (R(2, 1, i)*R(3, 3, i) + R(2, 3, i)*R(3, 1, i));
        T(6, 6, i) = (R(2, 2, i)*R(3, 3, i) + R(2, 3, i)*R(3, 2, i));
    else
        disp('order not supported');
        return
    end
end

% function T = VectorizedCOBMatrix(R)
% % VECTORIZEDCOBMATRIX - Generate array of 6 x 6 basis transformation
% %  matrices for the vectorized tensor representation in 3-D given by:
% %         
% %   [A] = [A_11, A_12, A_13; A_12, A_22, A_23; A_13, A_23, A_33]
% %       |
% %       |
% %       V
% %   {A} = [A_11, A_22, A_33, sqrt(2)*A_12, sqrt(2)*A_13, sqrt(2)*A_23]
% %           
% %   where the operation R*A*R' (in tensor notation) is obtained by the
% %   matrix-vector product [T]*{A}.
% %
% %   Usage:
% %   
% %       T = VectorizedCOBMatrix(R)
% %
% %   Inputs:
% %
% %       1. R: The 3 x 3 x n array of rotation matrices representing a
% %           change of basis
% %
% %   Outputs:
% %
% %       1. T: The 6 x 6 x n array transformation matrices as described above
% %   
% num_pts = size(R, 3);
% 
% T = zeros(6, 6, num_pts);
% for i = 1:num_pts
%     T(1, 1, i) = R(1, 1, i)^2;    
%     T(1, 2, i) = R(1, 2, i)^2;
%     T(1, 3, i) = R(1, 3, i)^2;
%     T(1, 4, i) = sqrt(2)*R(1, 1, i)*R(1, 2, i);
%     T(1, 5, i) = sqrt(2)*R(1, 1, i)*R(1, 3, i);
%     T(1, 6, i) = sqrt(2)*R(1, 2, i)*R(1, 3, i);
%     
%     T(2, 1, i) = R(2, 1, i)^2;
%     T(2, 2, i) = R(2, 2, i)^2;
%     T(2, 3, i) = R(2, 3, i)^2;
%     T(2, 4, i) = sqrt(2)*R(2, 1, i)*R(2, 2, i);
%     T(2, 5, i) = sqrt(2)*R(2, 1, i)*R(2, 3, i);
%     T(2, 6, i) = sqrt(2)*R(2, 2, i)*R(2, 3, i);
%     
%     T(3, 1, i) = R(3, 1, i)^2;
%     T(3, 2, i) = R(3, 2, i)^2;
%     T(3, 3, i) = R(3, 3, i)^2;
%     T(3, 4, i) = sqrt(2)*R(3, 1, i)*R(3, 2, i);
%     T(3, 5, i) = sqrt(2)*R(3, 1, i)*R(3, 3, i);
%     T(3, 6, i) = sqrt(2)*R(3, 2, i)*R(3, 3, i);
%     
%     T(4, 1, i) = sqrt(2)*R(2, 1, i)*R(1, 1, i);
%     T(4, 2, i) = sqrt(2)*R(1, 2, i)*R(2, 2, i);
%     T(4, 3, i) = sqrt(2)*R(1, 3, i)*R(2, 3, i);
%     T(4, 4, i) = (R(1, 1, i)*R(2, 2, i) + R(1, 2, i)*R(2, 1, i));
%     T(4, 5, i) = (R(1, 1, i)*R(2, 3, i) + R(1, 3, i)*R(2, 1, i));
%     T(4, 6, i) = (R(1, 2, i)*R(2, 3, i) + R(1, 3, i)*R(2, 2, i));
%     
%     T(5, 1, i) = sqrt(2)*R(1, 1, i)*R(3, 1, i);
%     T(5, 2, i) = sqrt(2)*R(1, 2, i)*R(3, 2, i);
%     T(5, 3, i) = sqrt(2)*R(1, 3, i)*R(3, 3, i);
%     T(5, 4, i) = (R(1, 1, i)*R(3, 2, i) + R(1, 2, i)*R(3, 1, i));
%     T(5, 5, i) = (R(1, 1, i)*R(3, 3, i) + R(1, 3, i)*R(3, 1, i));
%     T(5, 6, i) = (R(1, 2, i)*R(3, 3, i) + R(1, 3, i)*R(3, 2, i));
%     
%     T(6, 1, i) = sqrt(2)*R(2, 1, i)*R(3, 1, i);
%     T(6, 2, i) = sqrt(2)*R(2, 2, i)*R(3, 2, i);
%     T(6, 3, i) = sqrt(2)*R(2, 3, i)*R(3, 3, i);
%     T(6, 4, i) = (R(2, 1, i)*R(3, 2, i) + R(2, 2, i)*R(3, 1, i));
%     T(6, 5, i) = (R(2, 1, i)*R(3, 3, i) + R(2, 3, i)*R(3, 1, i));
%     T(6, 6, i) = (R(2, 2, i)*R(3, 3, i) + R(2, 3, i)*R(3, 2, i));
% end