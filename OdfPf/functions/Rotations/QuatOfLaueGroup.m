function qsym = QuatOfLaueGroup(schoenfliesTag)
% QUATOFLAUEGROUP - Generate quaternion representation for the specified
% Laue group.
%
% USAGE:
%
%      qsym = QuatOfLaueGroup(schoenfliesTag)
%
% INPUTS:
%
%      1) schoenfliesTag 1 x 1, a case-insensitive string representing the
%      Schoenflies symbol for the desired Laue group.  The 11 available 
%      choices are:
%
%           Class           Symbol      n
%          -------------------------------
%           Triclinic       Ci (S2)     1
%           Monoclinic      C2h         2
%           Orthorhombic    D2h (Vh)    4
%           Tetragonal      C4h         4
%                           D4h         8
%           Trigonal        C3i (S6)    3
%                           D3d         6
%           Hexagonal       C6h         6
%                           D6h         12
%           Cubic           Th          12
%                           Oh          24
%
% OUTPUTS:
%
%      1) qsym is 4 x n, where n is the number of symmetries (dep. on group 
%      -- see INPUTS list above).
%
% NOTES:
%
%      *) The conventions used for assigning a RHON basis, {x1, x2, x3}, to
%      each point group are consistent with those published in Appendix B
%      of [1].
%
%      *) This routine replaces the older (individual) symmetry group
%      routines, which are still available in the distribution, but made
%      redundant on 12/1/2006.
%
% REFERENCES:
%
%      [1] Nye, J. F., "Physical Properties of Crystals: Their
%      Representation by Tensors and Matrices", Oxford University Press,
%      1985. ISBN 0198511655

if strcmpi(schoenfliesTag, 'Ci') || strcmpi(schoenfliesTag, 'S2')
  % TRICLINIC
  angleAxis = [...
    0.0   1   0   0;			% identity
    ]';
elseif strcmpi(schoenfliesTag, 'C2h')
  % MONOCLINIC
  angleAxis = [...
    0.0   1   0   0;			% identity
    pi    0   1   0;			% twofold about 010 (x2)
    ]';
elseif strcmpi(schoenfliesTag, 'D2h') || strcmpi(schoenfliesTag, 'Vh')
  % ORTHORHOMBIC
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi       1    0    0;		% twofold about 100
    pi       0    1    0;		% twofold about 010
    pi       0    0    1;		% twofold about 001
    ]';
elseif strcmpi(schoenfliesTag, 'C4h')
  % TETRAGONAL (LOW)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi*0.5   0    0    1; 		% fourfold about 001 (x3)
    pi       0    0    1;		%
    pi*1.5   0    0    1;		%
    ]';
elseif strcmpi(schoenfliesTag, 'D4h')
  % TETRAGONAL (HIGH)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi*0.5   0    0    1; 		% fourfold about  0  0  1 (x3)
    pi       0    0    1;		%
    pi*1.5   0    0    1;		%
    pi       1    0    0; 		% twofold about  1  0  0 (x1)
    pi       0    1    0; 		% twofold about  0  1  0 (x2)
    pi       1    1    0; 		% twofold about  1  1  0
    pi      -1    1    0;		% twofold about -1  1  0
    ]';
elseif strcmpi(schoenfliesTag, 'C3i') || strcmpi(schoenfliesTag, 'S6')
  % TRIGONAL (LOW)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi*2/3   0    0    1; 		% threefold about 0001 (x3, c)
    pi*4/3   0    0    1;		%
    ]';
elseif strcmpi(schoenfliesTag, 'D3d')
  % TRIGONAL (HIGH)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi*2/3   0    0    1; 		% threefold about 0001 (x3, c)
    pi*4/3   0    0    1;		%
    pi       1    0    0;		% twofold about  2 -1 -1  0 (x1, a1)
    pi      -0.5  sqrt(3)/2  0;	% twofold about -1  2 -1  0 (a2)
    pi      -0.5 -sqrt(3)/2  0;	% twofold about -1 -1  2  0 (a3)
    ]';
elseif strcmpi(schoenfliesTag, 'C6h')
  % HEXAGONAL (LOW)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi/3     0    0    1; 		% sixfold about 0001 (x3, c)
    pi*2/3   0    0    1;
    pi       0    0    1;
    pi*4/3   0    0    1;
    pi*5/3   0    0    1;
    ]';
elseif strcmpi(schoenfliesTag, 'D6h')
  % HEXAGONAL (HIGH)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi/3     0    0    1; 		% sixfold about 0001 (x3, c)
    pi*2/3   0    0    1;
    pi       0    0    1;
    pi*4/3   0    0    1;
    pi*5/3   0    0    1;
    pi       1    0    0;		% twofold about  2 -1 -1  0 (x1, a1)
    pi      -0.5  sqrt(3)/2  0;	% twofold about -1  2 -1  0 (a2)
    pi      -0.5 -sqrt(3)/2  0;	% twofold about -1 -1  2  0 (a3)
    pi       sqrt(3)/2  0.5  0;	% twofold about  1  0 -1  0
    pi       0    1    0;		% twofold about -1  1  0  0 (x2)
    pi      -sqrt(3)/2  0.5  0;	% twofold about  0 -1  1  0
    ]';
elseif strcmpi(schoenfliesTag, 'Th')
  % CUBIC (LOW)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi       1    0    0; 		% twofold about  1  0  0 (x1)
    pi       0    1    0; 		% twofold about  0  1  0 (x2)
    pi       0    0    1; 		% twofold about  0  0  1 (x3)
    pi*2/3   1    1    1; 		% threefold about  1  1  1
    pi*4/3   1    1    1;		%
    pi*2/3  -1    1    1; 		% threefold about -1  1  1
    pi*4/3  -1    1    1;		%
    pi*2/3  -1   -1    1; 		% threefold about -1 -1  1
    pi*4/3  -1   -1    1;		%
    pi*2/3   1   -1    1; 		% threefold about  1 -1  1
    pi*4/3   1   -1    1;		%
    ]';
elseif strcmpi(schoenfliesTag, 'Oh')
  % CUBIC (HIGH)
  angleAxis = [...
    0.0      1    0    0; 		% identity
    pi*0.5   1    0    0; 		% fourfold about  1  0  0 (x1)
    pi       1    0    0;		%
    pi*1.5   1    0    0;		%
    pi*0.5   0    1    0; 		% fourfold about  0  1  0 (x2)
    pi       0    1    0;		%
    pi*1.5   0    1    0;		%
    pi*0.5   0    0    1; 		% fourfold about  0  0  1 (x3)
    pi       0    0    1;		%
    pi*1.5   0    0    1;		%
    pi*2/3   1    1    1; 		% threefold about  1  1  1
    pi*4/3   1    1    1;		%
    pi*2/3  -1    1    1; 		% threefold about -1  1  1
    pi*4/3  -1    1    1;		%
    pi*2/3  -1   -1    1; 		% threefold about -1 -1  1
    pi*4/3  -1   -1    1;		%
    pi*2/3   1   -1    1; 		% threefold about  1 -1  1
    pi*4/3   1   -1    1;		%
    pi       1    1    0; 		% twofold about  1  1  0
    pi      -1    1    0;		% twofold about -1  1  0
    pi       1    0    1;		% twofold about  1  0  1
    pi       0    1    1;		% twofold about  0  1  1
    pi      -1    0    1;		% twofold about -1  0  1
    pi       0   -1    1;		% twofold about  0 -1  1
    ]';
end
%
Angle = angleAxis(1,:);
Axis  = angleAxis(2:4,:);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
qsym = QuatOfAngleAxis(Angle, Axis);

%%%
%%%%%%
%%%%%%%%% NESTED FUNCTIONS
%%%%%%
%%%

  function quat = QuatOfAngleAxis(angle, rotaxis)
    % QUATOFANGLEAXIS - Quaternion of rotation of angle about axis.
    %
    %  quat = QuatOfAngleAxis(angle, rotaxis)
    %
    %  angle   is 1 x n, the list of rotation angles.
    %  rotaxis is 3 x n, the list of rotation axes, which need not
    %          be normalized (e.g. [1 1 1]'), but must be nonzero
    %
    % quat is 4 x n, the quaternion representations of the given
    %         rotations.  The first component of quat is nonnegative.
    %
    halfangle = 0.5*angle;
    cphiby2 = cos(halfangle);
    sphiby2 = sin(halfangle);
    %
    rescale = sphiby2 ./sqrt(dot(rotaxis, rotaxis, 1));
    %
    quat = [cphiby2; repmat(rescale, [3 1]) .* rotaxis ] ;
    %
    q1negative = (quat(1,:) < 0);
    quat(:, q1negative) = -1*quat(:, q1negative);

  end

%%%
%%%%%%
%%%%%%%%% NESTED FUNCTIONS
%%%%%%
%%%

end