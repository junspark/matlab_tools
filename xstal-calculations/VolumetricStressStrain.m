function tout = VolumetricStressStrain(tin)

if (size(tin,1) == 3) && (size(tin,2) == 3)
    h   = tin(1,1) + tin(2,2) + tin(3,3);
    h   = h/3;
    tout    = eye(3,3).*h;
else
    disp('input matrix size incorrect');
end