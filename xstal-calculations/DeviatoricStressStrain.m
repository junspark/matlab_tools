function tout = DeviatoricStressStrain(tin)

if (size(tin,1) == 3) && (size(tin,2) == 3)
    h   = tin(1,1) + tin(2,2) + tin(3,3);
    h   = h/3;
    h   = eye(3,3).*h;
    
    tout    = tin - h;
else
    disp('input matrix size incorrect');
end