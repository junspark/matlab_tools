function vm = VMStressStrain(tin)

if (size(tin,1) == 3) && (size(tin,2) == 3)
    vm  = (tin(1,1) - tin(2,2))^2 + (tin(2,2) - tin(3,3))^2 + (tin(3,3) - tin(1,1))^2 + 6*(tin(1,2)^2 + tin(1,3)^2 + tin(2,3)^2);
    vm  = sqrt(vm/2);
else
    disp('input matrix size incorrect');
end