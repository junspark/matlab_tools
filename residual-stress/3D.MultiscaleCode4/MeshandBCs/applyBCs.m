function [nps,numels,npsyms,numelsyms]=applyBCs(R,t,alpha)
noR = size(R,2);
noT = size(t,2);
noA = size(alpha,2);

% % % % traction-free
nps     = [];
% % % Bottom surface : traction-free
for i = 1:1:noR-1
    for j = 1:1:noA-1
        for k = 1
            n1  = k + noT*(j-1) + noA*noT*(i-1);
            n2  = k + noT*j + noA*noT*(i-1);
            n3  = k + noT*j + noA*noT*i;
            n4  = k + noT*(j-1) + noA*noT*i;
            
            nps = [nps; n1 n4 n3 n2];
        end
    end
end
% % % % Outer rim : traction-free
for i = noR
    for j = 1:1:noA-1
        for k = 1:1:noT-1
            n1  = k + noT*(j-1) + noA*noT*(i-1);
            n2  = k + noT*j + noA*noT*(i-1);
            n5  = k + 1 + noT*(j-1) + noA*noT*(i-1);
            n6  = k + 1 + noT*j + noA*noT*(i-1);
            
            nps = [nps; n1 n5 n6 n2];
        end
    end
end
% % % % % Top surface : traction-free
for i = 1:1:noR-1
    for j = 1:1:noA-1
        for k = noT
            n5  = k + noT*(j-1) + noA*noT*(i-1);
            n6  = k + noT*j + noA*noT*(i-1);
            n7  = k + noT*j + noA*noT*i;
            n8  = k + noT*(j-1) + noA*noT*i;
            
            nps = [nps; n5 n6 n7 n8];
        end
    end
end
numels  = size(nps,1);

% % % % In-plane surface tractions
npsyms  = [];
% % % % Inner rim : In-plane surface tractions vanish
for i = 1
    for j = 1:1:noA-1
        for k = 1:1:noT-1
            n1  = k + noT*(j-1) + noA*noT*(i-1);
            n2  = k + noT*j + noA*noT*(i-1);
            n5  = k + 1 + noT*(j-1) + noA*noT*(i-1);
            n6  = k + 1 + noT*j + noA*noT*(i-1);
            
            npsyms  = [npsyms; n1 n2 n6 n5];
        end
    end
end
% % % % LHS surface : In-plane surface tractions vanish
for i = 1:1:noR-1
    for j = 1
        for k = 1:1:noT-1
            n1  = k + noT*(j-1) + noA*noT*(i-1);
            n4  = k + noT*(j-1) + noA*noT*i;
            n5  = k + 1 + noT*(j-1) + noA*noT*(i-1);
            n8  = k + 1 + noT*(j-1) + noA*noT*i;
            
            npsyms  = [npsyms; n1 n5 n8 n4];
        end
    end
end
% % % % RHS surface : In-plane surface tractions vanish
for i = 1:1:noR-1
    for j = noA-1
        for k = 1:1:noT-1
            n2  = k + noT*j + noA*noT*(i-1);
            n3  = k + noT*j + noA*noT*i;
            n6  = k + 1 + noT*j + noA*noT*(i-1);
            n7  = k + 1 + noT*j + noA*noT*i;
            
            npsyms  = [npsyms; n2 n3 n7 n6];
        end
    end
end
numelsyms   = size(npsyms,1);