function [np] = findconnectivity(R,t,alpha)
%
noR = size(R,2);
noT = size(t,2);
noA = size(alpha,2);
c   = 1;
for i = 1:1:noR-1
    for j = 1:1:noA-1
        for k = 1:1:noT-1
            n1  = k + noT*(j-1) + noA*noT*(i-1);
            n2  = k + noT*j + noA*noT*(i-1);
            n3  = k + noT*j + noA*noT*i;
            n4  = k + noT*(j-1) + noA*noT*i;
            n5  = k + 1 + noT*(j-1) + noA*noT*(i-1);
            n6  = k + 1 + noT*j + noA*noT*(i-1);
            n7  = k + 1 + noT*j + noA*noT*i;
            n8  = k + 1 + noT*(j-1) + noA*noT*i;
            
            np(c,1:8)   = [n1 n2 n3 n4 n5 n6 n7 n8];
            c           = c + 1;
        end
    end
end





% 
% 
% no_t=size(t,2);
% no_R=size(R,2);
% no_alpha=size(alpha,2);
% c=0;
% np=zeros((no_alpha-1)*(no_t-1)*(no_R-1),8);
% for ai=1:1:no_alpha-1
%     for ti=1:1:no_t-1
%         for Ri=1:1:no_R-1
%             
%             c=c+1;
%             node1=Ri+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node2=Ri+(no_R*no_alpha*(ti-1))+1+(no_R*(ai-1));
%             node3=Ri+1+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node4=Ri+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node5=Ri+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node6=Ri+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node7=Ri+no_R+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node8=Ri+no_R+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             np(c,1:8)=[ node1 node2 node3 node4 ...
%                         node5 node6 node7 node8];
% %                     original
% %                     np(c,1:8)=[ node1 node4 node3 node2 ...
% %                         node5 node8 node7 node6];
%         end
%     end
% end