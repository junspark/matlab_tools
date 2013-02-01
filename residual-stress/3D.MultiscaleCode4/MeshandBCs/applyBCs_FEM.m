function [nps,numels,npsyms,numelsyms]=applyBCs_FEM(R,t,alpha,x,y,z)
% Number of divisions
no_t=size(t,2);
no_R=size(R,2);
no_alpha=size(alpha,2);
%
%
c=0;
nps=0;
% % Bottom surfaces are traction-free
% for ai=1:1:no_alpha-1
%     %for ti=1:1:no_t-1
%         ti=1;
%         for Ri=1:1:no_R-1
%             c=c+1;
%             node1=Ri+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node2=Ri+(no_R*no_alpha*(ti-1))+1+(no_R*(ai-1));
%             node3=Ri+1+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node4=Ri+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             nps(c,1:4)=[ node4 node3 node2 node1 ];
%         end
%     %end
% end
% 
% 
% 
% % top surfaces are traction-free
% for ai=1:1:no_alpha-1
%     %for ti=1:1:no_t-1
%         ti=no_t-1;
%         for Ri=1:1:no_R-1
%             c=c+1;
%             node5=Ri+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node6=Ri+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node7=Ri+no_R+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node8=Ri+no_R+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             nps(c,1:4)=[ node5 node6 node7 node8 ];
%         end
%     %end
% end
% 
% % Surfaces at the outer rim of the disk are traction-free
% for ai=1:1:no_alpha-1
%     for ti=1:1:no_t-1
%         %for Ri=1:1:no_R-1
%             Ri=no_R-1;
%             c=c+1;
%             node2=Ri+(no_R*no_alpha*(ti-1))+1+(no_R*(ai-1));
%             node3=Ri+1+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node6=Ri+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node7=Ri+no_R+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             nps(c,1:4)=[node2 node3 node7 node6];
%         %end
%     end
% end


npsyms=0;
% c=0;
% % In-plane surface tractions vanish at right-hand-side surfaces 
% %for ai=1:1:no_alpha-1
% ai=1;
%     for ti=1:1:no_t-1
%         for Ri=1:1:no_R-1
%             c=c+1;
%             node1=Ri+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node2=Ri+(no_R*no_alpha*(ti-1))+1+(no_R*(ai-1));
%             node5=Ri+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node6=Ri+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             npsyms(c,1:4)=[ node1 node2 node6 node5 ];
%         end
%     end
% %end
% 
% % In-plane surface tractions vanish at left-hand-side surfaces 
% %for ai=1:1:no_alpha-1
% ai=no_alpha-1;
%     for ti=1:1:no_t-1
%         for Ri=1:1:no_R-1
%             c=c+1;
%             node3=Ri+1+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node4=Ri+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
%             node7=Ri+no_R+1+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             node8=Ri+no_R+(no_R*no_alpha*ti)+(no_R*(ai-1));
%             npsyms(c,1:4)=[ node3 node7 node8 node4 ];
%         end
%     end
% %end
%
% In-plane surface: Apply pressure load
pres=linspace(0,400,no_t-1);
for ai=1:1:no_alpha-1
    for ti=1:1:no_t-1
        %for Ri=1:1:no_R-1
            Ri=1;
            c=c+1;
            node1=Ri+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
            node8=Ri+no_R+(no_R*no_alpha*ti)+(no_R*(ai-1));
            node4=Ri+no_R+(no_R*no_alpha*(ti-1))+(no_R*(ai-1));
            node5=Ri+(no_R*no_alpha*ti)+(no_R*(ai-1));
            nps(c,1:4)=[ node1 node5 node8 node4 ];
            % 
        %end
    end
end
%
numelsyms=size(npsyms,1);
numels=size(nps,1);


