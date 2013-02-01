function [Rnodes,tnodes,alphanodes] = findnodalvalues(R,t,dt,alpha,dalpha)

% Calculations start HERE!
% _______________________________
%
% Number of divisions
no_t=size(t,2);
no_R=size(R,1);
no_alpha=size(alpha,2);
%
% Create the nodes that are centered at the given 
% Radial locations of the nodes
c=0;
for i=1:1:no_R
    %
    if i==1
        dR=R(2)-R(1);
        c=c+1;
        Rnodes(c)=R(i)-(dR/2);
    end
    if i==no_R
        dR=R(i)-R(i-1);
    else
        dR=R(i+1)-R(i);
    end
    c=c+1;
    Rnodes(c)=R(i)+(dR/2);    
end
% thickness locations of the nodes
c=0;
if no_t>1
    for i=1:1:no_t
        %
        if i==1
            dt=t(2)-t(1);
            c=c+1;
            tnodes(c)=t(i)-(dt/2);
        end
        if i==no_t
            dt=t(i)-t(i-1);
        else
            dt=t(i+1)-t(i);
        end
        c=c+1;
        tnodes(c)=t(i)+(dt/2);    
    end
else
   tnodes(1)=t-dt;
   tnodes(2)=t+dt;
end
% Alpha locations of nodes
% If there is only one alpha
if no_alpha==1
    alphanodes(1)=alpha-dalpha;
    alphanodes(2)=alpha+dalpha;
else
    c=0;
    for i=1:1:no_alpha
        if i==1
            dalpha=alpha(2)-alpha(1);
            c=c+1;
            alphanodes(c)=alpha(i)-(dalpha/2);
        end
        if i==no_alpha
            dalpha=alpha(i)-alpha(i-1);
        else
            dalpha=alpha(i+1)-alpha(i);
        end
        c=c+1;
        alphanodes(c)=alpha(i)+(dalpha/2);    
    end
end
