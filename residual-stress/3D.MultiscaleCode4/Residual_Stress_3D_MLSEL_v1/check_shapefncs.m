clc
clear
close

% Elemen type
meltyp=20;

% % Enter the mapped coordinates
% q1=1;q2=1;q3=1;



% Enter the nodal coordinates
% Nodal coordinates (27 node brick)
x(1)=0;y(1)=0;z(1)=0;
x(2)=1;y(2)=0;z(2)=0;
x(3)=1;y(3)=1;z(3)=0;
x(4)=0;y(4)=1;z(4)=0;
x(5)=0;y(5)=0;z(5)=1;
x(6)=1;y(6)=0;z(6)=1;
x(7)=1;y(7)=1;z(7)=1;
x(8)=0;y(8)=1;z(8)=1;
x(9)=0.5;y(9)=0;z(9)=0;
x(10)=1;y(10)=0.5;z(10)=0;
x(11)=0.5;y(11)=1;z(11)=0;
x(12)=0;y(12)=0.5;z(12)=0;
x(13)=0.5;y(13)=0;z(13)=1;
x(14)=1;y(14)=0.5;z(14)=1;
x(15)=0.5;y(15)=1;z(15)=1;
x(16)=0;y(16)=0.5;z(16)=1;
x(17)=0;y(17)=0;z(17)=0.5;
x(18)=1;y(18)=0;z(18)=0.5;
x(19)=1;y(19)=1;z(19)=0.5;
x(20)=0;y(20)=1;z(20)=0.5;
% % Next nodes are for 27 node brick (nodes on the face center)
% x(21)=0;y(21)=0.5;z(21)=0.5;
% x(22)=1;y(22)=0.5;z(22)=0.5;
% x(23)=0.5;y(23)=1;z(23)=0.5;
% x(24)=0.5;y(24)=0;z(24)=0.5;
% x(25)=0.5;y(25)=0.5;z(25)=0;
% x(26)=0.5;y(26)=0.5;z(26)=1;
% x(27)=0.5;y(27)=0.5;z(27)=0.5;
% % Nodal coordinates (10 node tet)
% x(1)=0;y(1)=0;z(1)=0;
% x(2)=1;y(2)=0;z(2)=0;
% x(3)=0;y(3)=1;z(3)=0;
% x(4)=0;y(4)=0;z(4)=1;
% x(5)=0.5;y(5)=0;z(5)=0;
% x(6)=0;y(6)=0.5;z(6)=0;
% x(7)=0;y(7)=0;z(7)=0.5;
% x(8)=0.5;y(8)=0.5;z(8)=0;
% x(9)=0;y(9)=0.5;z(9)=0.5;
% x(10)=0.5;y(10)=0.5;z(10)=0;

% Plot the element using shape functions
figure
hold on
plot3(x,y,z,'bo')

% Loop for checking the nodal coordinates (27 node brick)
for q1=0:0.5:1
    for q2=0:0.5:1
        for q3=0:0.5:1

%            % necessary for tet element 
%            if q1+q2+q3<=1
            
                for n=1:meltyp
                    % Value of the shape function
                    sf=sfn3d(q1,q2,q3,n,meltyp);
                    N(n)=sf;
                end

                % Multiply the shape function values to find the coordinates of the
                % position corresponding to the inputted values of q1, q2, q3.
                x1=N*x';
                y1=N*y';
                z1=N*z';

                % Plot the coordinates for checking
                plot3(x1,y1,z1,'r*');
                pause
%             end
            
        end
    end
end
