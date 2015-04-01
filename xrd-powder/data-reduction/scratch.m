%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotation for GE1 is 152.5
% XRDIMAGE.Image.RotAngle     = 152.5;
% XRDIMAGE.Image.RotAngle     = 67.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         imgi    = imrotate(imgi, XRDIMAGE.Image.RotAngle);
%         imgi    = [nan(1040, size(imgi,2)); imgi];
%         imgi    = [nan(1040, size(imgi,2)); imgi];
%         imgi    = [nan(size(imgi,1), 520) imgi nan(size(imgi,1), 520)];

%         imgi    = [imgi, zeros(size(imgi,1), 1040)];
%         imgi    = [zeros(520, size(imgi,2)); imgi; zeros(520, size(imgi,2))];

switch XRDIMAGE.Image.IsHydra
            case 0
                disp('single panel')
                [Lx, Ly]    = size(imgi);

                % imgi(1024:end,:)            = 2000;
                % imgi(:,1024:end)            = 1000;
                % imgi(1024:end,1024:end)     = 500;
                % imgi(1:1024,1:1024)         = 3000;
            case 1
                disp('GE1')
                imgi        = imrotate(imgi, 62.5, 'bilinear');
                
                [Lx, Ly]    = size(imgi);
                imgi    = [imgi; zeros(600, Ly)];
                [Lx, Ly]    = size(imgi);
                imgi    = [imgi zeros(3365, 600)];
                
                [Lx, Ly]    = size(imgi);
                XRDIMAGE.Instr.numpixels    = size(imgi,1);

                % imgi(1324:end,:)            = 2000;
                % imgi(:,1324:end)            = 1000;
                % imgi(1325:end,1325:end)     = 500;
                % imgi(1:1324,1:1324)         = 3000;
        
            case 2
                disp('GE2')
            case 3
                disp('GE3')
            case 4
                disp('GE4')
            otherwise
                disp('Unknown configuration!!')
        end
