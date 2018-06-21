if ~strcmp(filetype,'cor2'); irow=intersect(find(strcmp(meta.det1_fname,fileroot)),find(meta.det1_fnum==ifile)); irow=irow(length(irow));%  find data row corresponding to current file number AND fileroot
                else irow=intersect(find(meta.det1_fnum==ifile),find(meta.det1_fnum==iframe)); %for cor2, match both file number and frame number - not implemented in current metadata!
                end
                i0(index1)=(meta.scaler4_val(irow)-i0_dark*meta.integ_time(irow)).*meta.scaler4_units(irow).*ic_cf; %i0 integrated flux in photons;f25 to compensate for fact this before the 50um horiz slit
                i1(index1)=(meta.scaler6_val(irow)-i1_dark*meta.integ_time(irow)).*meta.scaler6_units(irow).*ic_cf; %i1 is NOT in photons, just proportional to i0 and trans
                samX(index1)=meta.samX2(irow); samY(index1)=meta.samY(irow); samZ(index1)=meta.samZ(irow); phi(index1)=meta.aY(irow);samEta(index1)=meta.aZ(irow);
                t_exp(index1)=meta.det1_time_per_frame(irow).*meta.det1_frames_per_file(irow); t_scaler(index1)=meta.integ_time(irow);
                crosshead_mm(index1)=meta.encoder4(irow);load_N(index1)=meta.ev2(irow);
                temperature1(index1)=meta.ev3(irow);temperature2(index1)=meta.ev4(irow);temperature3(index1)=meta.ev5(irow);
                trans(index1)=(meta.scaler6_val(irow)-i0_dark*meta.integ_time(irow))./(meta.scaler4_val(irow)-i1_dark*meta.integ_time(irow))./(i1_over_i0_air); %calculated transmission, normalized using no-sample ion chamber readings


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
