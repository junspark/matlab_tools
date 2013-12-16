function [material,HKLS]=findModelPeaks(instr,material,polImg)
lambda =instr.wavelength;
TTH=[];
for ii=1:length(material.latParm)
    if material.structure(ii)==1
        if ispc
            hkls=load('bcc.hkls');
        else
            hkls=load('bcc.hkls');
        end
        if ii==1
            HKLS=hkls;
        end
        xtalStruc='cubic';
    elseif material.structure(ii)==2
        if ispc
            hkls=load('fcc.hkls');
        else
            hkls=load('fcc.hkls');
        end
        if ii==1
            HKLS=hkls;
        end
        xtalStruc='cubic';
    elseif material.structure(ii)==3
        if ispc
            hkls=load('hcp.hkls');
        else
            hkls=load('hcp.hkls');
        end
        if ii==1
            HKLS=hkls;
        end
        xtalStruc='hexagonal';
    else
        keyboard
    end
    
    
    if material.structure(ii)==3
        tth=calcPlaneSpacing(...
            xtalStruc,...
            material.latParm{ii},... % atish, it was material.latParm([ii ii+1])
            1,...
            hkls,...
            lambda,...
            instr.tthRng(1),...
            instr.tthRng(2)...
            );
    else
        tth=calcPlaneSpacing(...
            xtalStruc,...
            material.latParm{ii},...
            1,...
            hkls,...
            lambda,...
            instr.tthRng(1),...
            instr.tthRng(2)...
            );
    end
    TTH=[TTH; tth];
end

material.tth    = TTH;

figure(1);
clf;
r=instr.distance.*tand(material.tth{1});
plot(r/instr.pixelSize,-20*ones(size(r)),'bs');
hold on;
plot(material.rList(:,1)/instr.pixelSize,-20*ones(size(material.rList,1),1),'r<');
plot(material.rList(:,2)/instr.pixelSize,-20*ones(size(material.rList,1),1),'r>');
r=instr.distance.*tand(material.tth{1});
plot(r(material.modelPeaks)/instr.pixelSize,-100*ones(size(material.modelPeaks)),'gs');

if material.num>1
    for nnn = 2:1:material.num
        r   = instr.distance.*tand(material.tth{nnn});
        plot(r/instr.pixelSize,-20*ones(size(r)),'rs');
    end
    legend('CeO2','Bounds','Model Peaks','Mat');
else
    legend('CeO2','LB','UB','Model Peaks');
end
plot(polImg.radius',polImg.intensity','-.');
ylabel('Intensity','FontSize',14,'FontWeight','Bold');
xlabel('Radial Distance (pixels)','FontSize',14,'FontWeight','Bold');
A=gca;
set(A,'FontSize',14,'FontWeight','Bold');