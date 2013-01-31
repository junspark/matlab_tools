function funcOut=peakOnlyGeomModel(x)
instr       = evalin('base','instr');
peakInfo    = evalin('base','peakInfo');
material    = evalin('base','material');
cakeParms   = evalin('base','cakeParms');

centTrans=[x(1)    x(2)];
gammaYprime=x(3);
gammaXhatPrime=x(4);
D=x(5)*instr.distance;
if instr.DistType==0
    DistPars=[];
elseif instr.DistType==1
    a0  = x(6);
    DistPars=a0;
elseif instr.DistType==2
    a1  = x(6);
    a2  = x(7);
    n1  = x(8);
    n2  = x(9);
    DistPars=[a1 a2 n1 n2 instr.t1 instr.e1 instr.detectorSize];
elseif instr.DistType==3
    a1  = x(6);
    a2  = x(7);
    n1  = x(8);
    n2  = x(9);
    t1  = x(10);
    e1  = x(11);
    DistPars=[a1 a2 n1 n2 t1 e1 instr.detectorSize];
end

mappedData = XFormRadialSpectraPeaks(centTrans./1000,D, ...
    gammaYprime/1000, gammaXhatPrime/1000, ...
    peakInfo.loc'.*instr.pixelSize, ...
    cakeParms.azim, instr.tthRng,...
    DistPars,instr.DistType)';
tth_theo	= repmat(material.tth{1}(material.modelPeaks),size(mappedData,1),1);
% whos mappedData tth_theo
funcOut     = tth_theo - mappedData;