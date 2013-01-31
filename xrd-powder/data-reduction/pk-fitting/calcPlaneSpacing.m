function [tth]=calcPlaneSpacing(typeLat,latParams,numPhases,hkls,lambda,TTHMIN,TTHMAX)

ii  = 1;
for j = 1:numPhases
    if strcmp(typeLat, 'cubic')
        % get dhkl and tthhkl
        [dsp{j}, tmpTTH]    = PlaneSpacings(...
            latParams, ...
            typeLat, ...
            hkls', ...
            lambda);
        tmpTTH  = tmpTTH*2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identify dhkls that are the same with different hkl
        % example {511} and {333}
        % can be a separate function to kickout Duplicates
        % why CALIBFLAG needs to be tested?
        kickOut{j} = true(size(tmpTTH));
        % if ~isempty(CALIBFLAG)
        % if i == CALIBFLAG
        eqvDs   = DuplicateVectors(dsp{j});
        if ~isempty(eqvDs)
            tlen = length(eqvDs);
            for k = 1:tlen
                kickOut{j}(eqvDs{k}(2:end)) = 0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % marker for tmpTTH to keep in
        keepIn{j} = tmpTTH >= TTHMIN & tmpTTH <= TTHMAX;
        
        %         % xcutChk for checking XCUT with tmpTTH
        %         xcutChk{j}  = logical(ones(size(tmpTTH)));
        %         if ~isempty(XCUT)
        %             for k = 1:size(XCUT, 2)
        %                 xcutChk{j}  = ((tmpTTH <= XCUT(1, k) | tmpTTH >= XCUT(2, k))) & xcutChk{j};
        %             end
        %         end
        
        keepIn{j} = keepIn{j} & kickOut{j};% & xcutChk{j};
        
        % original TTHs
        fullTTH{j} = tmpTTH;
        
        % TTHs after filtering out duplicates, inside XBND, and outside XCUT
        tth{j}     = tmpTTH(keepIn{j});
        %latParams(j) = latParams(ii);
        %ii = ii + 1;
    elseif strcmp(typeLat, 'hexagonal')
        [dsp{j}, tmpTTH]    = PlaneSpacings(...
            latParams, ...
            typeLat, ...
            hkls', ...
            lambda);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identify dhkls that are the same with different hkl
        % example {511} and {333}
        % can be a separate function to kickout Duplicates
        % why CALIBFLAG needs to be tested?
        kickOut{j} = true(size(tmpTTH));
        % if ~isempty(CALIBFLAG)
        % if i == CALIBFLAG
        eqvDs   = DuplicateVectors(dsp{j});
        if ~isempty(eqvDs)
            tlen = length(eqvDs);
            for k = 1:tlen
                kickOut{j}(eqvDs{k}(2:end)) = 0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % marker for tmpTTH to keep in
        keepIn{j} = tmpTTH >= TTHMIN & tmpTTH <= TTHMAX;
        
        %         % xcutChk for checking XCUT with tmpTTH
        %         xcutChk{j}  = logical(ones(size(tmpTTH)));
        %         if ~isempty(XCUT)
        %             for k = 1:size(XCUT, 2)
        %                 xcutChk{j}  = ((tmpTTH <= XCUT(1, k) | tmpTTH >= XCUT(2, k))) & xcutChk{j};
        %             end
        %         end
        
        keepIn{j} = keepIn{j} & kickOut{j};% & xcutChk{j};
        
        % original TTHs
        fullTTH{j} = tmpTTH;
        
        % TTHs after filtering out duplicates, inside XBND, and outside XCUT
        tth{j}     = tmpTTH(keepIn{j});
        %latParams(j) = latParams(ii);
        %ii = ii + 1;

    end
end