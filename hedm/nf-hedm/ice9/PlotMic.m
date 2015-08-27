function PlotMic( varargin )
%
% PlotMic is an updated version of the legacy function plot_mic
% This is D. B. Menasche build on a S.F. Li Original
% 11/14/2014
% dmenasche@gmail.com 
%
% Usage: PlotMic( Mode , Keyword , 'keyvalue' , Keyword2 , 'keyvalue2' ...
% -> for SingleLayer mode at a minimium, you need to include 'mic', 'sw' 
% and 'plotType' keywords
% -> for Volume mode, you need a 2 x n cell with the mics in row 1 and the
% sidewidths in row 2
%
% Valid Modes:
%   'Legacy' - use old syntax after listing this string
%   'SingleLayer'- use keywords for this and modes below
%   'Volume'- beware of plotting too big an area
%   'Grains'- plot select grains from the volume
%   'Facets'- plot facet voxels, but you need the facet file
% 
%
% Valid Keywords:
%   'mic'                       - mic file single layer
%   'micvolume'                 - a cell with all mic files in volume the 
%                                   sidewidths should be in the second
%                                   row
%   'sw'                        - sidewidth of mic file
%   'layerHeight'               - vertical (z) height of layer plotted
%   'alpha'                     - controls opacity: defaults to 1
%   'plotType' 
%       2 - confidence map
%       3 - rodrigues vectors
%       9 - plot grain IDs from column 11
%   'colorMap'                  - colormap for confidence style plots
%       'Red' - original red-blue
%       'Jet' - MATLAB jet
%       'Hot' - MATLAB hot
%       'Bone'- MATLAB bone
%   'rodColorMap'               - colormap for rodrigues vectors
%       'unscaledCFZ' - unscaled scaling of the cubic fundemental zone
%       'periodic'    - orientations near FZ edges are close in color
%       'scaledCFZ'   - rescale by minimum of ith rodrigues component !! not
%       yet implemetned!!!
%   'ApplyPrerotationTranslation' translation mic by a vector [x y]
%   'ApplyPostRotationTranslation' translate mic by a vector [x y]
%   'TransformOrientations'     - transform orientations by a quaternion, q
%   'Rotation'                  - give angle in degrees to rotate
%   '3D'                        - plot in 3D ( 0 or 1 ) 
%   'TopFace'                   - top face to 3D plot ( 0 or 1, enabled by default)
%   'BottomFace'                - bottom face on 3D ( 0 or 1, enables by default)
%   'Edges'                     - add edges for 3D plotting ( 0 or 1 )
%   'Boundaries'                - if you have boundaryvoxels marked in col
%                                   17, you plot only boundary 
    hold on;

    %%% Defaults:
    
    sOptions.mode = 'Legacy';
    
    if nargin == 5
        sOptions.mic = varargin{2};
        sOptions.sideWidth = varargin{3};
        sOptions.plotType = varargin{4};
        sOptions.confidenceRange = varargin{5};       
    end    
    
    sOptions.bHasRadius = 0;
    sOptions.radius = [];
    sOptions.bBoundaries = 1;
    sOptions.fRotation = 0;
    sOptions.confidenceRange = 0.000;
    sOptions.colorMap = 'Red';
    sOptions.layerHeight = 0.000;
    sOptions.bEdges = 0;
    sOptions.bTopFace = 1;
    sOptions.bBottomFace = 1;
    sOptions.alpha = 1;
    sOptions.sRGB = 'unscaledCFZ'; % can be also 'periodic'  
    sOptions.voxelColors3D = [];
    sOptions.voxelColors2D = [];
    sOptions.b3D = 0;
    sOptions.bTransformOrientations = 0;
    sOptions.bApplyPrerotationTranslation = 0;
    sOptions.bApplyPostrotationTranslation = 0;
    sOptions.vPreTranslation = [0 0 0];
    sOptions.vGrainsToPlot = [];
    sOptions.minLabelSize = 100;
    sOptions.bTriangulation = 0;
    sOptions.layerSpacing = 0.004;
    sOptions.KAMtopScale = 2.5;
    sOptions.fConfThresh = 0.3;
    sOptions.fEdgeAlpha = 0.5;
    sOptions.bPlotTripleLines = 0;
    sOptions.sTripleLineCmap = 'Hot';
    sOptions.bTripleLineSingleColor = 0;
    sOptions.sTLColor = 'w';
    sOptions.plotType = 3;
    
    while ~isempty(varargin)
        
        switch upper(varargin{1})
            case 'LEGACY'
                plot_mic( sOptions.mic , sOptions.sideWidth , sOptions.plotType , sOptions.confidenceRange , sOptions.sRGB );
                break
            case 'SINGLELAYER'
                sOptions.mode = 'SingleLayer';             
                varargin(1) = [];
            case 'VOLUME'
                sOptions.mode = 'Volume';
                varargin(1) = [];
            case 'GRAINS'
                sOptions.mode = 'Grains';                                      
                varargin(1) = [];
            case 'FACETS'
                sOptions.mode = 'Facets';                                      
                varargin(1) = [];
            case 'LAYERHEIGHT'
                sOptions.layerHeight = varargin{2};
                varargin(1:2) = [];
            case 'ALPHA'
                if varargin{2} > 1 || varargin{2} < 0
                    ME = MException('PlotMic:AlphaError','Alpha field must be in range [0,1]');
                    throw(ME)
                end
                sOptions.alpha = varargin{2};
                varargin(1:2) = [];                                                
            case 'MIC'
                sOptions.mic = varargin{2};
                varargin(1:2) = [];
            case 'BOUNDARYPOINTS'
                sOptions.facets = varargin{2};
                varargin(1:2) = [];
            case 'SW' 
                sOptions.sideWidth = varargin{2};
                varargin(1:2) = [];
            case 'PLOTTYPE' 
                sOptions.plotType = varargin{2};
                varargin(1:2) = [];
            case 'COLORMAP'
                expectedMaps = {'Red','Hot','Jet','Bone'};
                mapName = validatestring(varargin{2},expectedMaps,'PlotMic','colorMap',1);
                sOptions.colorMap = varargin{2};                
                varargin(1:2) = [];
            case 'RODCOLORMAP'
                expectedMaps = {'unscaledCFZ','periodic','scaledCFZ'};
                mapName = validatestring(varargin{2},expectedMaps,'PlotMic','rodColorMap',1);
                sOptions.sRGB = varargin{2};
                varargin(1:2) = [];
            case 'ROTATION'
                sOptions.fRotation = varargin{2};
                varargin(1:2) = [];
            case 'EDGES'
                sOptions.bEdges = varargin{2};
                varargin(1:2) = [];
            case 'EDGEALPHA'
                sOptions.fEdgeAlpha = varargin{2};
                varargin(1:2) = [];
            case 'TOPFACE'
                sOptions.bTopFace = varargin{2};
                varargin(1:2) = [];
            case 'BOTTOMFACE' 
                sOptions.bBottomFace = varargin{2};
                varargin(1:2) = [];
            case '3D'
                sOptions.b3D = varargin{2};
                varargin(1:2) = [];
            case 'TRANSFORMORIENTATIONS'                
                sOptions.bTransformOrientations = 1;
                sOptions.qTransform = varargin{2};
                varargin(1:2) = [];          
            case 'APPLYPREROTATIONTRANSLATION'
                sOptions.bApplyPrerotationTranslation = 1;
                sOptions.vPreTranslation = varargin{2};
                varargin(1:2) = []; 
            case 'APPLYPOSTROTATIONTRANSLATION'
                sOptions.bApplyPostrotationTranslation = 1;
                sOptions.vPostTranslation = varargin{2};
                varargin(1:2) = [];
            case 'BOUNDARIES'
                sOptions.bBoundaries = varargin{2};
                varargin(1:2) = [];                 
            case 'GRAINSTOPLOT'
                sOptions.vGrainsToPlot = varargin{2};
                varargin(1:2) = []; 
            case 'MINLABELSIZE'
                sOptions.minLabelSize = varargin{2};
                varargin(1:2) = []; 
            case 'MICVOLUME'
                sOptions.volume = varargin{2};
                varargin(1:2) = []; 
            case 'RADIUS'
                sOptions.bHasRadius = 1;
                sOptions.radius = varargin{2};
                varargin(1:2) = []; 
            case 'TRIPLELINES'
                sOptions.bPlotTripleLines = 1;
                sOptions.TripleLines = varargin{2};
                varargin(1:2) = []; 
            case 'TRIPLELINECOLORMAP' 
                expectedMaps = {'Bone','Jet','Hot'};
                mapName = validatestring(varargin{2},expectedMaps,'PlotMic','tripleLineColorMap',1);
                sOptions.sTripleLineCmap = varargin{2};
                varargin(1:2) = []; 
            case 'TRIPLELINESINGLECOLOR'
                sOptions.bTripleLineSingleColor = 1;
                sOptions.sTLColor = varargin{2};
                varargin(1:2) = [];   
            case 'TRIANGULATION'
                sOptions.bTriangulation = varargin{2};
                varargin(1:2) = []; 
            case 'LAYERSPACING'
                sOptions.layerSpacing = varargin{2};
                varargin(1:2) = [];
            case 'KAMTOPSCALE'
                sOptions.KAMtopScale = varargin{2};
                varargin(1:2) = [];           
            case 'COLORSCALARS'
                sOptions.colorScalars = varargin{2};
                varargin(1:2) = [];
            otherwise
                 ME = MException('PlotMic:UnrecognizedOption','Check your keyword inputs.');
                 throw(ME)
        end        
    end    
        
    if strcmp( sOptions.mode , 'SingleLayer' ) == 1
        
        if sOptions.bHasRadius 
            ind = GetMicIndicesInRadius( sOptions.mic , sOptions.radius , 0 , 0 );
            sOptions.mic = sOptions.mic(ind,:);
        end        
        
        if sOptions.bTransformOrientations 
            disp('REDUCING TO THE CUBIC FUNDEMENTAL ZONE')
            qs = QuatOfRMat( RMatOfBunge( sOptions.mic(:,7:9)' ,'degrees' ) );
            q0 = repmat( sOptions.qTransform , 1 , size( sOptions.mic , 1 ) );
            sOptions.mic(:,7:9) = BungeOfRMat( RMatOfQuat( ToFundamentalRegionQ( QuatProd( q0 , qs ) , CubSymmetries ) ) , 'degrees' )';
        end
        if sOptions.bApplyPrerotationTranslation
            sOptions.mic(:,1) = sOptions.mic(:,1) + sOptions.vPreTranslation(1);
            sOptions.mic(:,2) = sOptions.mic(:,2) + sOptions.vPreTranslation(2);
            %sOptions.mic(:,3) = sOptions.mic(:,3) +
            %sOptions.vPreTranslation(3); z not yet implemented.
        end
        PlotSingleLayer( sOptions )    
    end
    
    if strcmp( sOptions.mode , 'Grains' ) == 1
        PlotGrainsFromVolume( sOptions )   
        if sOptions.bPlotTripleLines
            PlotTripleLines( sOptions )
        end
    end
    
    if strcmp( sOptions.mode , 'Volume' ) == 1

        if sOptions.bHasRadius
            for i = 1:size( sOptions.volume , 2 ) 
                mic = sOptions.volume{1,i};
                ind = GetMicIndicesInRadius( mic , sOptions.radius , 0 , 0 );                
                sOptions.volume{1,i} = mic(ind,:);
            end
        end     
        if sOptions.bTransformOrientations 
            disp('REDUCING TO THE CUBIC FUNDEMENTAL ZONE')
            qs = QuatOfRMat( RMatOfBunge( sOptions.mic(:,7:9)' ,'degrees' ) );
            q0 = repmat( sOptions.qTransform , 1 , size( sOptions.mic , 1 ) );
            sOptions.mic(:,7:9) = BungeOfRMat( RMatOfQuat( ToFundamentalRegionQ( QuatProd( q0 , qs ) , CubSymmetries ) ) , 'degrees' )';
        end
        if sOptions.bApplyPrerotationTranslation
            sOptions.mic(:,1) = sOptions.mic(:,1) + sOptions.vPreTranslation(1);
            sOptions.mic(:,2) = sOptions.mic(:,2) + sOptions.vPreTranslation(2);
            %sOptions.mic(:,3) = sOptions.mic(:,3) +
            %sOptions.vPreTranslation(3); z not yet implemented.
        end
        PlotVolume( sOptions )    
    end
    
    if strcmp( sOptions.mode , 'Facets' ) == 1
        PlotGrainFacetsFromPoints(  sOptions );
        if sOptions.bPlotTripleLines
            PlotTripleLines( sOptions )
        end
    end        
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      New 'SingleLayer' wrapper on plot_mic - features added.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotSingleLayer( varargin )
    
    sOptions = varargin{1};
    snp = sOptions.mic;
    sidewidth = sOptions.sideWidth;
    plotType = sOptions.plotType;
    ConfidenceRange = sOptions.confidenceRange;
    bEdges = sOptions.bEdges;
        
    [tri_x, tri_y, snp_Euler, snp_Conf , snp_Euler_Unordered ] = GetTriangleVertices( sOptions );
    
    if (plotType == 2 || plotType == 9 ) % plot confidence map
        if strcmp( sOptions.colorMap , 'Red' ) 
            tmp = [ snp_Conf ];
            tmp = [ snp_Conf ] - min( snp_Conf );   % original auto rescale
            tmp = tmp /( max(snp_Conf) - min(snp_Conf) );  % original auto rescale
            tmp = [tmp, zeros(length(tmp), 1), 1-tmp];
            sOptions.voxelColors2D = tmp;
            xConf = [0:0.01:1]';
            colormap( [ xConf, xConf * 0, 1- xConf]  );
            caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
            colorbar( 'location', 'eastoutside');
            sOptions.voxelColors2D = tmp;
        else
            if strcmp( upper(sOptions.colorMap) , 'HOT' ) 
                j = hot;
            elseif strcmp( upper(sOptions.colorMap) , 'BONE' ) 
                j = bone;
            elseif strcmp( upper(sOptions.colorMap) , 'JET' )
                j = jet;
            end
            tmp = [ snp_Conf ];
            tmp = tmp - min( snp_Conf );
            tmp = tmp /( max(snp_Conf) - min(snp_Conf) );  % original auto rescale            
            %j = DivideColormapNTimes( j , 4) ;
            m = size(j,1);
            tmp = j(floor((tmp*(m-1)))+1,:); 
            sOptions.voxelColors2D = tmp;
            colormap( j );
            caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
            colorbar
        end
                                        
    elseif(plotType==3) % plot rodrigues vectors            
                
        sOptions.voxelColors2D = RodOfQuat( QuatOfRMat( RMatOfBunge(snp_Euler', 'degrees') ) )';
        
        if sOptions.b3D
            sOptions.voxelColors3D = RodOfQuat( QuatOfRMat( RMatOfBunge(snp_Euler_Unordered', 'degrees') ) )';
        end        
        
        if strcmp(sOptions.sRGB,'periodic') == 1
            tmp = sOptions.voxelColors2D;
            x = 0.8284/2;   
            RGB = [abs(sin(pi*tmp(:,1)/x)) abs(cos(pi*tmp(:,2)/x)) abs(sin(pi*tmp(:,3)/x))];
            sOptions.voxelColors2D = RGB;
            if ~isempty( sOptions.voxelColors3D )
                tmp = sOptions.voxelColors3D;
                x = 0.8284/2;   
                RGB = [abs(sin(pi*tmp(:,1)/x)) abs(cos(pi*tmp(:,2)/x)) abs(sin(pi*tmp(:,3)/x))];
                sOptions.voxelColors3D = RGB;
            end
        elseif strcmp(sOptions.sRGB,'unscaledCFZ') == 1

            tmp = sOptions.voxelColors2D;
            x = 0.8284/2;   
            RGB = [tmp(:,1)/(2*x)+0.5 tmp(:,2)/(2*x)+0.5 tmp(:,3)/(2*x)+0.5];
            RGB( RGB(:,1) < 0 ,1) = 0;
            RGB( RGB(:,2) < 0 ,2) = 0;
            RGB( RGB(:,3) < 0 ,3) = 0;
            sOptions.voxelColors2D = RGB;
            
            if ~isempty( sOptions.voxelColors3D )
                tmp = sOptions.voxelColors3D;
                x = 0.8284/2;   
                RGB = [tmp(:,1)/(2*x)+0.5 tmp(:,2)/(2*x)+0.5 tmp(:,3)/(2*x)+0.5];
                RGB( RGB(:,1) < 0 ,1) = 0;
                RGB( RGB(:,2) < 0 ,2) = 0;
                RGB( RGB(:,3) < 0 ,3) = 0;
                sOptions.voxelColors3D = RGB;
            end
            
        end
                        
    elseif(plotType == 4 || plotType == 5 ) % plot confidence map, no scaling
        if strcmp( sOptions.colorMap , 'Red' ) 
            tmp = [ snp_Conf ];
            tmp = [tmp, zeros(length(tmp), 1), 1-tmp];
            sOptions.voxelColors2D = tmp;
            xConf = [0:0.01:1]';
            %colormap( [ xConf, xConf * 0, 1- xConf]  );
            %caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
            %colorbar( 'location', 'eastoutside');
        else
            if strcmp( upper(sOptions.colorMap) , 'HOT' ) 
                j = hot;
            elseif strcmp( upper(sOptions.colorMap) , 'BONE' ) 
                j = bone;
            elseif strcmp( upper(sOptions.colorMap) , 'JET' )
                j = jet;
            end
            
            tmp = [ snp_Conf ];
            
%            if length( tmp == 0 ) == length( tmp )
%                 ME = MException('PlotMic:UndefinedQuantity','Undefined quanity - You havent yet computed what youre trying to plot.');
%                 throw(ME)
%            end
            %j = DivideColormapNTimes( j , 4) ;
            m = size(j,1);
            if plotType==4
                tmp = j(floor((tmp/max(tmp)*(m-1)))+1,:); 
            elseif plotType==5
                % hardcode 4 degrees
                tmp( tmp > sOptions.KAMtopScale ) = sOptions.KAMtopScale;
                tmp = j(floor(((tmp/sOptions.KAMtopScale)*(m-1)))+1,:); 
            end

            sOptions.voxelColors2D = tmp;
            colormap( j );
            %caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
            %colorbar
        end
        
    elseif (plotType == 6)
        if strcmp( upper(sOptions.colorMap) , 'HOT' ) 
            j = hot;
        elseif strcmp( upper(sOptions.colorMap) , 'BONE' ) 
            j = bone;
        elseif strcmp( upper(sOptions.colorMap) , 'JET' )
            j = jet;
        else
            j = jet;
        end
        vScalars = [ snp_Conf ];;
        %fMin = min( vScalars );
        fMax = max( vScalars );
        vScalars = floor((vScalars/fMax)*63)+1;
        for i = 1:length( vScalars )
            vColors(i,:) = j(vScalars(i),:);
        end
        sOptions.voxelColors2D = vColors;
                
    end
  
    if strcmp( upper(sOptions.mode) , 'SINGLELAYER')

        if plotType == 4 || plotType == 2
            colorbar( 'location', 'eastoutside');
        elseif plotType == 5
            colorbar( 'location', 'eastoutside');
            caxis([0 sOptions.KAMtopScale])
        elseif plotType == 6
            colorbar( 'location', 'eastoutside');
            caxis([0 fMax])
        end         
    end
        
    if ~(sOptions.b3D)    
        tmp = sOptions.voxelColors2D;
        %size(tmp)
        tri_color = zeros(1, size(tmp,1), size(tmp,2));
        tri_color(1, :, :) = tmp;    
        tri_z = zeros(3,size(tri_x,2)) + sOptions.layerHeight;
        h = patch(tri_x, tri_y, tri_z, tri_color, 'EdgeColor', 'none');  

        alpha(h,sOptions.alpha);

        if plotType == 9
            nMinSize = sOptions.minLabelSize;
            gids = unique(snp(:,11));
            for i=1:length(gids)
                ind = find( snp_Euler_Unordered(:,1) == gids(i) );                
                if length(ind) > nMinSize   
                    x_avg = mean( snp_Euler_Unordered(ind,2)  );                    
                    y_avg = mean( snp_Euler_Unordered(ind,3)  );
                    text(x_avg,y_avg, num2str(gids(i)), 'Color', 'w','FontSize',12);
                end
            end
        end
        
    elseif sOptions.b3D 
       Plot3DVoxels( sOptions )         
       if sOptions.bTopFace
            tmp = sOptions.voxelColors2D;
            size_vec = size(tmp);
            tri_color = zeros(1, size_vec(1), size_vec(2));
            tri_color(1, :, :) = tmp;    
            tri_z = zeros(3,size(tri_x,2)) + sOptions.layerHeight+0.002;
            h = patch(tri_x, tri_y, tri_z, tri_color, 'EdgeColor', 'none');
            alpha(h,sOptions.alpha);
       end
       if sOptions.bBottomFace 
            tmp = sOptions.voxelColors2D;
            size_vec = size(tmp);
            tri_color = zeros(1, size_vec(1), size_vec(2));
            tri_color(1, :, :) = tmp;    
            tri_z = zeros(3,size(tri_x,2)) + sOptions.layerHeight+0.002;
            h = patch(tri_x, tri_y, tri_z, tri_color, 'EdgeColor', 'none');
            alpha(h,sOptions.alpha);
       end
    end
    box on;
    axis square;
    axis equal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      PlotVolume
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotVolume( varargin )
    % IF YOU HAVE NO BOUNDARIES AN YOU ASK FOR BOUNDARIES THIS FIALS
    sOptions = varargin{1};
    DMic = sOptions.volume;
    for i = 1:size( DMic , 2 )
        if sOptions.bBoundaries
            ind = find( DMic{1,i}(:,17) == 1 );
        else
            ind = 1:size( DMic{1,i} , 1 );
        end
        sOptions.mic = DMic{1,i}(ind,:); 
        sOptions.sideWidth = DMic{2,i};
        sOptions.layerHeight = DMic{1,i}(1,3);
        PlotSingleLayer( sOptions )
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      PlotGrainsFromVolume
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotGrainsFromVolume( varargin )

    sOptions = varargin{1};
    DMic = sOptions.volume;        
    
    if ~sOptions.bTriangulation && ~strcmp( sOptions.sRGB , 'scaledCFZ' )
    
        indCell = cell(1,size(DMic,2));
        for grains = sOptions.vGrainsToPlot

            for i = 1:size(DMic,2)        
                if sOptions.bBoundaries
                    ind = find( DMic{1,i}(:,11) == grains & DMic{1,i}(:,17) == 1);
                else
                    ind = find( DMic{1,i}(:,11) == grains );
                end
                indCell{i} = ind;
            end
            disp('REDUCING TO THE CUBIC FUNDEMENTAL ZONE')
            for i = 1:size(indCell,2)            
                if ~isempty(indCell{i})                                   
                    sOptions.mic = DMic{1,i}(indCell{i},:); 
                    sOptions.sideWidth = DMic{2,i};
                    sOptions.layerHeight = DMic{1,i}(1,3);
                    if sOptions.bTransformOrientations             
                        qs = QuatOfRMat( RMatOfBunge( sOptions.mic(:,7:9)' ,'degrees' ) );
                        q0 = repmat( sOptions.qTransform , 1 , size( sOptions.mic , 1 ) );
                        sOptions.mic(:,7:9) = BungeOfRMat( RMatOfQuat( ToFundamentalRegionQ( QuatProd( q0 , qs ) , CubSymmetries ) ) , 'degrees' )';
                    end
                    if sOptions.bApplyPrerotationTranslation
                        sOptions.mic(:,1) = sOptions.mic(:,1) + sOptions.vPreTranslation(1);
                        sOptions.mic(:,2) = sOptions.mic(:,2) + sOptions.vPreTranslation(2);                       
                    end
                    PlotSingleLayer( sOptions )
                end
            end

        end
        
    elseif strcmp( sOptions.sRGB , 'scaledCFZ' ) && ~sOptions.bTriangulation
        %%%%%%%%%%%%%%%%%%%%%%%% NOT IMPLEMENTED
        for grains = sOptions.vGrainsToPlot
            thisGrain = [];
            for i = 1:size(DMic,2)        
                if sOptions.bBoundaries
                    ind = find( DMic{1,i}(:,11) == grains & DMic{1,i}(:,17) == 1 & DMic{1,i}(:,4) == 1);                 
                else
                    ind = find( DMic{1,i}(:,11) == grains );
                end
                if ~isempty( ind )
                    thisGrain = [thisGrain; DMic{1,i}(ind,:) ];      
                end
            end
            
            tmp = RodOfQuat( QuatOfRMat( RMatOfBunge( thisGrain(:,7:9) , 'degrees' ) ) )'; 
            r1 = min( tmp(:,1) );
            r2 = min( tmp(:,2) );
            r3 = min( tmp(:,3) );
            tmp(:,1) = tmp(:,1) - r1 ;
            tmp(:,2) = tmp(:,2) - r2 ;
            tmp(:,3) = tmp(:,3) - r3 ;                
            RGB = [tmp(:,1)/(2*x)+0.5 tmp(:,2)/(2*x)+0.5 tmp(:,3)/(2*x)+0.5];
            RGB( RGB(:,1) < 0 ,1) = 0;
            RGB( RGB(:,2) < 0 ,2) = 0;
            RGB( RGB(:,3) < 0 ,3) = 0;
        
        end
    else
         
        for grains = sOptions.vGrainsToPlot
            thisGrainsPoints = [];
            or = [0 0 0];
            bFlag = 1;
            for i = 1:size(DMic,2)        
                if sOptions.bBoundaries
                    ind = find( DMic{1,i}(:,11) == grains & DMic{1,i}(:,17) == 1 & DMic{1,i}(:,4) == 1);
                    if ~isempty( ind ) && bFlag
                        or = DMic{1,i}(ind(1),7:9);
                        bFlag = 0;
                    end                    
                else
                    ind = find( DMic{1,i}(:,11) == grains );
                end
                if ~isempty( ind )
                    thisGrainsPoints = [thisGrainsPoints; DMic{1,i}(ind,1) ...
                        DMic{1,i}(ind,2) DMic{1,i}(ind,3)];                    
                end
            end
            
            DT = DelaunayTri(thisGrainsPoints(:,1),thisGrainsPoints(:,2),thisGrainsPoints(:,3));                        
            
            if strcmp(sOptions.sRGB,'unscaledCFZ') == 1

                tmp = RodOfQuat( QuatOfRMat( RMatOfBunge( or' , 'degrees' ) ) )';                
                x = 0.8284/2;   
                RGB = [tmp(:,1)/(2*x)+0.5 tmp(:,2)/(2*x)+0.5 tmp(:,3)/(2*x)+0.5];
                RGB( RGB(:,1) < 0 ,1) = 0;
                RGB( RGB(:,2) < 0 ,2) = 0;
                RGB( RGB(:,3) < 0 ,3) = 0;
                faceColor = RGB;
            
            elseif strcmp(sOptions.sRGB,'periodic') == 1
                
                tmp = RodOfQuat( QuatOfRMat( RMatOfBunge( or' , 'degrees' ) ) )'; 
                x = 0.8284/2;   
                RGB = [abs(sin(pi*tmp(:,1)/x)) abs(cos(pi*tmp(:,2)/x)) abs(sin(pi*tmp(:,3)/x))];
                faceColor = RGB;                                       
                         
            end
            
            
            if sOptions.bEdges
                tetramesh(DT,'FaceColor',faceColor,'EdgeColor',[0 0 0],'FaceAlpha',sOptions.alpha)
            else
                tetramesh(DT,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',sOptions.alpha)
            end
                        
        end
    end
    if sOptions.plotType == 4 || sOptions.plotType == 2
            colorbar( 'location', 'eastoutside');
    elseif sOptions.plotType == 5
            colorbar( 'location', 'eastoutside');
            caxis([0 sOptions.KAMtopScale])
    end        
    axis equal;
    %axis square;
end

function PlotGrainFacetsFromPoints( varargin )

    sOptions = varargin{1};

    % need a sw defined. 

    sw = sOptions.sideWidth;
    mic = sOptions.mic;
    gen = mic(1,5);
    Facets = sOptions.facets;        
    
    if sOptions.bTriangulation 
        disp('Cannot do triangulation of a facet!')
        return
    end

    disp('Coloring by facets, ignoring sRGB option. Use ColorMap option.');
    
    for grains = sOptions.vGrainsToPlot

        ind = find( cell2mat(Facets(1,:)) == grains );
        
        if strcmp(sOptions.colorMap,'Red')
            sOptions.colorMap = 'Jet';
        end
        
        if strcmp(sOptions.colorMap,'Jet')
            cmap = jet;
        end
        
        if strcmp(sOptions.colorMap,'Hot')
            cmap = hot;
        end
        
        if strcmp(sOptions.colorMap,'Bone')
            cmap = bone;
        end
        
        nUniqueFaces = unique( Facets{2,ind}(:,5) );
        
        for i = 1:length( nUniqueFaces )
            ind2 = find( Facets{2,ind}(:,5) == nUniqueFaces(i) );
    
            cindx = random('uniform',1,62);       
            sOptions.voxelColors2D = cmap(floor(cindx),:);
       
            sOptions.XYZpoints = Facets{2,ind}(ind2,1:3);
            sOptions.Kvec = Facets{2,ind}(ind2,4);                        
            
            if sOptions.bApplyPrerotationTranslation
                    sOptions.XYZpoints(:,1) = sOptions.XYZpoints(:,1) + sOptions.vPreTranslation(1);
                    sOptions.XYZpoints(:,2) = sOptions.XYZpoints(:,2) + sOptions.vPreTranslation(2);                       
            end

            fAngle = sOptions.fRotation;
            if fAngle ~= 0 
                % not yet tested - actually shouldn't work properly.
                th = fAngle*pi/180;
                R = [cos(th) -sin(th) 0 ;sin(th) cos(th) 0; 0 0 1];
                sOptions.XYZpoints(:,1:2) = (R*sOptions.XYZpoints(:,1:2)')';
            end

            if sOptions.bApplyPostrotationTranslation
                sOptions.XYZpoints(:,1) = sOptions.XYZpoints(:,1)+sOptions.vPostTranslation(1);
                sOptions.XYZpoints(:,2) = sOptions.XYZpoints(:,2)+sOptions.vPostTranslation(2);
            end


            dS=sw*(1/2)^gen;    
            sOptions.XYZpoints(:,1) = sOptions.XYZpoints(:,1)-0.5*dS;
            % lost info, use only +oriented triangles
            sOptions.XYZpoints(:,2) = sOptions.XYZpoints(:,2)-0.5*0.5774*dS*(-1).^(sOptions.Kvec-1);

            PlotVoxelsFromPoints( sOptions )

        end
    end

end

function PlotTripleLines( varargin )

    sOptions = varargin{1};
    
    tl_organized = sOptions.TripleLines;
    
    if ~sOptions.bTripleLineSingleColor
        if strcmp( sOptions.sTripleLineCmap , 'Jet' )
            cmap = jet;            
        end
        if strcmp( sOptions.sTripleLineCmap , 'Hot' )
            cmap = hot;            
        end
        if strcmp( sOptions.sTripleLineCmap , 'Bone' )
            cmap = bone;            
        end        
    end
    
    
    for i = sOptions.vGrainsToPlot
        ind = find( cell2mat( tl_organized(1,:) ) == i );        
        if sOptions.bTripleLineSingleColor
            hold on;
            plot3( tl_organized{3,ind}(:,1) , tl_organized{3,ind}(:,2) , ...
                tl_organized{3,ind}(:,3) , 'LineStyle' , '.' , 'color' , ...
                sOptions.sTLColor )                
        else
            lineIds = unique( tl_organized{3,ind}(:,10) );
            for i = 1:size(lineIds,1) 
                 ind2 = find( tl_organized{3,ind}(:,10) == lineIds(i) );
                 cindx = random('uniform',1,62);
                 plot3( tl_organized{3,ind}(ind2,1) ,  tl_organized{3,ind}(ind2,2) , tl_organized{3,ind}(ind2,3) , 'LineStyle' , '.' , 'color' , cmap(floor(cindx),:),'MarkerSize',6.5 )
            end
        end                  
    end


    
    

end

function PlotVoxelsFromPoints( varargin )
    
    sOptions = varargin{1};
    sw = sOptions.sideWidth;
    sw = sw*(1/2).^(sOptions.mic(1,5));
    dZ = sOptions.layerSpacing;
    ColorVec = sOptions.voxelColors2D;
    bPlotEdges = sOptions.bEdges;
    opacity = sOptions.alpha;
    
    points = zeros(size(sOptions.XYZpoints,1),3);
    points(:,1:2) = sOptions.XYZpoints(:,1:2);
    points(:,3) = sOptions.XYZpoints(:,3)-0.002;
    points(:,4:6) = points(:,1:3)+repmat(sw*[1 0 0],size(sOptions.XYZpoints,1),1);
    sgnUp = (-1).^(sOptions.Kvec-1);
    insert = repmat(sw*[0.5 sin(pi/3) 0],size(sOptions.XYZpoints,1),1);
    
    insert(:,2) = insert(:,2).*sgnUp;
    points(:,7:9) = points(:,1:3)+insert;
    points(:,10:12) = points(:,1:3)+repmat([0 0 dZ],size(sOptions.XYZpoints,1),1);
    points(:,13:15) = points(:,4:6)+repmat([0 0 dZ],size(sOptions.XYZpoints,1),1);
    points(:,16:18) = points(:,7:9)+repmat([0 0 dZ],size(sOptions.XYZpoints,1),1); 

    verts = [points(:,1:3);points(:,4:6);points(:,7:9);points(:,10:12);...
            points(:,13:15); points(:,16:18)];
    
    p1_indx = 1;
    p2_indx = 1+size(points(:,1:3),1);
    p3_indx = 1+2*size(points(:,1:3),1);
    p4_indx = 1+3*size(points(:,1:3),1);
    p5_indx = 1+4*size(points(:,1:3),1);
    p6_indx = 1+5*size(points(:,1:3),1);
    
    faces = zeros(size(sOptions.XYZpoints,1)*3,4);
  
    for i = 1:3:size( sOptions.XYZpoints,1 )*3
        faces(i,:) = [p1_indx p4_indx p5_indx p2_indx];
        faces(i+1,:) = [p2_indx p5_indx p6_indx p3_indx];
        faces(i+2,:) = [p1_indx p4_indx p6_indx p3_indx];
        p1_indx = p1_indx+1;
        p2_indx = p2_indx+1;
        p3_indx = p3_indx+1;
        p4_indx = p4_indx+1;
        p5_indx = p5_indx+1;
        p6_indx = p6_indx+1;
    end
        
    if sOptions.bTopFace
        topFaces = zeros(size(sOptions.XYZpoints,1),3);
        p1_indx = 1;
        p2_indx = 1+size(points(:,1:3),1);
        p3_indx = 1+2*size(points(:,1:3),1);
        for i = 1:size( topFaces , 1 );
            topFaces(i,:) = [p1_indx p2_indx p3_indx];
            p1_indx = p1_indx+1;
            p2_indx = p2_indx+1;
            p3_indx = p3_indx+1;
        end
        topFaceColors = zeros(size(sOptions.XYZpoints,1),3);
        topFaceColors = repmat( ColorVec , size(sOptions.XYZpoints,1) ,1);
    end       
    
    if sOptions.bBottomFace
        botFaces = zeros(size(sOptions.XYZpoints,1),3);
        p4_indx = 1+3*size(points(:,1:3),1);
        p5_indx = 1+4*size(points(:,1:3),1);
        p6_indx = 1+5*size(points(:,1:3),1);
        for i = 1:size( botFaces , 1 );
            botFaces(i,:) = [p4_indx p5_indx p6_indx];
            p4_indx = p4_indx+1;
            p5_indx = p5_indx+1;
            p6_indx = p6_indx+1;
        end
        botFaceColors = zeros(size(sOptions.XYZpoints,1),3);
        botFaceColors = repmat( ColorVec , size(sOptions.XYZpoints,1) ,1);
    end  
    
    faceColors = zeros(size(sOptions.XYZpoints,1)*3,3);
    indx=1;
    for i = 1:3:size( faceColors , 1) 
        faceColors(i,:) = ColorVec;
        faceColors(i+1,:) = ColorVec;
        faceColors(i+2,:) = ColorVec;
        indx=indx+1;
    end       

    if bPlotEdges        
        s = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','FaceAlpha',opacity,'EdgeAlpha',sOptions.fEdgeAlpha);
        if sOptions.bTopFace            
            s = patch('Faces',topFaces,'Vertices',verts,'FaceVertexCData',topFaceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
        end
        if sOptions.bBottomFace            
            s = patch('Faces',botFaces,'Vertices',verts,'FaceVertexCData',botFaceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
        end
    else                
        s = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
        if sOptions.bTopFace            
            s = patch('Faces',topFaces,'Vertices',verts,'FaceVertexCData',topFaceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
        end
        if sOptions.bBottomFace            
            s = patch('Faces',botFaces,'Vertices',verts,'FaceVertexCData',botFaceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
        end
    end
    
    
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      3D Voxel plotting function ( PlotManyVoxels ) 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot3DVoxels( varargin )
    
    sOptions = varargin{1};
    mic = sOptions.mic;
    sw = sOptions.sideWidth;
    sw = sw*(1/2).^(mic(1,5));
    fZ = sOptions.layerHeight;
    bPlotEdges = sOptions.bEdges;
    % ColorVec has to match the number of faces you want to plot
    ColorVec = sOptions.voxelColors3D;
    opacity = sOptions.alpha;
    fAngle = sOptions.fRotation;
    dZ = sOptions.layerSpacing;
    
    points = zeros(size(mic,1),18);
    points(:,1:2) = mic(:,1:2);
    points(:,3) = fZ-0.002;
    points(:,4:6) = points(:,1:3)+repmat(sw*[1 0 0],size(mic,1),1);
    sgnUp = (-1).^(mic(:,4)-1);
    insert = repmat(sw*[0.5 sin(pi/3) 0],size(mic,1),1);
    insert(:,2) = insert(:,2).*sgnUp;
    points(:,7:9) = points(:,1:3)+insert;
    points(:,10:12) = points(:,1:3)+repmat([0 0 dZ],size(mic,1),1);
    points(:,13:15) = points(:,4:6)+repmat([0 0 dZ],size(mic,1),1);
    points(:,16:18) = points(:,7:9)+repmat([0 0 dZ],size(mic,1),1); 

    verts = [points(:,1:3);points(:,4:6);points(:,7:9);points(:,10:12);...
            points(:,13:15); points(:,16:18)];
    
    if fAngle ~= 0 
        th = fAngle*pi/180;
        R = [cos(th) -sin(th) 0 ;sin(th) cos(th) 0; 0 0 1];
        verts = (R*verts')';
    end
    
    if sOptions.bApplyPostrotationTranslation
        verts(:,1) = verts(:,1)+sOptions.vPostTranslation(1);
        verts(:,2) = verts(:,2)+sOptions.vPostTranslation(2);
    end

    p1_indx = 1;
    p2_indx = 1+size(points(:,1:3),1);
    p3_indx = 1+2*size(points(:,1:3),1);
    p4_indx = 1+3*size(points(:,1:3),1);
    p5_indx = 1+4*size(points(:,1:3),1);
    p6_indx = 1+5*size(points(:,1:3),1);
    
    faces = zeros(size(mic,1)*3,4);
  
    for i = 1:3:size( mic,1 )*3
        faces(i,:) = [p1_indx p4_indx p5_indx p2_indx];
        faces(i+1,:) = [p2_indx p5_indx p6_indx p3_indx];
        faces(i+2,:) = [p1_indx p4_indx p6_indx p3_indx];
        p1_indx = p1_indx+1;
        p2_indx = p2_indx+1;
        p3_indx = p3_indx+1;
        p4_indx = p4_indx+1;
        p5_indx = p5_indx+1;
        p6_indx = p6_indx+1;
    end
        
    faceColors = zeros(size(ColorVec,1)*3,3);
    indx=1;
    for i = 1:3:size( faceColors , 1) 
        faceColors(i,:) = ColorVec(indx,:);
        faceColors(i+1,:) = ColorVec(indx,:);
        faceColors(i+2,:) = ColorVec(indx,:);
        indx=indx+1;
    end
    

    if bPlotEdges        
        s = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','FaceAlpha',opacity,'EdgeAlpha',sOptions.fEdgeAlpha);
    else                
        s = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceColors(:,:),'DiffuseStrength',0.75,'FaceColor','flat','FaceLighting','gouraud','EdgeColor','none','FaceAlpha',opacity);
    end

end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Legacy plot_mic, mostly untouched.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_mic( varargin )
    
    
    sRGB = varargin{5};     
    snp       = varargin{1};  
    sidewidth = varargin{2};
    plotType  = varargin{3};
    ConfidenceRange = varargin{4};    
    
    if (plotType ~= 9)
        [tri_x, tri_y, snp_Euler, snp_Conf ] = GetTriangleVerticesLegacy( snp, sidewidth );
    else
        [tri_x, tri_y, snp_Euler, snp_Conf ] = GetTriangleVerticesLegacy( snp, sidewidth ,1);
    end
    
    if(plotType == 2) % plot confidence map
        tmp = [ snp_Conf ];
        tmp = [ snp_Conf ] - min( snp_Conf );   % original auto rescale
        tmp = tmp /( max(snp_Conf) - min(snp_Conf) );  % original auto rescale
        tmp = [tmp, zeros(length(tmp), 1), 1-tmp];
        xConf = [0:0.01:1]';
        colormap( [ xConf, xConf * 0, 1- xConf]  );
        caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
        colorbar( 'location', 'eastoutside');

    elseif(plotType == 3) % plot rodrigues vectors            
        tmp = RodOfQuat( QuatOfRMat( RMatOfBunge(snp_Euler', 'degrees') ) )';
        if strcmp(sRGB,'periodic') == 1
            x = 0.8284/2;   
            RGB = [abs(sin(pi*tmp(:,1)/x)) abs(cos(pi*tmp(:,2)/x)) abs(sin(pi*tmp(:,3)/x))];
            tmp = RGB;
        elseif strcmp(sRGB,'unscaledCFZ') == 1
            x = 0.8284/2;   
            RGB = [tmp(:,1)/(2*x)+0.5 tmp(:,2)/(2*x)+0.5 tmp(:,3)/(2*x)+0.5];
            RGB( RGB(:,1) < 0 ,1) = 0;
            RGB( RGB(:,2) < 0 ,2) = 0;
            RGB( RGB(:,3) < 0 ,3) = 0;
            tmp = RGB;
        end  
        
    elseif(plotType == 9)
        tmp = [ snp_Conf ];
        tmp = [ snp_Conf ] - min( snp_Conf );   % original auto rescale
        tmp = tmp /( max(snp_Conf) - min(snp_Conf) );  % original auto rescale
        tmp = [tmp, zeros(length(tmp), 1), 1-tmp];
        xConf = [0:0.01:1]';
        colormap( [ xConf, xConf * 0, 1- xConf]  );
        caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
        caxis( [min(snp_Conf), max(snp_Conf) ] );
        colorbar( 'location', 'eastoutside');      

    elseif( plotType == 11 ) 
        tmp = [ snp_Conf ];
        tmp = tmp - min( snp_Conf );
        tmp = tmp /( max(snp_Conf) - min(snp_Conf) );  % original auto rescale
        j = hot;
        j = DivideColormapNTimes( j , 4) ;
        m = size(j,1);
        tmp = j(floor((tmp*(m-1)))+1,:);        
        colormap( j );
        caxis( [min(snp_Conf), max(snp_Conf) ] );   % original, auto rescale
        colorbar
    end
  
    size_vec = size(tmp);
    tri_color = zeros(1, size_vec(1), size_vec(2));
    tri_color(1, :, :) = tmp;
    
    
    h = patch(tri_x, tri_y, tri_color, 'EdgeColor', 'none');

    if plotType == 9
        nMinSize = 1200;
        gids = unique(snp(:,11));
        for i=1:length(gids)
            ind = find( snp(:,11) == gids(i) );
            if size(ind,1) > nMinSize
                x_avg = mean( snp(ind,1) );
                y_avg = mean( snp(ind,2) );
                text(x_avg,y_avg, num2str(gids(i)), 'Color', 'w','FontSize',12);
            end
        end
    end
    box on;
    axis square;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Get Verticies for the patches
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tri_x, tri_y, snp_Euler, snp_Conf , snp_Euler_Unordered ] = GetTriangleVertices( varargin )

    sOptions = varargin{1};    
    snp = sOptions.mic;
    if sOptions.plotType == 5
        bKAM = 1;        
        ind = find( snp(:,10) > sOptions.fConfThresh );
        snp = snp(ind,:);
    else
        bKAM = 0;
    end
    snp_Euler_Unordered = snp( :, 7:9 );
    sidewidth = sOptions.sideWidth;
        
    if sOptions.plotType == 9
        bGID = 1;
    else 
        bGID = 0;
    end
    
    if sOptions.plotType == 5
        bKAM = 1;           
    else
        bKAM = 0;
    end
    
    if sOptions.plotType == 6
        bScalars = 1;
        snp(:,10) = sOptions.colorScalars;
    else
        bScalars = 0;
    end
    
    if bGID
        snp_Euler_Unordered = repmat(snp( :, 11 ),1,3);                
        if sOptions.fRotation ~= 0
            th = sOptions.fRotation*pi/180;
            R = [cos(th) -sin(th);sin(th) cos(th)];
            xys = (R*snp(:,1:2)')';
            snp_Euler_Unordered(:,2:3) = xys;
        else
            snp_Euler_Unordered(:,2:3) = snp(:,1:2);
        end
        if sOptions.bApplyPostrotationTranslation
            snp_Euler_Unordered(:,2) = snp_Euler_Unordered(:,2) + ...
                sOptions.vPostTranslation(1);
            snp_Euler_Unordered(:,3) = snp_Euler_Unordered(:,3) + ...
                sOptions.vPostTranslation(2);
        end        
    end
    %{
    if bScalars
        snp_Euler_Unordered = repmat( sOptions.colorScalars , 1 , 3) ;
        if sOptions.fRotation ~= 0
            th = sOptions.fRotation*pi/180;
            R = [cos(th) -sin(th);sin(th) cos(th)];
            xys = (R*snp(:,1:2)')';
            snp_Euler_Unordered(:,2:3) = xys;
        else
            snp_Euler_Unordered(:,2:3) = snp(:,1:2);
        end
        if sOptions.bApplyPostrotationTranslation
            snp_Euler_Unordered(:,2) = snp_Euler_Unordered(:,2) + ...
                sOptions.vPostTranslation(1);
            snp_Euler_Unordered(:,3) = snp_Euler_Unordered(:,3) + ...
                sOptions.vPostTranslation(2);
        end        
    end
    %}
    snp(:,5) = 2.^ snp(:, 5);
    
    % find range
    max_x = 1*sidewidth;
    max_y = 1*sidewidth;
    min_x = -1*sidewidth;
    min_y = -1*sidewidth;

    % find all fitted triangles
    snp = sortrows(snp, 6);
    findvec = find(snp(:, 6) > 0);
    snp = snp(findvec, :);

    % sort by triangle
    snp = sortrows(snp, 4);

    % find triangles pointing up
    downsIndex = find(snp(:, 4) > 1);
    upsIndex = find(snp(:, 4) <= 1);

    ups = snp(upsIndex, :);
    downs = snp(downsIndex, :);

    % overly general:

    ups_sides = sidewidth ./ ups(:, 5);
    downs_sides = sidewidth ./ downs(:, 5);

    % calculate ups v1, ups v2, and ups v3
    ups_v1 = ups(:, 1:2);      % (x, y)
    ups_v2 = ups(:, 1:2);
    ups_v2(:, 1) = ups_v2(:, 1) + ups_sides;  % (x+s, y) direction
    ups_v3 = ups(:, 1:2);
    ups_v3(:, 1) = ups_v3(:, 1) + ups_sides/2; % (x+s/2, y)
    ups_v3(:, 2) = ups_v3(:, 2) + ups_sides/2 * sqrt(3); % (x+s/2, y+s/2 *sqrt(3));


    % calculate downs v1, downs v2, and downs v3
    downs_v1 = downs(:, 1:2);      % (x, y)
    downs_v2 = downs(:, 1:2);
    downs_v2(:, 1) = downs_v2(:, 1) + downs_sides;  % (x+s, y) direction
    downs_v3 = downs(:, 1:2);
    downs_v3(:, 1) = downs_v3(:, 1) + downs_sides/2; % (x+s/2, y)
    downs_v3(:, 2) = downs_v3(:, 2) - downs_sides/2 * sqrt(3); % (x+s/2, y - s/2 *sqrt(3));
    
    if sOptions.fRotation ~= 0
        th = sOptions.fRotation*pi/180;
        R = [cos(th) -sin(th);sin(th) cos(th)];
        c = R*ups_v1';
        ups_v1(:,1) = c(1,:)';
        ups_v1(:,2) = c(2,:)';
        c = R*ups_v2';
        ups_v2(:,1) = c(1,:)';
        ups_v2(:,2) = c(2,:)';

        c = R*ups_v3';
        ups_v3(:,1) = c(1,:)';
        ups_v3(:,2) = c(2,:)';

        c = R*downs_v1';
        downs_v1(:,1) = c(1,:)';
        downs_v1(:,2) = c(2,:)';

        c = R*downs_v2';
        downs_v2(:,1) = c(1,:)';
        downs_v2(:,2) = c(2,:)';

        c = R*downs_v3';
        downs_v3(:,1) = c(1,:)';
        downs_v3(:,2) = c(2,:)';
    end
    
    % format is in [v1; v2; v3], where v1, v2, and v3 are rol vectors
    tri_x = [ [ups_v1(:, 1); downs_v1(:, 1)]'; [ups_v2(:, 1); downs_v2(:, 1)]'; [ups_v3(:, 1); downs_v3(:, 1)]'];
    tri_y = [ [ups_v1(:, 2); downs_v1(:, 2)]'; [ups_v2(:, 2); downs_v2(:, 2)]'; [ups_v3(:, 2); downs_v3(:, 2)]'];

    if sOptions.bApplyPostrotationTranslation
        tri_x = tri_x+sOptions.vPostTranslation(1);
        tri_y = tri_y+sOptions.vPostTranslation(2);
    end
    
    snp_Reordered = [ ups; downs];
    snp_Euler = snp_Reordered( :, 7:9 );
    
    if bGID 
        snp_Conf  = snp_Reordered( :, 11  );                  
    elseif bKAM
        snp_Conf = snp_Reordered( : , 18 );
    elseif bScalars 
        snp_Conf = snp_Reordered( :, 10 );
    else
        snp_Conf  = snp_Reordered( :, 10  );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Get Verticies (Legacy) for the patches
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tri_x, tri_y, snp_Euler, snp_Conf ] = GetTriangleVerticesLegacy( varargin )

    snp = varargin{1};
    sidewidth = varargin{2};

    if nargin == 3
        bGID = varargin{3};
    else 
        bGID = 0;
    end

    snp(:,5) = 2.^ snp(:, 5);
    
    % find range
    max_x = 1*sidewidth;
    max_y = 1*sidewidth;
    min_x = -1*sidewidth;
    min_y = -1*sidewidth;

    % find all fitted triangles
    snp = sortrows(snp, 6);
    findvec = find(snp(:, 6) > 0);
    snp = snp(findvec, :);

    % sort by triangle
    snp = sortrows(snp, 4);

    % find triangles pointing up
    downsIndex = find(snp(:, 4) > 1);
    upsIndex = find(snp(:, 4) <= 1);

    ups = snp(upsIndex, :);
    downs = snp(downsIndex, :);

    % overly general:

    ups_sides = sidewidth ./ ups(:, 5);
    downs_sides = sidewidth ./ downs(:, 5);

    % calculate ups v1, ups v2, and ups v3
    ups_v1 = ups(:, 1:2);      % (x, y)
    ups_v2 = ups(:, 1:2);
    ups_v2(:, 1) = ups_v2(:, 1) + ups_sides;  % (x+s, y) direction
    ups_v3 = ups(:, 1:2);
    ups_v3(:, 1) = ups_v3(:, 1) + ups_sides/2; % (x+s/2, y)
    ups_v3(:, 2) = ups_v3(:, 2) + ups_sides/2 * sqrt(3); % (x+s/2, y+s/2 *sqrt(3));


    % calculate downs v1, downs v2, and downs v3
    downs_v1 = downs(:, 1:2);      % (x, y)
    downs_v2 = downs(:, 1:2);
    downs_v2(:, 1) = downs_v2(:, 1) + downs_sides;  % (x+s, y) direction
    downs_v3 = downs(:, 1:2);
    downs_v3(:, 1) = downs_v3(:, 1) + downs_sides/2; % (x+s/2, y)
    downs_v3(:, 2) = downs_v3(:, 2) - downs_sides/2 * sqrt(3); % (x+s/2, y - s/2 *sqrt(3));

    % format is in [v1; v2; v3], where v1, v2, and v3 are rol vectors
    tri_x = [ [ups_v1(:, 1); downs_v1(:, 1)]'; [ups_v2(:, 1); downs_v2(:, 1)]'; [ups_v3(:, 1); downs_v3(:, 1)]'];
    tri_y = [ [ups_v1(:, 2); downs_v1(:, 2)]'; [ups_v2(:, 2); downs_v2(:, 2)]'; [ups_v3(:, 2); downs_v3(:, 2)]'];

    snp_Reordered = [ ups; downs];
    snp_Euler = snp_Reordered( :, 7:9 );
    if bGID 
        snp_Conf  = snp_Reordered( :, 11  );
    else
        snp_Conf  = snp_Reordered( :, 10  );
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Helpers
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function j2 = DivideColormapNTimes(j,n) 
    i = 0;
    while i~=n
        j = divide_colormap(j);
        i = i+1;
    end
    j2=j;
end

function j2 = divide_colormap(j)
   j2 = zeros( 2*size(j,1) - 1, 3 );
    for i=1:(size(j,1)-1)
            j2(2*i-1,:) = j(i,:);
            j2(2*i,:) = 0.5*(j(i,:)+j(i+1,:));
    end
end


function ind = GetMicIndicesInRadius( mic, fRadius , x0 , y0 )
    ind = find( sqrt( (mic(:,1)-x0).^2 + (mic(:,2)-y0).^2 ) < fRadius );
end
  