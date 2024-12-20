function microstructure = parseGrainData_MultiLayer_ff(pname_all, qsym, varargin)
% parseGrainData_MultiLayer_ff
%   Parse the hedm grain file for multiple layers from ge-linescan.
%   Currently, works for FF-HEDM data reduced using MIDAS only.
%
%   microstructure = parseGrainData_MultiLayer(pname_all, qsym) reads the grain
%   log file and returns the microstructure information
%
%   INPUT:
%
%   pname_all
%       full file path of the grain log file generated by FF-HEDM MIDAS
%       code
%
%   qsym
%       Symmetry operators in quaternions
%
%   These arguments can be followed by a list of
%   parameter/value pairs. Options are:
%
%   'Technique'     far-field (FF) or near-field (NF). default is FF.
%   
%   IN FF-HEDM mode:
%   'CrdSystem'     coordinate system in the log file (default is APS).
%	'LabToSample'	rigid body rotation to bring the sample frame coincide
%					with the laboratory frame. this option only works with
%					APS CrdSystem for now.
%	'C_xstal' 		single crystal stiffness matrix / tensor to compute the
%					stresses.
%   'OffsetValue'           offset step size for multiple layers. The
%                           offset direction is always in the vertical
%                           direction. Default is zero.
%	'OutputReflectionTable' outputs reflection table if requested. default is false.
%   'ComputeSelfMisoTable'  computes misorientation angle between the
%                           constituent grains. default is false.
%   'ComputeSelfDistTable'  computes COM distance table between the
%                           constituent grains. default is false.
%   'Verbose'               outputs more information on parser progress if
%                           true. default is false.
%   'MIDAS_version'         MIDAS version (default is v7)
%
%   OUTPUT:
%   microstrcture with the following sub-fields
%       grains      = constituent grain information
%       nGrains     = number of grains found
%       miso_table  = self misorientation table (if requested)
%       dist_table  = self COM distance table (if requested)
%
%   nGrains
%       number of grains found in the volume
%
%   miso_table
%       misorientation table in between all the grains found
%
%   dist_table
%       distance table in between all the grains found
%
%   grains field has the following sub-fields
%       GrainID - grain id number designated by MIDAS
%       R - rotation matrix
%       quat - corresponding quaternion
%       rod - corresponding rodrigues vector
%       COM - center of mass position
%       lattprms - lattice parameter
%       DiffPos - positional discrepancy / uncertainty
%       DiffOme - omega discrepancy / uncertainty
%       DiffAngle - angular discrepancy / uncertainty
%       GrainRadius - grain radius
%       Completeness - completeness
%       StrainFab - strain in FABLE formulation
%       Strain - strain in strain gauge formulation
%       StrainFab_vec - vectorized StrainFab
%       Strain_vec - vectorized Strain
%       PhaseNumber - phase number
%       V - []
%       Esam - []
%       Ecry - []
%       F - []
%       ReflectionTable - reflection table with substructure
%       CrdSys - coordinate system
%       StrainFabUnits - StrainFab units *
%       StrainUnits - Strain units *
%       StressFab - stress computed from StrainFab *
%       Stress - stress computed from Strain *
%       StressFab_h - hydrostatic stress from StressFab *
%       StressFab_d - deviatoric stress from StressFab *
%       StressFab_vm - von Mises stress from StressFab *
%       Stress_h - hydrostatic stress from Stress *
%       Stress_d - deviatoric stress from Stress *
%       Stress_vm - von Mises stress from Stress *
%       StrainRMS - strain discrepancy / error from projection
%       C_xstal - single crystal stiffness *
%       * - only computed if requested
%
%   ReflectionTable has the following sub-fields
%       * all fields with Spots_csv prefix is from the SpotsMatrix.csv file
%       Spots_csv_SpotID = spot id designated by MIDAS after spot
%       consolidation
%       Spots_csv_ome = omega
%       Spots_csv_DetHCrd = detector H crd
%       Spots_csv_DetVCrd = detector V crd
%       Spots_csv_ome_raw = omega raw
%       Spots_csv_eta = eta
%       Spots_csv_RingNum = ring unmber
%       Spots_csv_YLab = Y crd of the spot in lab frame
%       Spots_csv_ZLab = Z crd of the spot in lab frame
%       Spots_csv_th = theta
%       Spots_csv_strain_error = strain error per spot
%       Spots_csv_derived_ome_Aero = omega angle in aero convention
%       Spots_csv_derived_tth = tth calculated from th
%       Spots_csv_derived_eta_vff = eta in vff convention
%       Spots_csv_derived_eta_hexrd = eta in HEXRD convention
%       * all fields with RingNr_csv prefix is from the Radius_StartNr_XXX.csv file
%       RingNr_csv_I_integrated = integrated intensity of the spot
%       RingNr_csv_ome = omega
%       RingNr_csv_YCen = Y position of the spot
%       RingNr_csv_ZCen = Z position of the spot
%       RingNr_csv_I_max = max intensity
%       RingNr_csv_ome_min = omega min where the spot was seen
%       RingNr_csv_ome_max = omega max where the spot was seen
%       RingNr_csv_RingRadius = ring radius in um
%       RingNr_csv_th = theta
%       RingNr_csv_eta = eta
%       RingNr_csv_dome = delta omega
%       RingNr_csv_nimg = number of images
%       RingNr_csv_grain_volume = grain volume
%       RingNr_csv_grain_radius = grain radius
%       RingNr_csv_I_pwdr = powder intensity
%       RingNr_csv_sig_r = spot width in radius
%       RingNr_csv_sig_eta = spot width in eta
%       RingNr_csv_derived_ome_Aero = omega in aero convention
%       RingNr_csv_derived_ome_min_Aero = min omega in aero convention
%       RingNr_csv_derived_ome_max_Aero = max omega in aero convention
%       RingNr_csv_derived_tth = 2theta
%       RingNr_csv_derived_eta_vff = vff eta computed from eta
%       RingNr_csv_derived_eta_hexrd = eta in HEXRD convention
%       * all fields with hkls_csv prefix is from the hkls.csv file
%       hkls_csv_hkls = family of crystallographic plane hkl
%       hkls_csv_dspacing = nominal dspacing
%       hkls_csv_qvec = scattering vector for the family of
%       crystallographic plane
%       hkls_csv_th = theta
%       hkls_csv_tth = 2tehta
%       hkls_csv_RingRadius = ring radius
%       num_spots = number of spots in the ReflectionTable
%
%   The columns of the Grains.csv file are:
%       Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2]
%       X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence
%
%       where each row describes a grain
%
%       O[row][col] is the orientation matrix of the grain that takes crystal frame to
%       ESRF lab coordinate system.
%       X,Y,Z define the center of mass coordinate of the grain in ESRF
%       lab coordinate system
%       a, b, c, alpha, beta, gamma are the crystal lattice
%       parameters of the grain (NEED TO DESCRIBE HOW THESE ARE DEFINED)
%       Err1, Err2, Err3
%       MeanRadius is the size of the grain
%       Confidence is the completeness of the grain (number of g-vectors
%       found / number of g-vectors anticipated)
%
%   In the case of nf-HEDM, input file is from Ice9 file postprocessed with
%   segmentation routine. This functionality was originally provided by
%   Dave Menasche at Carnegie Mellon University.
%
%   The columns of the input file are:
%       grain id
%       Center of mass (x = along beam, y = OB, z = up)
%       Average Orientation in Bunge convention that transforms a vector
%       in crystal frame to the laboratory frame
%       Volume
%       AverageConfidence (bug in the segmentation routine)
%       NumberNeighbors
%       IDsofNeighbors
%       MisorsWithNeighbors
%
%   Outstanding issues:
%   a. Rmat - do we need transformation or not?
%
%   b. Additional needs for grain informatoin
%       SHOULD GET THESE AS WELL FROM MIDAS
%       nExpGvec = Number of expected G vectors
%       nMeasGvec = Number of measured G vectors
%       nMeasOnce = Number of G vectors measured once
%       nMeasMore = Number of G vectors measured more than once
%       meanIA = Average internal angle between prediced and measured
%       gvec = G vector table
%       hkl = 3 hkl values

% default options
optcell = {...
    'Technique', 'ff-midas', ...
    'NumFrames', 1440, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'OffsetValue', 0, ...
    'OutputReflectionTable', false, ...
    'ComputeSelfMisoTable', false, ...
    'ComputeSelfDistTable', false, ...
    'Verbose', false, ...
    'MIDAS_version', 'v7', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
switch lower(opts.CrdSystem)
    case 'aps'
        disp('COM / orientations / strains will be in the APS coordinate system.')
    case 'esrf'
        disp('COM / orientations / strains will be in the ESRF coordinate system')
	otherwise
        error('LAB coordinate system unknown')
end

disp(sprintf('The LAB FRAME and SAMPLE FRAME are IDENTICAL WHEN OMEGA = %2.1f deg', opts.LabToSample))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

switch lower(opts.Technique)
    case 'ff-midas'
        num_pname_all   = length(pname_all);
        
        ct  = 1;
        for iii = 1:1:num_pname_all
            microstructure_iii  = parseGrainData_OneLayer_ff(pname_all{iii,1}, qsym, ...
                'NumFrames', opts.NumFrames, ...
                'CrdSystem', opts.CrdSystem, ...
                'LabToSample', opts.LabToSample, ...
                'C_xstal', opts.C_xstal, ...
                'OutputReflectionTable', opts.OutputReflectionTable, ...
                'ComputeSelfMisoTable', opts.ComputeSelfMisoTable, ...
                'ComputeSelfDistTable', opts.ComputeSelfDistTable, ...
                'Verbose', opts.Verbose);
            
            for jjj = 1:1:microstructure_iii.nGrains
                microstructure_iii.grains(jjj).GrainID_AllLayers    = 1e6*iii + microstructure_iii.grains(jjj).GrainID;
                microstructure_iii.grains(jjj).COM_AllLayers        = microstructure_iii.grains(jjj).COM;
                switch lower(opts.CrdSystem)
                    case 'aps'
                        microstructure_iii.grains(jjj).COM_AllLayers(2)      = microstructure_iii.grains(jjj).COM_AllLayers(2) + opts.OffsetValue(iii);
                    case 'esrf'
                        microstructure_iii.grains(jjj).COM_AllLayers(3)      = microstructure_iii.grains(jjj).COM_AllLayers(3) + opts.OffsetValue(iii);
                end
                
                microstructure.grains(ct)   = microstructure_iii.grains(jjj);
                ct  = ct + 1;
            end
            
            microstructure.nGrains(iii) = microstructure_iii.nGrains;
            if opts.ComputeSelfMisoTable
                microstructure.miso_table{iii,1}    = microstructure_iii.miso_table;
            end
            
            if opts.ComputeSelfDistTable
                microstructure.dist_table{iii,1}    = microstructure_iii.dist_table;
            end
        end
        microstructure.nGrains_AllLayers    = sum(microstructure.nGrains);
        microstructure.Technique    = 'ff-midas';
        microstructure.CrdSystem    = opts.CrdSystem;
        microstructure.LabToSample  = opts.LabToSample;
        microstructure.C_xstal      = opts.C_xstal;
    otherwise
        error('Technique unknown.')
end
