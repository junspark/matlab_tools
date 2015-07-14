% stock_feb14_saxs
% use to normalize data, cartesian-polar transform and output data for
% reconstructions
clear all;

%% user defined parameters
fileroot = 'pixiradtest';  file_numbers=248926; %full set: 4409 to 9465 files to be evaluated

use_metadata=0; %import metadata associated with files, saved in seperate file
metadata_file='fastpar_stock_feb14_SX.par';
samp_thick_cm=0.2; %sample thickness in cm (1 is default)
correct_atten=1;
atten.x = 1037-3; atten.y = 678; % in pix,pix - for atten
atten.height_mm=2; %attenuator apex, in mm
atten.radius_mm=2.1; %attenuator radius, in mm
atten.material='Cu';
subtract_back=1;
backfile='por_sp_ref2_04410.tif';
darkfile='por_sp_ref2_04406.tif';
scale_abs=0;
trans=0.75; %sample transmission, for background subtraction
abs_int_factor=0.31; %data multiplier to get data on same scale as glassyC - don't change this!
plot_raw=1; %set to 1 to plot raw 2d data
size_plot_raw=300; %size of centered image to plot
plot_corr=0; %set to 1 to plot raw 2d data
size_plot_corr=300; %size of centered image to plot
plot_transform=0 ; %set to 1 to plot 2d transformed image
plot_lineouts=0;
q_min=0.007; %min value in q to plot/output, in 1/A, 0.005 for normal segments
calc_lineouts=1;
azi_mins=[0:20:340]; azi_maxs=[20:20:360];%first/last azimuth values, in degrees, to calculate output for later reconstructions

overlay_lineouts=0;
iquad_overlay=0; %not used in present form!
output_lineouts=0; %set to 1 to output data from a single azimuth
fit_saxs_data=0;
output_fits=0;% set to 1 to output fitted data, added output if metadata used
save_workspace=0; %set to 1 to output resulting workspace variables in fileroot- & date-specific mat file, for retreival in subsequent matlab sessions
symbols={'bo-' 'ro-' 'go-' 'ko-' 'b--' 'r--' 'g--' 'k--'};
% call program for transformation lookup table or load this information if it already exists
if ~exist('stock_feb14_saxs_trans_info.mat'); stock_feb14_set_saxs; else load('stock_feb14_saxs_trans_info.mat'); end;
%define lineout angles to plot; needs to be done after loading mat file b/c needs eta.N variable
%lineout_azi_bins=eta.N/8+(-1:1); %45 deg segments
lineout_azi_bins=eta.N/4+(-2:2); %normal segments
%% main program

%compute different x-axis values
radius_mm=pix.size*(rho.min+0.5:rho.N/(rho.max-rho.min):rho.max-0.5);
two_theta=atan(radius_mm/center.z); d_A=lambda./2./sin(two_theta./2); q_1_A=2*pi./d_A;
xcen=center.x; ycen=pix.y-center.y; %center regions for plotting
lineout_azi_deg=lineout_azi_bins.*360/eta.N;
azi_index_mins=max(1,round(azi_mins/eta.int)); %starting azi bin to plot/output; if 0 or negative vals input default to 1
azi_index_maxs=max(azi_index_mins,round(azi_maxs/eta.int)); %ending azi bin to plot/output; if less than azi_index_min set to that val (1bin avg)

if correct_atten;
    atten_trans=calc_conical_atten_rect(pix.x,pix.y,pix.size,energy,atten.x,atten.y,atten.height_mm,atten.radius_mm,atten.material);
else
    atten_trans=ones(pix.x,pix.y);
end

%% set up metadata information, if selected
if use_metadata
    % read parameter data which can be used to correlate x-ray data with sample transmission, keyence, loading, position, etc
   fid=fopen(metadata_file);  
   met=textscan(fid,'%s %s %s %s %s %s %f %f %f %f %f %s %f %f %s %f %f %f %f %f %f %f %f %f %s %f\n');
   meta0.samX=met{(7)}; meta0.samZ=met{(11)}; meta0.samY=met{(9)}; meta0.phi=met{(18)}; 
   meta0.i0=met{(20)}; meta0.i1=met{(21)}; meta0.ifile=met{(16)};meta0.GEfile=met{(26)};
end;


if overlay_lineouts; figure(30);clf; end;

if output_fits; fout_fits=[fileroot '_' num2str(min(file_numbers)) '_' num2str(max(file_numbers)) '_' date]; fid_out2=fopen(fout_fits,'w+'); end;

% read in, transform ('cake') and fit selected data files
iifile=0;
for ifile=file_numbers; iifile=iifile+1;
    n5 = '00000'; n5(5-length(num2str(ifile))+1:5) = num2str(ifile);
    filename=[fileroot '_'  n5  '.tif'];
    if exist(filename)
        fprintf('------Reading in, transforming and plotting/fitting/outputting data from file %s------\n',filename);
        if output_lineouts;fout_lineouts=[fileroot '_' n5 '.txt']; fopen(fout_lineouts,'w+'); end;
        
        if use_metadata; %assign metadata for a given file to arrays
            irow=find(meta0.ifile==ifile); irow=irow(length(irow));%  find data row corresponding to current file number
            meta.i0(iifile)=meta0.i0(irow); meta.i1(iifile)=meta0.i1(irow); meta.samX(iifile)=meta0.samX(irow);
            meta.samY(iifile)=meta0.samY(irow); meta.samZ(iifile)=meta0.samZ(irow); meta.phi(iifile)=meta0.phi(irow);
        else
            meta.etime(ifile)=1; %needed for abs_int calculation (default - though obviously incorrect in this case)
        end;
        
        im_raw=double(imread(filename)); %read in sample file
        im_raw=im_raw'; %puts image in correct x-y coordinates, relative to Fit2d notation
        im=unzing(im_raw); %filter raw data
        im=im./atten_trans; %apply attenuator correction, if user-specified
        
        if subtract_back;
            im_back=unzing(double(imread(backfile)));
            im_dark=unzing(double(imread(darkfile)));
            im_back_sub=abs((im-im_dark')-(im_back'-im_dark')*trans);
        else
            im_back_sub=im;
        end;
        
        if scale_abs;
            im_back_sub_abs=im_back_sub.*abs_int_factor./meta.etime(ifile)./samp_thick_cm; %put data on absolute scale
        else
            im_back_sub_abs=im_back_sub;
        end
        
        if plot_raw;
            figure(3);clf;imagesc(log(im_raw(round(xcen-size_plot_raw:xcen+size_plot_raw),round(ycen-size_plot_raw:ycen+size_plot_raw))));
            title(sprintf('raw, unzinged image from %s',filename),'Interpreter','None');colorbar;
        end
        if plot_corr;
            figure(4);clf;imagesc(log(im_back_sub_abs(round(xcen-size_plot_corr:xcen+size_plot_corr),round(ycen-size_plot_corr:ycen+size_plot_corr))));
            title(sprintf('attenuator corrected, unzinged image from %s',filename),'Interpreter','None');colorbar;
        end;
        imm = cake_trafo3(im_back_sub_abs(mask), iisparse); % transform data
        imm = reshape(imm, eta.N, rho.N); % transform data
        if plot_transform; figure(6);clf;imagesc(log(imm));title(sprintf('transformed & unzinged image from %s',filename),'Interpreter','None');
            xlabel(sprintf('radial bins, covering %i to %i pixels',rho.min,rho.max));
            ylabel(sprintf('azimuthal bins, covering %i to %i deg',eta.min,eta.max)); colorbar; grid on; pause; end;
        
        if plot_lineouts; figure(78);clf; azi_all=[];
            for i=1:4;   [dum,i_st]=min(abs(q_1_A-q_min));
                x_plot=i_st:length(q_1_A); azi_plot=round((i-1)*round(eta.N/4)+lineout_azi_bins);
                ovflow=find(azi_plot>eta.N); azi_plot(ovflow)=azi_plot(ovflow)-eta.N;
                azi_all=[azi_all, azi_plot];
                plot_handle(i)=semilogy(q_1_A(x_plot),median(imm(azi_plot,x_plot),1),symbols{i}); hold on; end;
            %            plot_handle(i)=semilogy(q_1_A(x_plot),median(imm(:,x_plot),1),symbols{i}); hold on; end;
            xlabel('q, 1/A'); ylabel('I, 1/cm, all azimuths'); grid on;
            %    legend_text={'1', '2', '3', '4'};
            if correct_atten;
                title(sprintf('HE-SAXS lineouts, attenuator corrected, for file %s', filename),'Interpreter','None');
            else
                title(sprintf('HE-SAXS lineouts, no atten correction, for file %s', filename),'Interpreter','None');
            end;
            legend(plot_handle(:),'q1', 'q2','q3','q4');
            text(0.1,0.35,sprintf('HE-SAXS PARAMETERS USED:'),'sc');
            text(0.1,0.3,sprintf('time=%4.1f sec',meta.etime(ifile)),'sc');
            text(0.1,0.25,sprintf('thickness=%4.3f cm',samp_thick_cm),'sc');
            text(0.1,0.2,sprintf('transmission=%4.3f',trans),'sc');
            text(0.1,0.15,sprintf('normalization factor=%3.2f',abs_int_factor),'sc');
            text(0.1,0.1,sprintf('azi range= %3.2f deg',max(lineout_azi_deg)-min(lineout_azi_deg)),'sc');
            text(0.1,0.05,sprintf('background file == %s',backfile),'sc','Interpreter','None');
        end; %if plot_lineouts
        
        if overlay_lineouts; azi_plot2=[];
            [dum,i_st]=min(abs(q_1_A-q_min)); x_plot=i_st:length(q_1_A);
            azi_plot=round((iquad_overlay-1)*round(eta.N/4)+lineout_azi_bins);
            for ii=0:2;
                azi_plot2=[azi_plot2,(ii*round(eta.N/4)+lineout_azi_bins)];
            end
            ovflow=find(azi_plot>eta.N); azi_plot(ovflow)=azi_plot(ovflow)-eta.N;
            lineout(ifile,:)=mean(imm(azi_plot2,x_plot),1);
            %           lineout(ifile,:)=mean(imm(:,x_plot),1);
            figure(30); plot_handle2(ii)=semilogy(q_1_A(x_plot),lineout(ifile,:),'b-');hold on;
            %        text(0.1,0.35,sprintf('HE-SAXS PARAMETERS USED:'),'sc');
            %        text(0.1,0.3,sprintf('time=%4.1f sec',meta.etime(ifile)),'sc');
            %        text(0.1,0.25,sprintf('thickness=%4.3f cm',samp_thick_cm),'sc');
            %        text(0.1,0.2,sprintf('transmission=%4.3f',trans),'sc');
            %        text(0.1,0.15,sprintf('normalization factor=%3.2f',abs_int_factor),'sc');
            %        text(0.1,0.1,sprintf('data from azi bins %3.2f to %3.2f deg',min(lineout_azi_deg),max(lineout_azi_deg)),'sc');
            %        text(0.1,0.05,sprintf('background file == %s',backfile),'sc','Interpreter','None');
        end; %if overlay_lineouts
        
        if fit_saxs_data; stock_feb14_fit_saxs; end;
        
        if calc_lineouts; %these can be used for later reconstructions
            for iazi=1:length(azi_mins);
                azi_indeces=azi_index_mins(iazi):azi_index_maxs(iazi);
                azi_indeces(find(azi_indeces>round(360/eta.int)))=azi_indeces(find(azi_indeces>round(360/eta.int)))-round(360/eta.int); %account for values >360 degrees
                meanint(iifile,:,iazi)=mean(imm(azi_indeces,:),1); %average over user-defined azimuthal bins
                normint(iifile,:,iazi)=meanint(iifile,:,iazi)./max(max(meanint(iifile,:,iazi)));
            end;
        end; %if calc_lineouts
        
        if output_lineouts;
            for ii=1:length(lineout);
                fprintf(fid_out,'%f %f \n',q_1_A(ii+i_st-1),lineout(ifile,ii));
            end;
        end; %if output_lineouts==1
        if output_fits;
            if use_metadata;
                fprintf(fid_out2,'%i %f %f %f %f %f %f %f %f %f\n',ifile,meta.samX(ifile),meta.crosshead(ifile),meta.load(ifile),...
                    mean(saxs_fit(ifile,1).cen),mean(saxs_fit(ifile,1).int),mean(saxs_fit(ifile,1).fwhm),...
                    mean(saxs_fit(ifile,2).cen),mean(saxs_fit(ifile,2).int),mean(saxs_fit(ifile,2).fwhm));
            else
                fprintf(fid_out2,'%i %f %f %f %f %f %f\n',ifile,...
                    mean(saxs_fit(ifile,1).cen),mean(saxs_fit(ifile,1).int),mean(saxs_fit(ifile,1).fwhm),...
                    mean(saxs_fit(ifile,2).cen),mean(saxs_fit(ifile,2).int),mean(saxs_fit(ifile,2).fwhm));
            end; %if use_metadata
        end; %if output_fits
    else
        fprintf('could not find specified filename %s',filename);
    end % if filename exists
    leg_text2(iifile)={num2str(ifile)};
end % looping over ifiles

if overlay_lineouts; figure(30); xlabel('q,1/A'); ylabel('absolute intensity, 1/cm'); grid on;
    title(sprintf('Attenuator corrected data for quadrants 1-3'),'Interpreter','None');
    %legend(plot_handle2(:),leg_text2(:));
    axis tight; end;

if output_fits;   fprintf('Wrote saxs results  to file %s \n', fout_fits); clear fout_fits; end;

if output_lineouts;
    for ii=1:length(lineout);
        fprintf(fid_out,'%f %f \n',q_1_A(ii),lineout_corr(ifile,ii));
    end;
    fprintf('Wrote saxs lineout to file %s \n', fout_lineouts); clear fout_lineouts;
end; %if output_lineouts==1

if (output_lineouts||output_fits); fclose('all'); end;


%% save the current workspace for future sessions; file and date stamped (clear memory-intense variables first)
if save_workspace;
    clear im imm iisparse mask im_raw atten_trans im_back_sub im_back_sub_abs; %do not save memory-intense variables
    fout=[fileroot num2str(min(file_numbers)) '_' num2str(max(file_numbers)) '_' date];
    save(fout);
    fprintf('Wrote workspace to mat file %s \n', fout); clear fout;
end;