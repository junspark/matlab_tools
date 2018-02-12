function edd_calibrant(varargin)
% GUI to calibrate detector
%
%   + 1.0  2017/11/03
%          - initial version
%
% Copyright 2017 Andrew Chuang (chuang.cp@gmail.com)
% $Revision: 1.00 $  $Date: 2017/11/03 $

% choice = questdlg('Which detector do you want to calibrate?', ...
%     'Choose Detector', ...
%     'Det-1 (V)','Det-2 (H)','None','Det-1 (V)');
% % Handle response
% switch choice
%     case 'Det-1 (V)'
%         disp([choice ', will be calibrated!!'])
%         detno = 1;
%     case 'Det-2 (H)'
%         disp([choice ', will be calibrated!!'])
%         detno = 2;
%     case 'None'
%         disp('Okay!! Bye....')
%         return;
% end

choice = questdlg('Which type of calibration do you want to do?', ...
    'Choose Calibration', ...
    'Energy Calibration','TOA Calibration','None','Energy Calibration');
% Handle response
switch choice
    case 'Energy Calibration'
        disp(['Do ' choice '!!'])
        h = varargin{1};
        if ~isfield(h,'data')
            warndlg({'Please read the data first!!'});
            return;
        end
        opt.data = h.data.data;
        opt.detno = 1;
        opt.source = 1;
        edd_cali_energy(opt)
    case 'TOA Calibration'
        warndlg({'This is a paid function!!','Please upgrade to the full version.'},'!! Warning !!')
        return;
        disp(['Do ' choice '!!'])
        h = varargin{1};
        if ~isfield(h,'data')
            warndlg({'Please read the data first!!'});
            return;
        end
        opt.data = h.data.data;
        h.config.visopt
        h.config.visopt.Inst(1)
        h.config.visopt.Inst(2)
        opt.detno = 1;
        opt.source = 1;
        edd_cali_TOA(opt)
    case 'None'
        disp('Okay!! Bye....')
        return;
end


% %%%% select which type of calibration
% str = {'Co-57','Cd-109'};
% [selection,~] = listdlg('PromptString','Select type',...
%                 'SelectionMode','single',...
%                 'ListString',str);
% if isempty(selection)            
%     disp("No selection was made\n")
%     return
% end
%         
% switch selection
%     case 1
% %         prompt = {'Enter Channel Number','Enter colormap name:'};
% %         dlg_title = 'Energy = a * Channel + b';
% %         num_lines = 1;
% %         defaultans = {'20','20'};
% %         answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%     case 2
%         disp("Only work for Co-57 as of today!!")
%         return
% end




