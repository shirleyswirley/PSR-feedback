%-------------------------------
% Set up paths
%-------------------------------
utils_path = '/graid1/shirlleu/PSRfeedback/utils/';
data_path = '/graid1/shirlleu/PSRfeedback/data/';
fig_save_path = '/graid1/shirlleu/PSRfeedback/pdfs_pngs/';

%-------------------------------
% Set up warnings, defaults, utils
%-------------------------------
warning('off','all');
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');
addpath(genpath(utils_path));
