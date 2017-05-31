% FRAP_analysis_Octave.m v1.0.0
% This script was developed by Agata Witkowska for the study
% Witkowska A. & Jahn R., Biophysical Journal (2017)
% "Rapid SNARE-mediated fusion of liposomes and chromaffin granules with giant unilamellar vesicles"
% http://dx.doi.org/10.1016/j.bpj.2017.03.010


% This Script was developed for the analysis of FRAP data with a circular bleach spot generated with ZEN 2010 software (Zeiss) and is optimized to work with anative file format for this software ".czi" and is intended to use with the accompanying Fiji macro
% It groups multiple FRAP traces into groups with the same bleach parameters and analyses groups of data instead of individual files.
% Required packages: optim, io.

clear all;
close all;

pkg load optim;
pkg load io;

% Import of multiple files version
[DIR, NAME, EXT] = fileparts(mfilename("fullpathext")); % directory, name and extension of this script file
 filelist = textread(fullfile(DIR, 'Log.txt'), '%s', 'Delimiter', ',');
norm2sort = struct;

for f = 1:length(filelist)
  file = fopen(filelist{f}); % open file with results
  [DIRDATA, NAMEDATA, EXTDATA] = fileparts(filelist{f});
  nlines = fskipl(file, Inf); % count number of lines in file
  frewind(file); % set the file position to the beginning
  data = dlmread(file, '\t', [1 2 nlines - 3 3]); % load "Mean(FRAP)" and "Mean(background)" columns to data
  frewind(file); % set the file position to the beginning
  params = dlmread(file, '\t', [nlines - 1 1 nlines 4]); % load "buffer frame", "frame time", "pixel size" and "FRAP radius" to params
  params = num2cell(params); % make separate cells in params
  [pbf ft px fr] = deal(params{:}); % pbf - "prebleach frames"; ft - "frame time"; px - "pixel size"; fr - "FRAP radius"
  clear params; % remove params array
  fclose (file); % close file

  % Two-step normalisation accoriding to Miura, K. (2012) "Analysis of FRAP Curves", European Advanced Light Microscopy Network
  tbleach = pbf * ft;
  Ifrap_pre = 0;
  n = length(data);
  data = [data zeros(n, 2)];
  for t = 1:pbf;
    Ifrap_pre += (data(t, 1) - data(t, 2)) / pbf;
  end
  for t = 1:n;
    data(t, 3) = (data(t, 1) - data(t, 2)) / Ifrap_pre; % first normalisation norm1
  end
  Ifrap_bleach = data(pbf + 1, 3);
  Ifrap_pre_norm = mean(data(1:pbf, 3));
  for t = 1:n;
    data(t, 4) = (data(t, 3) - Ifrap_bleach) / (Ifrap_pre_norm - Ifrap_bleach); % second normalisation norm2
  end

  % Generation of structured array with experiments grouped according to pbf, ft, px, fr and no of timepoints
  grname = strrep(["gr__" num2str(pbf) "__" num2str(ft) "__" num2str(px) "__" num2str(fr) "__" num2str(length(data))], ".", "_"); % group name
  if (!isfield(norm2sort, grname))
    norm2sort.(grname).params = [pbf ft px fr]; % group params pbf ft px fr
    norm2sort.(grname).time = colon(0, n - 1) * ft; % group name time
  endif
  if (isfield(norm2sort.(grname), "all"))
    norm2sort.(grname).all = [norm2sort.(grname).all; transpose(data(:, 4))]; % adds norm2 to the array with all norm2 data from one group
  else
    norm2sort.(grname).all = transpose(data(:, 4)); % creates array with all norm2 data from agroup
  endif
  if (!isfield(norm2sort.(grname), "selected"))
    norm2sort.(grname).selected.data = [];
    norm2sort.(grname).selected.names = {};
  endif
  if (mean(data(pbf + 0.5 * (end - pbf):end, 4)) > 0.15) % select only data with recovery, can be adjusted
    norm2sort.(grname).selected.data = [norm2sort.(grname).selected.data; transpose(data(:, 4))];
    norm2sort.(grname).selected.names{end + 1, 1} = NAMEDATA;
  endif
  norm2sort.(grname).(NAMEDATA) = transpose(data(:, 4));

  data = [transpose(colon(0, n - 1) * ft) data]; % add time to data [time, mean FRAP, mean background, norm 1, norm 2]

  dlmwrite(fullfile(DIRDATA, ["out " NAMEDATA ".txt"]), data, 'delimiter', '\t'); % saves data with normalisation to separate file for each experiment
end

grnames = fieldnames(norm2sort); % structural array with group names

splitpath = strrep(strjoin(strsplit(DIRDATA, filesep)(end - 1:end), "_"), " ", "_");
xlsname = ["Fit_log_" splitpath ".xlsx"];
header = {"Group Name" "A" "err_A" "tau" "err_tau" "D" "chi2" "deg_f" "Q"};
xlsDIRDATA = "C:/Octave"; % workaround the bug in io package
xlswrite(fullfile(xlsDIRDATA, xlsname), header, "Log");

% Analysis on groups
gg = 1;
for g = 1:length(grnames)
  pbf = norm2sort.(grnames{g}).params(1);
  time = norm2sort.(grnames{g}).time;
  norm2 = norm2sort.(grnames{g}).selected.data;
  all = norm2sort.(grnames{g}).all;
  dlmwrite(fullfile(DIRDATA, [grnames{g} "_norms_mean.txt"]), [transpose(time) transpose(mean(all)) transpose(std(all))], 'delimiter', '\t');

  % Plots of groups
  figure('Name', grnames{g})
  subplot (2, 1, 1)
  plot(time, all)
  grimgnames = fieldnames(norm2sort.(grnames{g}))(5:end);
  legend(grimgnames)

  if (length(norm2) ! = 0)
    if (length(norm2(:, 1)) == 1)
      y = norm2;
    else
      y = mean(norm2);
    endif
    sigma = std(y(1:pbf)); % measurement error
 
    % Fitting
    tfit = time(pbf + 1:pbf + 150);
    y = y(pbf + 1:pbf + 150);
    t0 = tfit(1) - 0.0000001;
    % model function:
    f = @ (p, tfit) p(1) * exp(- 2 * p(2) ./ (tfit - t0)) .* (besselj(0, 2 * p(2) ./ (tfit - t0)) + besselj(1, 2 * p(2) ./ (tfit - t0)));
    % initial values:
    init = [0.5; 0.1];
 
    % Linear constraints
    A = [1; - 1]; B = 0;
    settings = optimset ("inequc", {A, B});
 
    % Start optimisation
    [p, model_values, cvg, outp] = nonlin_curvefit (f, init, tfit, y);
 
    % Plot mean data and fit
    subplot (2, 1, 2)
    plot(tfit, y, tfit, model_values)
 
    % Statistics of fit
    chi2 = sum(((y - model_values) / sigma) .^ 2);
    Q = gammainc(0.5 * chi2, 0.5 * (length(y) - length(p)), "upper");
    D = fr ^ 2 / (4 * p(2));
 
    settings = optimset ("ret_covp", true, "objf_type", "wls");
 
    p_err_ar = curvefit_stat (f, p, tfit, y, settings);
 
    p_err = sqrt(diag(p_err_ar.covp));
 
    % Generation of output xls file
    LOG = {grnames{g} p(1) p_err(1) p(2) p_err(2) D chi2 (length(y) - length(p)) Q};
     LOG_cell = sprintf('A%s:end%s', num2str(gg + 1), num2str(gg + 1));
    xlswrite(fullfile(xlsDIRDATA, xlsname), LOG, "Log", LOG_cell);
 
     xlswrite(fullfile(xlsDIRDATA, xlsname), {"Time (s)" "Intensity" "Fit"}, sprintf('group %s', num2str(g)));
     xlswrite(fullfile(xlsDIRDATA, xlsname), transpose([tfit; y; model_values]), sprintf('group %s', num2str(g)), sprintf('A2:end%s', num2str(1 + length(tfit))));
     xlswrite(fullfile(xlsDIRDATA, xlsname), {"Group Name" "All" "Selected"}, sprintf('group %s', num2str(g)), 'E1:end1');
     xlswrite(fullfile(xlsDIRDATA, xlsname), {grnames{g}}, sprintf('group %s', num2str(g)), 'E2');
     xlswrite(fullfile(xlsDIRDATA, xlsname), grimgnames, sprintf('group %s', num2str(g)), sprintf('F2:F%s', num2str(1 + length(grimgnames))));
     xlswrite(fullfile(xlsDIRDATA, xlsname), norm2sort.(grnames{g}).selected.names, sprintf('group %s', num2str(g)), sprintf('G2:G%s', num2str(1 + length(norm2sort.(grnames{g}).selected.names))));
    gg + +;
  else
     xlswrite(fullfile(xlsDIRDATA, xlsname), {"Time (s)" "Intensity"}, sprintf('group %s (-)', num2str(g)));
     xlswrite(fullfile(xlsDIRDATA, xlsname), transpose([tfit; y]), sprintf('group %s (-)', num2str(g)), sprintf('A2:end%s', num2str(1 + length(tfit))));
     xlswrite(fullfile(xlsDIRDATA, xlsname), {"Group Name" "All"}, sprintf('group %s (-)', num2str(g)), 'E1:end1');
     xlswrite(fullfile(xlsDIRDATA, xlsname), {grnames{g}}, sprintf('group %s (-)', num2str(g)), 'E2');
     xlswrite(fullfile(xlsDIRDATA, xlsname), grimgnames, sprintf('group %s (-)', num2str(g)), sprintf('F2:F%s', num2str(1 + length(grimgnames))));
  endif

  % Saving of the output figure in .svg format
  saveas(g, fullfile(xlsDIRDATA, ["Group_no" num2str(g) "_" splitpath ".svg"]));
end