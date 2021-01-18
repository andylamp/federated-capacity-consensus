function [params] = setup_vars(params)
%SETUP_VARS Function responsible for initialising the correct 
% dataset and graph output paths based on OS type while printing
% our current execution configuration.
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 17/11/2020
% 
% License: GPLv3
% 

% check if we have a struct supplied, otherwise create one
if nargin < 1
  params = struct;
end

if ispc
    fprintf("\n !! Detected Windows PC !!\n");
    params.graph_path = ".\graphs\";
else
    fprintf("\n !! Detected Unix-like PC !!\n");
    params.graph_path = "./graphs/";
end

%% Execution parameters

fprintf("\n !! Execution parameters\n");

% check if we print figures
if ~isfield(params, "pflag")
  params.pflag = 1;
end

% check if we print figs
if ~isfield(params, "fig_print")
  params.fig_print = 1;
end

% check if we print pdf's
if ~isfield(params, "pdf_print")
  params.pdf_print = 0;
end

% check if we print png images
if ~isfield(params, "png_print")
  params.png_print = 0;
end

% check if we have an execution type
if ~isfield(params, "type")
  fprintf("\n\t ** WARN: No execution type detected, setting to unknown");
  params.type = "unknown";
end

% spacing
fprintf("\n\n");


%% Figure printing configuration

if params.pflag == 1
  fprintf("\n !! Printing functionality is: ENABLED\n");
  % if we do, make sure the directory is created please note that the name
  % adheres to the ISO-8601 timestamp format tagged with the type of run.
  fp = sprintf("%s%s-%s", params.graph_path, ...
    datestr(now, 'yyyymmddTHHMMSS'), params.type);
  fprintf("\n\t** Trying to create save folder in: %s", fp);
  [s, ~, ~] = mkdir(char(fp));
  if s ~= 1
    fprintf("\n\t!! Error, could not create the folder; disabling figure export\n");
    params.pflag = 0;
  else
    fprintf("\n\t** Folder %s created successfully", fp);
    % update the path with the specific runtime folder
    if ispc
        params.graph_path = strcat(fp, "\");
    else
        params.graph_path = strcat(fp, '/');
    end
    fprintf("\n\t** Target graph dir: %s", params.graph_path);
  end
  
  % check if we are printing a .pdf or .png
  if params.pdf_print == 1
      fprintf("\n\t** Printing output is set to: PDF");
  elseif params.png_print == 1
      fprintf("\n\t** Printing output is set to: PNG");
  end
  
  % check if we are printing .fig
  if params.fig_print == 1
      fprintf("\n\t** Printing to .fig format is: ENABLED\n");
  else
      fprintf("\n\t** Printing to .fig format is: DISABLED\n");
  end
else
  fprintf("\n !! Printing functionality is: DISABLED\n");
end

end

