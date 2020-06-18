% CODEGEN_SCRIPT   Generate MEX-function dz_free_space_mex from dz_free_space.
% 
% Script generated from project 'dz_free_space.prj' on 18-Jun-2020.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;
cfg.IntegrityChecks = false;

%% Define argument types for entry-point 'dz_free_space'.
ARGS = cell(1,1);
ARGS{1} = cell(6,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[3 3]);

%% Invoke MATLAB Coder.
codegen -config cfg dz_free_space -args ARGS{1}

