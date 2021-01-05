% CODEGEN_DZ_FREE_SPACE   Generate MEX-function dz_free_space_mex from
%  dz_free_space.
% 
% Script generated from project 'dz_free_space.prj' on 04-Nov-2020.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;
cfg.SaturateOnIntegerOverflow = false;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;

%% Define argument types for entry-point 'dz_free_space'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[3  3 Inf],[0 0 1]);
ARGS{1}{7} = coder.typeof(false);

%% Invoke MATLAB Coder.
codegen -config cfg dz_free_space -args ARGS{1}

