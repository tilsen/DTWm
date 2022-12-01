
%{
The following is the workflow used to conduct analyses reported in the 
manuscript and to generate all figures.
%}

%The dtw_more toolbox is required (https://github.com/tilsen/DTWm.git):
addpath('/home/tilsen/Projects/toolboxes/st_tools/dtw_more/'); 

%%

%defines paths (must be edited)
tw_helpers;

%generate the datasets:
tw_gen_datasets;

%conduct the LRR analyses:
tw_lrr_batch;

%collect the LRR analyses, exclude outliers, calculate geometric means and
%standard deviations:
tw_collect_lrr;

%%

%makes a dataset of average trajectories and confints for each subject/condition
tw_summarize_traj;

%generates datasets for gam fitting
tw_gen_gamdata;

%{ 
Run gam models in R with: 
    tw_gamfit_batch.R
%}

%extract model predictions (necessary for further plotting)
tw_process_gampred;

%model diagnostics and diagnostic plots:
tw_gam_diagnostics;


