
%in order of appearance in manuscript:

%examples of local rate calculation with known rates
fig_process_localrate;

%examples of dissimilar signals
fig_examples_incommensurate;

%comparison of LTW vs. DTW
fig_LTW_DTW_comparison;

%LTW method
fig_LTW;

%examples of distance matrices
fig_distance_matrices;

%examples of warping curves as binary relations
fig_warping_curve;

%DTW alignment methods
fig_alignment_methods;

%local slope constraints
fig_slope_constraints;

%local slope/edge effects
fig_edge_effect;

%global window constraint
fig_global_window;

%just copy to images folder
copyfile(['.' filesep 'figure_code' filesep 'lrr_overview.png'], ...
    ['.' filesep 'images' filesep 'fig_lrr_overview.png']);

%local slope/LRR example
fig_local_slope;

%LRR example 1
fig_temporal_inferences;

%LRR example 2
fig_temporal_structure;

%generates an example of exclusion
fig_anomalous_traj;

%average LRR-gsd
fig_LRR_CVC_both;

%GAM LRR-gsd predictions and differences
fig_LRR_gam_diffs;

%examples of GAM LRR-gsd differences with events
fig_LRR_gam_events;

%alignment of features with gestural/segmental events
fig_LRR_feature_alignment;

%comparison of articulatory vs. acoustic fits
fig_acoustic_artic_comp;

%(appendix) generates a figure with plots of all trajectories for all subjects
fig_subj_traj; 

%(appendix) all alpha/beta features:
fig_features_bysubj;

%renames figures in order of data modified:


