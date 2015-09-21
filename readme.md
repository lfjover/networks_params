#Scripts to generate the figures in: "Multiple regimes of robust patterns between network structure and biodiversity"

figs_paper.m generates all the figures in the paper using results precomputed and saved in the folder data.

To show how the data was created an example is provided in example_one_run.m. It calculates the dynamics for  100 different matrices in the  matrix ensemble using the same parameter set for all of them. This results in the following plot

![alt tag](https://github.com/lfjover/networks_params/blob/master/bio_one_set.png)

The plots in Figure 3 were created in this way (with a single run for each parameter set). For Figure 2 and 4 each plot is an average of 100 diferent runs where each run uses a different parameters set. To calculate the data in Figures 2 and 4 we used a computer cluster to compute each run in parallel.

Each run is computed using the function pers_one_set.m  which calls pp_integrator_stop_heu.m to integrate  the system based on the stopping time heuristic described in the Methods section.

For Figure 7 pp_integrator_fixed_time.m was used to integrate the system for a fixed time.

The code that generates the paramters used in all the runs is in generate_paramters.m
