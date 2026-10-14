The file `Human-to-mosquito.R` includes the within-human model used to simulate trajectories of parasitemia and male and female gametocytemia, PK–PD models for dihydroartemisinin (DHA), piperaquine (PPQ), and primaquine (PQ), and a transmission model to calculate an individual's infectiousness to mosquitoes.

The file `MDA_ACT.R` contains a PK–PD model that simulates three rounds of mass drug administration (MDA) with DHA and piperaquine. This file must be sourced before running the MDA simulations.

The file `Agent based transmission model.R` contains simulations of 100 MDA trajectories under two treatment regimens (ACT alone or ACT plus 0.25 mg/kg primaquine), as well as scenarios in which the routine treatment strategy switches from ACT alone to ACT plus 0.25 mg/kg primaquine.

The file 'Visualization.R' contains code for visualizing the simulation results.

This code was developed in R (version ≥ 4.2.0). A paired simulation comparing two treatment regimens (ACT alone and ACT plus 0.25 mg/kg primaquine), under either the MDA or routine treatment strategy, typically takes approximately 24 hours to complete on a workstation equipped with an Intel Xeon E5-2673 v4 processor (20 cores, 40 threads) and 64 GB of RAM.

This project is licensed under the BSD 3-Clause License - see the LICENSE file for details.

To reproduce the Figure 2, run the script Human-to-mosquito.R. To reproduce the Figures 4　and　５, first run the script Human-to-mosquito.R and MDA_ACT.R, then run Agent based transmission model.R to generate the simulation results, finally run Visualization.R to produce the figures.
