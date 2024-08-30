# COMP 683 Final Project

This is my final project for COMP683, a graduate level CS course at UNC that explores the latest applications of machine learning and data science in biomedical research.

## Background
Trajectory inference is a powerful method used to order cells along developmental trajectories in pseudotime, an axis representing differentiation of the gene profile.
Differential expression analysis applies statistical tests to find genes that are expressed at different levels between two lineages or branches in a developmental trajectory.  
Combined, these methods are powerful for highlighting the role of genes of interest in biological processes and visualizing and understanding the relationships between different cell types

When represented as graphs, biological processes can take the form of several different topologies, such as bifurcating, tree and cyclical.  Creating methods that can accurately represent each of the topologies is an ongoing struggle, and there is a high degree of variation between output trajectories inferred by each model.  Oftentimes, some models infer the wrong topology.  Therefore, it is extremely important to quantify which methods are most useful in which situations and to what extent their predictions can be trusted.

The Dynverse, a family of trajectory inference packages developed by researchers, has made running and evaluating trajectory inference methods painless, but it does not include some of the newest cutting edge methods such as scTEP and Monocle3.  To accurately compare new and future methods, the Dynverse needs to be kept up to date which is one of the challenges I try to tackle here.

TradeSeq is a “modular” trajectory based DE analysis framework that uses an input of cell’s lineages and pseudotimes, but it has not yet been publicly integrated with scTEP or more generally, with Dynverse wrapped trajectories.  This is another challenge I attempt to tackle.

## Contributions

I wrapped 3 trajectory inference methods using Dynwrap: scTEP pathways step included, scTEP without the pathways step, and Monocle3 (I edited an available rough draft implementation of Monocle3 in Dynwrap). Using the Dyno pipeline,  I inferred trajectories using Slingshot, scTEP,  Monocle3, PAGA and PAGA-Tree for 3 datasets: my silver standard, and 2 gold standards (Saelens et. al, described in  more detail later).  I evaluated these trajectories using the Dyneval package and 4 metrics: HIM, F1_branches, F1_milestones and correlation.  I created two methods for using a Dynwrap trajectory upstream of TradeSeq, one which assigned cells to lineages based on the closest vertex, and one which assigned cells to lineages based on the edge the cell was projected onto.  I took implementations of TradeSeq workflows with Slingshot and Monocle3 upstream from the package documentation.  Using TradeSeq’s association test, end test, and pattern test, I calculated differential expressed genes for a synthetic trifurcating dataset produced by Dyngen downstream of Slingshot, Monocle3, MFA and scTEP without pathways.  I then compared these results to the ground truth (provided by Dyngen) to quantify the best upstream trajectory inference method.
