# DP.RST

**DP.RST** is an R package implementing **DP-RST (Dirichlet Process Mixture of Random Spanning Trees)** - a Bayesian nonparametric method for clustering **spatial transcriptomics** data when (i) the tissue boundary is **non-convex/irregular** and (ii) the **number of spatial domains is unknown**.  

The method is designed for challenging geometries such as “Swiss-roll” tissues, where Euclidean proximity can incorrectly connect spatially disjoint layers.  

---

## Method overview (what DP-RST does)

DP-RST performs clustering in **two layers**:

1. **Spatially contiguous partitioning (graph + tree cutting).**  
   A spatial graph is constructed to respect the tissue’s intrinsic geometry (via **constrained Delaunay triangulation**), and spatial clusters are formed by **randomly cutting a random spanning tree** induced from that graph. This tree-based construction does **not constrain cluster shapes** and naturally supports non-convex domains.  

2. **Refinement via a Dirichlet process mixture.**  
   The initial spatial regions are **merged probabilistically** based on transcriptomic similarity using a **Dirichlet process (DP)**, enabling **data-driven inference of the number of tissue domains**.  

### Inference
Posterior inference is performed via MCMC; to improve mixing in a multi-modal posterior, DP-RST uses **parallel tempering** and retains samples from the “cold” chain. 
A point estimate of the partition can be obtained using a **least-squares criterion** based on posterior co-clustering probabilities.  

---

## Installation

### From GitHub
```r
# install.packages("remotes")
remotes::install_github("rvaleriya/DP.RST")
```

## Implementation Guide (Vignette)

For a detailed walkthrough of the method and its parameters, please refer to the package vignette. This guide demonstrates how to apply DP-RST to simulated data (such as the “Swiss-roll” example).

To view the vignette after installation, run:

```r
vignette("DP-RST", package = "DP.RST")
```
