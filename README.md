# 🧬Effective RNA Velocity Analysis: Mathematics, Implementation, and Application
> 💡 **Tip: [Please find **👉 MY BLOG** for an introduction and complete view of the project behind the code in this repository.](https://myhugoblog)**

## Introduction
Main content of this repository is to provide a comprehensive understanding of the mathematics, implementation, and application of RNA velocity analysis. It includes detailed explanations of the mathematical foundations, practical coding examples, and real-world applications in single-cell RNA sequencing data analysis.

This repository is designed to be a valuable resource for researchers and practitioners in the field of computational biology, particularly those interested in RNA velocity analysis. It aims to bridge the gap between theoretical concepts and practical implementation, making it easier for users to apply RNA velocity analysis in their own research.

In addition to the main content listed below, the repository includes supplementary markdown and script files that delve into even deeper, more technical aspects of RNA velocity analysis. They’re a valuable resource for exploring advanced concepts.

## 📚 Table of Contents

> 💡 Click on 🧭Intro for a quick overview, 📝Blog for an intuitive walkthrough, or 🐙Git for code and math details.

| # | Title | 🧭 Intro | 📝 Blog | 🐙 GitHub |
|--:|-------|:--------:|:--------:|:----------:|
| 1 | **Math Derivation of CME-defined Stochastic Model of RNA Velocity** | [🧭](#1-math-derivation-of-cme-defined-stochastic-model-of-rna-velocity) | [📝](https://myhugoblog) | [🐙](./11a.velocity_unraveled/) |
| 2 | **Math Derivation for Steady-State RNA Velocity Model** | [🧭](#2-math-derivation-for-steady-state-rna-velocity-model) | [📝](https://myhugoblog) | [🐙](./11c.velocity_steady_state/) |
| 3 | **Dynamic Model — (1) Math Solutions** | [🧭](#3-dynamic-model--1-math-solutions) | [📝](https://myhugoblog) | [🐙](./11d.velocity_dynamic_model_derivation/) |
| 4 | **Dynamic Model — (2) Parameter Inference (EM Algorithm)** | [🧭](#4-dynamic-model--2-parameter-inference-em-algorithm) | [📝](https://myhugoblog) | [🐙](./11e.velocity_dynamic_model_inference/) |
| 5 | **Dynamic Model — (3) Post Hoc Velocity Graph** | [🧭](#5-dynamic-model--3-post-hoc-velocity-graph) | [📝](https://myhugoblog) | [🐙](./11f.velocity_dynamic_model_posthoc_velocity-graph/) |
| 6 | **Dynamic Model — (4) Latent Time** | [🧭](#6-dynamic-model--4-latent-time) | [📝](https://myhugoblog) | [🐙](./11f1.velocity_dynamic_model_posthoc_latent-time/) |
| 7 | **Dynamic Model — (5) Global Time Normalization** | [🧭](#7-dynamic-model--5-global-time-normalization) | [📝](https://myhugoblog) | [🐙](./11f2.velocity_dynamic_model_posthoc_global-time/) |
| 8 | **Dynamic Model — (6) Computational Implementation & Reproducibility** | [🧭](#8-dynamic-model--6-computational-implementation--reproducibility) | [📝](https://myhugoblog) | [🐙](./11g.velocity_dynamic_model_implement/) |
| 9 | **Dynamic Model — (7) Gillespie Stochastic Simulation Algorithm (SSA)** | [🧭](#9-dynamic-model--7-gillespie-stochastic-simulation-algorithm-ssa) | [📝](https://myhugoblog) | [🐙](./11h.velocity_gillespie_ssa_dyngen/) |
| 10 | **Dynamic Model — (8) Effective scVelo Analysis** | [🧭](#10-dynamic-model--8-effective-scvelo-analysis) | [📝](https://myhugoblog) | [🐙](./11i.Velocity_pipeline/) |


## Contents

💡 **Tip: Please click [📝BLOG](only-for-demo) in the section titles for intuitive walkthrough of the content in each section's [🐙Git](only-for-demo).**

Images below are credited to [**Logan Voss on Unsplash**](https://unsplash.com) unless otherwise noted.

### 1. Math Derivation of CME-defined Stochastic Model of RNA Velocity [📝Blog](https://myhugoblog) | [🐙Git](11a.velocity_unraveled) 

![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11a.velocity_unraveled/featured.jpg)


The stochastic model defined by the Chemical Master Equation (CME) outperforms deterministic ODE models in capturing the inherent stochasticity of single-cell RNA sequencing (scRNA-seq) data. It is actively developed to provide a more accurate representation of feature counts and their underlying biological processes. And it has also enabled the generation of simulated data to evaluate deterministic ODE models and associated data processing methods commonly used in scRNA-seq analysis. Thus I derive the key equations from the paper [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492), a pivotal paper demonstrating the transformative potential of stochastic approaches. 

<br>

### 2. Math derivation for steady-state RNA velocity mode [📝Blog](https://myhugoblog) | [🐙Git](11c.velocity_steady_state)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11c.velocity_steady_state/featured.jpg)

The steady‑state model was the first to enable a mathematical estimation of RNA velocity, and most subsequent methods are modified versions of it or its generalization (the dynamic model in `scVelo`; see our other blogs). It has its limitations and a solid understanding of its underlying mathematics is needed to apply the model effectively. Here, we derive the steady-state model in `scVelo` and `velocyto`.


<br>



### 3. Dynamic RNA velocity model-- (1) math solutions [📝Blog](https://myhugoblog) | [🐙Git](11d.velocity_dynamic_model_derivation)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11d.velocity_dynamic_model_derivation/featured.png)



The steady‑state model’s reliance on true steady states is at odds with [known biophysical behavior](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492). The dynamic model removes this requirement to broaden RNA velocity’s applicability but inevitably introduces new assumptions that may not hold for every dataset. Effective use of the dynamic model therefore demands a clear understanding of its strengths and limitations. Our blog series toward this goal begins by delving into the mathematical foundations of the dynamic model.

<br>

### 4. Dynamic RNA velocity model-- (2) parameter inference [📝Blog](https://myhugoblog) | [🐙Git](11e.velocity_dynamic_model_inference)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11e.velocity_dynamic_model_inference/featured.png)



To effectively apply the dynamic model for revealing RNA velocity in single-cell RNA-seq data, this second installment of our blog series takes a deep dive into its parameter inference using a two-stage EM algorithm. In this approach, latent time is initially assigned using an explicit formula, and then refined through standard optimization during the "Expectation" step of the final EM iteration.

<br>


### 5. Dynamic RNA velocity model-- (3) post hoc velocity graph [📝Blog](https://myhugoblog) | [🐙Git](11f.velocity_dynamic_model_posthoc_velocity-graph)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11f.velocity_dynamic_model_posthoc_velocity-graph/featured.png)



In this third installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq, we start our deep dive into the post hoc computations in scVelo that shape both the visualization and interpretation of RNA velocity.

Here specifically looks into the two key components that are computed post hoc in scVelo to derive the velocity graph:
- cosine similarity in section A)
- exponential kernel transformation in section B)

And how velocity graph are projected onto an embedding, such as UMAP or t-SNE, which is the common way to visualize the inferred RNA velocity in section C)

Finally, the reconstructability score {{< math >}} $r$ {{< /math >}} in section D), which quantifies how well a subset of genes can recapitulate the overall RNA velocity dynamics inferred from the full set of genes.

<br>


### 6. Dynamic RNA velocity model-- (4) latent time [📝Blog](https://myhugoblog) | [🐙Git](11f1.velocity_dynamic_model_posthoc_latent-time)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11f1.velocity_dynamic_model_posthoc_latent-time/featured.png)



In this fourth installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq, we reveal the mathematical foundations of **latent time** in the context of RNA velocity analysis, which enables in-depth grasp and interpretation of the latent time of RNA velocity.

Latent time is a crucial interpretation of RNA velocity independent of pre-defined (never perfect) dimensional reduction and it embedding. And latent time provide expressive visualization of the **temporal dynamics of cell differentiation** according to RNA velocity. 

Here specifically focuses on:
- Identification of root cells in section A)
- Discussion and numeric examples of transport maps in section B)
- Gene and cell dependence of latent time in section C)
- Neighborhood-convolved latent time in section D)

<br>


### 7. Dynamic RNA velocity model-- (5) Global time normalization [📝Blog](https://myhugoblog) | [🐙Git](11f2.velocity_dynamic_model_posthoc_global-time)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11f2.velocity_dynamic_model_posthoc_global-time/featured.png)



In scVelo paper, the first step to calculate the Gene-shared latent time, after inferring parameters of kinetic rates, is that gene-specific time points of well-fitted genes (with a likelihood of at least 0.1) are normalized to a common global (overall) time scale. Global time normalization is to address the challenge that different genes may have different intrinsic timescales for their kinetic processes. Without normalization, genes with faster kinetics would dominate the velocity field, while slower genes would contribute less to trajectory inference.

However, the method of Global Time Normalization is **not described** in the original scVelo paper. Its source codes [`latent_time()`](https://github.com/theislab/scvelo/blob/f89590b5912ff3c47c7a486fd02e45589632e766/scvelo/tools/_em_model_core.py#L779) and [`compute_shared_time()`](https://github.com/theislab/scvelo/blob/f89590b5912ff3c47c7a486fd02e45589632e766/scvelo/tools/_em_model_utils.py#L82) shows a multi-step ad-hoc voting method to calculate a global latent time for the cells by considering the fitted latent times from multiple high-likelihood genes. **[Newer papers](https://pmc.ncbi.nlm.nih.gov/articles/PMC9550177/) identify Global Time Normalization as a key weakness** of the original scVelo implementation. Below are the alternative methods to calculate the global time normalization, which I want to try some time later to potentially address the relative scale of different genes and avoid the assumption of equal full cycle time for all genes.

<br>


### 8. Dynamic RNA velocity model-- (6) Computational handling in implementation [📝Blog](https://myhugoblog) | [🐙Git](11g.velocity_dynamic_model_implement)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11g.velocity_dynamic_model_implement/featured.png)



Although not explicitly stated in the scVelo paper, the dynamic model of RNA velocity relies heavily on the concept of **neighbors**. This is not merely a computational convenience-- it’s a core principle in how scVelo interprets and visualizes RNA velocity data. In this blog, we unveil how scVelo leverages neighborhood information in its computations in section A).

In scVelo, certain functions are stochastic and depend on random seeds, while others are deterministic. Understanding the **seed dependencies** is crucial for reproducibility in RNA velocity analyses.
Section B) provides a comprehensive overview of which scVelo functions rely on random seeds and which do not, essential to navigate the stochastic nature of scVelo functions for better reproducibility.

The **object structure** of scVelo illustrated in section C) focus on how gene-specific parameters and cell-and-gene-specific data are organized, the key to understanding how scVelo models RNA velocity at gene and cell levels.

Together, the insights presented here not only improve the accuracy of RNA velocity interpretation, but also empower users to identify and thoughtfully fine-tune key (hyper)parameters, enabling a more systematic and biologically meaningful analysis of single-cell datasets.


<br>


### 9. Dynamic RNA velocity model-- (7) Gillespie Stochastic Simulation Algorithm [📝Blog](https://myhugoblog) | [🐙Git](11h.velocity_gillespie_ssa_dyngen)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11h.velocity_gillespie_ssa_dyngen/featured.png)



The dynamic RNA velocity model makes assumptions that are not always satisfied in real-world data, which helps explain why default pipelines often fail to capture true RNA velocity. Our previous six blog posts on effectively applying this model have equipped us to identify and thoughtfully fine-tune key (hyper)parameters and processing steps in scVelo, allowing for more accurate and meaningful analyses.

However, to rigorously validate these adjustments, we need high-quality simulated data with known ground truth. That’s why this seventh blog delves into the Gillespie Stochastic Simulation Algorithm (SSA)-- a Monte Carlo method well suitable to simulate biochemical reaction systems where molecule numbers are low and stochastic effects are significant.

SSA is particularly well-suited for single-cell RNA-seq applications because it:
- Models intrinsic noise at low molecule counts
- Reflects realistic cell-to-cell variability
- Handles complex, nonlinear networks
- Captures rare, probabilistic events
- Aligns with the biological reality of gene expression

In the context of RNA velocity, stochastic modeling:
- Better estimate splicing rates and latent time
- Model uncertainty in velocity vectors
- Reflect non-smooth, noisy transitions across cell states

Specifically, this post illustrates Gillespie’s Stochastic Simulation Algorithm (SSA) through mathematical derivation, intuitive explanation, and concrete numerical examples-- providing a solid foundation to grasp the algorithm and use its simulated datasets effectively in RNA velocity analysis.

<br>


### 10. Dynamic RNA velocity model-- (8) Effective scVelo analysis [📝Blog](https://myhugoblog) | [🐙Git](11i.Velocity_pipeline)
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11i.Velocity_pipeline/featured.png)



scVelo, like many other tools, like many other tools, computes mathematical models based on assumptions that are often unmet by real-world single-cell RNA-seq datasets.

Empowered by the in-depth conceptual and implementation foundations unveiled in my earlier blog posts, this post continues to explore effective applications of the dynamic RNA velocity model by:

- Identifying key steps and parameters in scVelo pipeline and math models that are tunable to align with the mathematical foundations of dynamic RNA velocity model and to improve the estimation accuracy
- Benchmarking strategies using simulated datasets with ground-truth velocities generated through state-of-the-art stochastic simulations
- Revealing biologically significant tumor development trajectory by applying the strategies on real-world scRNAseq dataset
- Providing tools to visualize the results and to generate heatmap that labels genes with known biological significance and highly correlated with latent time

The full code is available [here](11i.Velocity_pipeline/scripts).