# 🧬Effective RNA Velocity Analysis: Mathematics, Implementation, and Application
> 💡 **Tip: [Please find **👉 MY BLOG** for an introduction and complete view of the project behind the code in this repository.](https://myhugoblog)**

## Introduction
Main content of this repository is to provide a comprehensive understanding of the mathematics, implementation, and application of RNA velocity analysis. It includes detailed explanations of the mathematical foundations, practical coding examples, and real-world applications in single-cell RNA sequencing data analysis.

This repository is designed to be a valuable resource for researchers and practitioners in the field of computational biology, particularly those interested in RNA velocity analysis. It aims to bridge the gap between theoretical concepts and practical implementation, making it easier for users to apply RNA velocity analysis in their own research.

In addition to the main content listed below, the repository includes supplementary markdown and script files that delve into even deeper, more technical aspects of RNA velocity analysis. They’re a valuable resource for exploring advanced concepts.

## Contents
### Math Derivation of CME-defined Stochastic Model of RNA Velocity (Blog [Github](11a.velocity_unraveled))

![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11a.velocity_unraveled/featured.jpg)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

The stochastic model defined by the Chemical Master Equation (CME) outperforms deterministic ODE models in capturing the inherent stochasticity of single-cell RNA sequencing (scRNA-seq) data. It is actively developed to provide a more accurate representation of feature counts and their underlying biological processes. And it has also enabled the generation of simulated data to evaluate deterministic ODE models and associated data processing methods commonly used in scRNA-seq analysis. Thus I derive the key equations from the paper [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492), a pivotal paper demonstrating the transformative potential of stochastic approaches. 

### Math derivation for steady-state RNA velocity mode (Blog [Github](11c.velocity_steady_state))
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11c.velocity_steady_state/featured.jpg)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

The steady‑state model was the first to enable a mathematical estimation of RNA velocity, and most subsequent methods are modified versions of it or its generalization (the dynamic model in `scVelo`; see our other blogs). It has its limitations and a solid understanding of its underlying mathematics is needed to apply the model effectively. Here, we derive the steady-state model in `scVelo` and `velocyto`.

### Dynamic RNA velocity model-- (1) math solutions (Blog [Github](11d.velocity_dynamic_model_derivation))
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11d.velocity_dynamic_model_derivation/featured.png)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

The steady‑state model’s reliance on true steady states is at odds with [known biophysical behavior](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492). The dynamic model removes this requirement to broaden RNA velocity’s applicability but inevitably introduces new assumptions that may not hold for every dataset. Effective use of the dynamic model therefore demands a clear understanding of its strengths and limitations. Our blog series toward this goal begins by delving into the mathematical foundations of the dynamic model.

### Dynamic RNA velocity model-- (2) parameter inference (Blog [Github](11e.velocity_dynamic_model_parameter_inference))
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11e.velocity_dynamic_model_inference/featured.png)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

To effectively apply the dynamic model for revealing RNA velocity in single-cell RNA-seq data, this second installment of our blog series takes a deep dive into its parameter inference using a two-stage EM algorithm. In this approach, latent time is initially assigned using an explicit formula, and then refined through standard optimization during the "Expectation" step of the final EM iteration.

### Dynamic RNA velocity model-- (3) post hoc velocity graph (Blog [Github](11f.velocity_dynamic_model_posthoc_velocity-graph))
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11f.velocity_dynamic_model_posthoc_velocity-graph/featured.png)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

In this third blog on effectively applying the dynamic model of RNA velocity, we look into post hoc computed cosine similarity and the exponential kernel that shape the RNA velocity graph and embedding. This begins our deep dive into scVelo’s post hoc computations that determine visualization and interpretation.

### Dynamic RNA velocity model-- (4) latent time (Blog [Github](11f1.velocity_dynamic_model_posthoc_latent-time))
![Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)](11f1.velocity_dynamic_model_posthoc_latent-time/featured.png)

*Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)*

In this fourth installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq, we reveal the mathematical foundations of **latent time** in the context of RNA velocity analysis, which enables in-depth grasp and interpretation of the latent time of RNA velocity.

Latent time is a crucial interpretation of RNA velocity independent of pre-defined (never perfect) dimensional reduction and it embedding. And latent time provide expressive visualization of the **temporal dynamics of cell differentiation** according to RNA velocity. 

Here specifically focuses on:
- Identification of root cells in section A)
- Discussion and numeric examples of transport maps in section B)
- Gene and cell dependence of latent time in section C)
- Neighborhood-convolved latent time in section D)

