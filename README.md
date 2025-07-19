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

