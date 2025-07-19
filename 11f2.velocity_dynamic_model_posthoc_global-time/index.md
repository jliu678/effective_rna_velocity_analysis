---
title: 🧬 Dynamic RNA velocity model-- (5) Global time normalization 
summary: scVelo source codes shows a ad-hoc voting method to calculate a global latent time, which was identified as key weaknees. Here present alternative methods to potentially address the relative scale of different genes and avoid the assumption of equal full cycle time for all genes. 
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, Global time normalization
  - Dynamic model
  - scVelo
  - Differential equations
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Introduction

In scVelo paper, the first step to calculate the Gene-shared latent time, after inferring parameters of kinetic rates, is that gene-specific time points of well-fitted genes (with a likelihood of at least 0.1) are normalized to a common global (overall) time scale. Global time normalization is to address the challenge that different genes may have different intrinsic timescales for their kinetic processes. Without normalization, genes with faster kinetics would dominate the velocity field, while slower genes would contribute less to trajectory inference.

However, the method of Global Time Normalization is **not described** in the original scVelo paper. Its source codes [`latent_time()`](https://github.com/theislab/scvelo/blob/f89590b5912ff3c47c7a486fd02e45589632e766/scvelo/tools/_em_model_core.py#L779) and [`compute_shared_time()`](https://github.com/theislab/scvelo/blob/f89590b5912ff3c47c7a486fd02e45589632e766/scvelo/tools/_em_model_utils.py#L82) shows a multi-step ad-hoc voting method to calculate a global latent time for the cells by considering the fitted latent times from multiple high-likelihood genes. **[Newer papers](https://pmc.ncbi.nlm.nih.gov/articles/PMC9550177/) identify Global Time Normalization as a key weakness** of the original scVelo implementation. Below are the alternative methods to calculate the global time normalization, which I want to try some time later to potentially address the relative scale of different genes and avoid the assumption of equal full cycle time for all genes.

## Problem Statement

Each gene {{< math >}} $g$ {{< /math >}} has its own characteristic timescale determined by its kinetic parameters:
- Transcription rate: {{< math >}} $\alpha_g$ {{< /math >}}
- Splicing rate: {{< math >}} $\beta_g$ {{< /math >}}
- Degradation rate: {{< math >}} $\gamma_g$ {{< /math >}}

The raw latent time {{< math >}} $t_{c,g}$ {{< /math >}} for cell {{< math >}} $c$ {{< /math >}}and gene {{< math >}} $g$ {{< /math >}}exists in gene-specific units, making cross-gene comparisons problematic.

## Alternative 1
### Step 1: Gene-Specific Characteristic Timescales

First, we need to identify the characteristic timescale for each gene. This is typically determined by the dominant kinetic process:

**Splicing Timescale**
{{< math >}} $$\tau_{\text{splice},g} = \frac{1}{\beta_g}$$ {{< /math >}}

**Degradation Timescale**  
{{< math >}} $$\tau_{\text{degrade},g} = \frac{1}{\gamma_g}$$ {{< /math >}}

**Combined Characteristic Timescale**

The effective timescale combines both processes:

{{< math >}} $$\tau_g = \frac{1}{\beta_g + \gamma_g}$$ {{< /math >}}

Alternatively, some implementations use:

{{< math >}} $$\tau_g = \frac{1}{\max(\beta_g, \gamma_g)}$$ {{< /math >}}

### Step 2: Raw Latent Time Normalization

The raw latent time {{< math >}} $t_{c,g}$ {{< /math >}} is normalized by the gene's characteristic timescale:

{{< math >}} $$\tilde{t}_{c,g} = \frac{t_{c,g}}{\tau_g}$$ {{< /math >}}

This gives dimensionless, normalized latent times where:
- {{< math >}} $\tilde{t} = 0$ {{< /math >}}: Beginning of kinetic process
- {{< math >}} $\tilde{t} = 1$ {{< /math >}}: One characteristic time unit has passed

### Step 3: Global Time Coordinate

To obtain a single global time coordinate for each cell, we aggregate across genes. Several approaches are used:

#### Method 1: Weighted Average
{{< math >}} $$T_c = \frac{\sum_{g} w_g \cdot \tilde{t}_{c,g}}{\sum_{g} w_g}$$ {{< /math >}}

where weights {{< math >}} $w_g$ {{< /math >}} are typically based on gene reliability:

{{< math >}} $$w_g = \text{fit\_likelihood}_g \cdot \mathbb{I}(\text{fit\_likelihood}_g > \theta)$$ {{< /math >}}

with {{< math >}} $\mathbb{I}$ {{< /math >}} being the indicator function and {{< math >}} $\theta$ {{< /math >}} a threshold (e.g., 0.1).

#### Method 2: Median-Based Approach
{{< math >}} $$T_c = \text{median}_g(\tilde{t}_{c,g})$$ {{< /math >}}

This is more robust to outlier genes with extreme kinetic parameters.

#### Method 3: Principal Component Approach
For genes with high likelihood scores, compute the first principal component:

{{< math >}} $$T_c = \text{PC1}(\{\tilde{t}_{c,g}\}_{g \in G_{\text{high}}})$$ {{< /math >}}

where {{< math >}} $G_{\text{high}} = \{g : \text{fit\_likelihood}_g > \theta\}$ {{< /math >}}

### Step 4: Final Normalization

The global time is often rescaled to {{< math >}} $[0,1]$ {{< /math >}} interval:

{{< math >}} $$T_c^{\text{norm}} = \frac{T_c - \min_c(T_c)}{\max_c(T_c) - \min_c(T_c)}$$ {{< /math >}}

## Alternative 2: Velocity-Based Normalization

We may normalize based on velocity magnitudes:

{{< math >}} $$\sigma_g = \sqrt{\text{Var}(v_{c,g})}$$ {{< /math >}}

{{< math >}} $$\tilde{v}_{c,g} = \frac{v_{c,g}}{\sigma_g}$$ {{< /math >}}

Then use normalized velocities to infer global time through:

{{< math >}} $$T_c = \arg\max_{t} \sum_g \mathcal{L}(v_{c,g}(t) | \tilde{v}_{c,g})$$ {{< /math >}}

## Practical Considerations

### Robust Estimation
To handle outliers in kinetic parameters:

{{< math >}} $$\tau_g^{\text{robust}} = \text{median}_{g'} \left(\frac{1}{\beta_{g'} + \gamma_{g'}}\right) \cdot \frac{\beta_g + \gamma_g}{\text{median}_{g'}(\beta_{g'} + \gamma_{g'})}$$ {{< /math >}}

### Confidence Weighting
Incorporate parameter uncertainty:

{{< math >}} $$w_g = \frac{\text{fit\_likelihood}_g}{\text{std}(\alpha_g, \beta_g, \gamma_g)} \cdot \mathbb{I}(\text{fit\_likelihood}_g > \theta)$$ {{< /math >}}




