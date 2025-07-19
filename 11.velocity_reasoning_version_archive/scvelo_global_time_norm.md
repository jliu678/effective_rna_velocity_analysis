---
title: "Scvelo Global Time Norm"
date: 2025-07-10
draft: True
---

# Global Time Normalization in scVelo
This is **not described** in the original scVelo paper, and **newer papers identify this as a key weakness** of the original scVelo implementation and propose improvements.


## Overview

Global time normalization in scVelo addresses the challenge that different genes may have different intrinsic timescales for their kinetic processes. Without normalization, genes with faster kinetics would dominate the velocity field, while slower genes would contribute less to trajectory inference.

## Problem Statement

Each gene {{< math >}} $g$ {{< /math >}} has its own characteristic timescale determined by its kinetic parameters:
- Transcription rate: {{< math >}} $\alpha_g$ {{< /math >}}
- Splicing rate: {{< math >}} $\beta_g$ {{< /math >}}
- Degradation rate: {{< math >}} $\gamma_g$ {{< /math >}}

The raw latent time {{< math >}} $t_{c,g}$ {{< /math >}} for cell {{< math >}} $c$ {{< /math >}}and gene {{< math >}} $g$ {{< /math >}}exists in gene-specific units, making cross-gene comparisons problematic.

## Step 1: Gene-Specific Characteristic Timescales

First, we need to identify the characteristic timescale for each gene. This is typically determined by the dominant kinetic process:

### Splicing Timescale
{{< math >}} $$\tau_{\text{splice},g} = \frac{1}{\beta_g}$$ {{< /math >}}

### Degradation Timescale  
{{< math >}} $$\tau_{\text{degrade},g} = \frac{1}{\gamma_g}$$ {{< /math >}}

### Combined Characteristic Timescale
The effective timescale combines both processes:

{{< math >}} $$\tau_g = \frac{1}{\beta_g + \gamma_g}$$ {{< /math >}}

Alternatively, some implementations use:

{{< math >}} $$\tau_g = \frac{1}{\max(\beta_g, \gamma_g)}$$ {{< /math >}}

## Step 2: Raw Latent Time Normalization

The raw latent time {{< math >}} $t_{c,g}$ {{< /math >}} is normalized by the gene's characteristic timescale:

{{< math >}} $$\tilde{t}_{c,g} = \frac{t_{c,g}}{\tau_g}$$ {{< /math >}}

This gives dimensionless, normalized latent times where:
- {{< math >}} $\tilde{t} = 0$ {{< /math >}}: Beginning of kinetic process
- {{< math >}} $\tilde{t} = 1$ {{< /math >}}: One characteristic time unit has passed

## Step 3: Global Time Coordinate

To obtain a single global time coordinate for each cell, we aggregate across genes. Several approaches are used:

### Method 1: Weighted Average
{{< math >}} $$T_c = \frac{\sum_{g} w_g \cdot \tilde{t}_{c,g}}{\sum_{g} w_g}$$ {{< /math >}}

where weights {{< math >}} $w_g$ {{< /math >}} are typically based on gene reliability:

{{< math >}} $$w_g = \text{fit\_likelihood}_g \cdot \mathbb{I}(\text{fit\_likelihood}_g > \theta)$$ {{< /math >}}

with {{< math >}} $\mathbb{I}$ {{< /math >}} being the indicator function and {{< math >}} $\theta$ {{< /math >}} a threshold (e.g., 0.1).

### Method 2: Median-Based Approach
{{< math >}} $$T_c = \text{median}_g(\tilde{t}_{c,g})$$ {{< /math >}}

This is more robust to outlier genes with extreme kinetic parameters.

### Method 3: Principal Component Approach
For genes with high likelihood scores, compute the first principal component:

{{< math >}} $$T_c = \text{PC1}(\{\tilde{t}_{c,g}\}_{g \in G_{\text{high}}})$$ {{< /math >}}

where {{< math >}} $G_{\text{high}} = \{g : \text{fit\_likelihood}_g > \theta\}$ {{< /math >}}

## Step 4: Final Normalization

The global time is often rescaled to {{< math >}} $[0,1]$ {{< /math >}} interval:

{{< math >}} $$T_c^{\text{norm}} = \frac{T_c - \min_c(T_c)}{\max_c(T_c) - \min_c(T_c)}$$ {{< /math >}}

## Mathematical Implementation in scVelo

The complete normalization procedure:

### Step A: Compute Characteristic Times
{{< math >}} $$\tau_g = \frac{1}{\beta_g + \gamma_g} \quad \forall g$$ {{< /math >}}

### Step B: Normalize Gene-Specific Latent Times
{{< math >}} $$\tilde{t}_{c,g} = \frac{t_{c,g}}{\tau_g} \quad \forall c,g$$ {{< /math >}}

### Step C: Compute Weights
{{< math >}} $$w_g = \begin{cases} 
\text{fit\_likelihood}_g & \text{if } \text{fit\_likelihood}_g > 0.1 \\
0 & \text{otherwise}
\end{cases}$$ {{< /math >}}

### Step D: Global Time Aggregation
{{< math >}} $$T_c = \frac{\sum_{g} w_g \cdot \tilde{t}_{c,g}}{\sum_{g} w_g}$$ {{< /math >}}

### Step E: Final Rescaling
{{< math >}} $$T_c^{\text{global}} = \frac{T_c - \min_c(T_c)}{\max_c(T_c) - \min_c(T_c)}$$ {{< /math >}}

## Alternative Formulation: Velocity-Based Normalization

Some implementations normalize based on velocity magnitudes:

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

## Result

The final global time {{< math >}} $T_c^{\text{global}}$ {{< /math >}} provides:
- Consistent temporal ordering across all cells
- Gene-independent time units
- Robust aggregation of temporal information
- Foundation for accurate velocity field computation

This normalization ensures that the RNA velocity analysis captures the true developmental progression rather than being dominated by genes with particular kinetic timescales.