---
title: 🧬 Dynamic RNA velocity model-- (6) Computational handling in implementation  
summary: Here begins the deep dive into scVelo’s computational handling that is only unveiled in its implementation but reshapes visualization and interpretation. We look into neighbor reliance, seed dependency and object structure of scVelo in this sixth blog on effectively applying the dynamic model of RNA velocity. 
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, Neighbor reliance, Seed dependency
  - Dynamic model
  - scVelo
  - Implementation
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Introduction

Although not explicitly stated in the scVelo paper, the dynamic model of RNA velocity relies heavily on the concept of **neighbors**. This is not merely a computational convenience-- it’s a core principle in how scVelo interprets and visualizes RNA velocity data. In this blog, we unveil how scVelo leverages neighborhood information in its computations in section A).

In scVelo, certain functions are stochastic and depend on random seeds, while others are deterministic. Understanding the **seed dependencies** is crucial for reproducibility in RNA velocity analyses.
Section B) provides a comprehensive overview of which scVelo functions rely on random seeds and which do not, essential to navigate the stochastic nature of scVelo functions for better reproducibility.

The **object structure** of scVelo illustrated in section C) focus on how gene-specific parameters and cell-and-gene-specific data are organized, the key to understanding how scVelo models RNA velocity at gene and cell levels.

Together, the insights presented here not only improve the accuracy of RNA velocity interpretation, but also empower users to identify and thoughtfully fine-tune key (hyper)parameters, enabling a more systematic and biologically meaningful analysis of single-cell datasets.

## A) Steps in scVelo that Use Neighbors

### 1. scv.pp.moments(adata, ...)

**Uses:** Nearest neighbors ( `.obsp['distances']`  and  `.obsp['connectivities']` ). It internally calls  `scv.pp.neighbors`  which in turn will call  `scanpy.pp.pca`  to compute PCA when  `adata.obsm['X_pca']`  does not exist.

**Purpose:** To compute first- and second-order moments (mean and uncentered variance) of spliced and unspliced counts for each cell based on its local neighborhood.

The moments computation can be mathematically expressed as:

{{< math >}} 
$$M_s^{(i)} = \sum_{j \in N(i)} w_{ij} \cdot s_j \tag{1}$$
{{< /math >}}

{{< math >}} 
$$M_u^{(i)} = \sum_{j \in N(i)} w_{ij} \cdot u_j \tag{2}$$
{{< /math >}}

where {{< math >}} $M_s^{(i)}$ {{< /math >}} and {{< math >}} $M_u^{(i)}$ {{< /math >}} are the neighborhood-smoothed spliced and unspliced counts for cell {{< math >}} $i$ {{< /math >}}, {{< math >}} $N(i)$ {{< /math >}} represents the neighborhood of cell {{< math >}} $i$ {{< /math >}}, and {{< math >}} $w_{ij}$ {{< /math >}} are the connectivity weights.

These moments are directly used in:
-  `scv.tl.velocity()`  – fits the model per gene
-  `scv.tl.recover_dynamics()`  – if using dynamical model

Thus indirectly in all downstream velocity plotting/inference.

The  `use_raw`  parameter of  `scvelo.tl.recover_dynamics`  controls whether to use neighbors-smoothed moments or raw data as shown in:  `use_raw (bool or None (default: None)) -- if True, use .layers['spliced'], else use moments from .layers['Ms']` 

### 2. scv.tl.velocity_graph(adata)

**Uses:**  `.obsp['connectivities']`  or recomputes internally if not present (deprecated).

**Purpose:** Computes a transition probability matrix (cosine similarity of velocity vectors) between neighboring cells.

The velocity graph computation involves:

{{< math >}} 
$$V_{ij} = \frac{\vec{v}_i \cdot (\vec{x}_j - \vec{x}_i)}{|\vec{v}_i| \cdot |\vec{x}_j - \vec{x}_i|} \tag{3}$$
{{< /math >}}

where {{< math >}} $V_{ij}$ {{< /math >}} represents the transition probability from cell {{< math >}} $i$ {{< /math >}} to cell {{< math >}} $j$ {{< /math >}}, {{< math >}} $\vec{v}_i$ {{< /math >}} is the velocity vector of cell {{< math >}} $i$ {{< /math >}}, and {{< math >}} $\vec{x}_j - \vec{x}_i$ {{< /math >}} is the displacement vector between cells.

**Depends on:** cell-cell proximity and velocity vectors -- neighbor graph ensures computational efficiency and biological relevance.

### 3. scv.tl.latent_time(adata)

**Uses:** Neighbor-smoothed velocity graph to determine pseudotime-like ordering.

**Purpose:** Computes the global latent time based on inferred dynamics and the transition probabilities (again, neighbor-dependent).

The latent time computation can be expressed as:

{{< math >}} 
$$t_{\text{latent}}^{(i)} = \sum_{j} V_{ij} \cdot t_j \tag{4}$$
{{< /math >}}

where {{< math >}} $t_{\text{latent}}^{(i)}$ {{< /math >}} is the latent time for cell {{< math >}} $i$ {{< /math >}}, and the sum is over neighboring cells weighted by their transition probabilities {{< math >}} $V_{ij}$ {{< /math >}}.

### 4. scv.pl.velocity_embedding_stream(...)

**Indirectly uses:**  `velocity_graph` , which depends on the neighbor graph.

**Purpose:** Projects streamlines based on velocity vectors, assuming a local neighborhood flow.

The streamline projection follows the differential equation:

{{< math >}} 
$$\frac{d\vec{x}}{dt} = \vec{v}(\vec{x}) \tag{5}$$
{{< /math >}}

where {{< math >}} $\vec{x}$ {{< /math >}} represents the position in embedding space and {{< math >}} $\vec{v}(\vec{x})$ {{< /math >}} is the interpolated velocity field based on the neighbor-dependent velocity graph.




## B) Random Seed Dependencies of scVelo Functions

### Steps that **depend** on random seed

{{% callout note %}}
These functions contain stochastic elements and will produce different results with different random seeds.
{{% /callout %}}

| Function | Description | Seed Dependency |
|----------|-------------|-----------------|
| `scv.pp.neighbors()` | Builds KNN graph, often used for velocity graph | ✅ **Yes** - Approximate neighbors algorithms (e.g., `method='umap'` or `'hnsw'`) use randomization |
| `scv.tl.velocity_graph()` | Computes transition probabilities based on velocity vectors | ✅ **Yes** - Uses neighbors graph and can have nondeterministic steps |
| `scv.tl.recover_dynamics()` | Fits the dynamical model to splicing kinetics (nonlinear optimization) | ✅ **Yes** - Random initialization in curve fitting procedures |
| `scv.tl.latent_time()` | Estimates latent time based on gene likelihoods | ✅ **Yes** - Depends on recovered dynamics which are stochastic |
| `scv.pl.velocity_embedding_stream()` | Stream plot of velocities on embedding | ✅ **Yes** - Depends on UMAP randomness and interpolation grid seed |
| `scv.tl.velocity_pseudotime()` | Alternative to latent time, based on velocity graph | ✅ **Yes** - Velocity graph contains randomness |
| `sc.tl.umap()` (from Scanpy) | Embedding used for velocity visualization | ✅ **Yes** - UMAP uses random initialization |

{{% callout warning %}}
**Reproducibility Tip**: Set random seeds before running these functions to ensure reproducible results across different runs.
{{% /callout %}}

### Steps that **do not depend** on random seed

{{% callout note %}}
These functions are deterministic and will produce identical results given the same input data.
{{% /callout %}}

| Function | Description | Seed Dependency |
|----------|-------------|-----------------|
| `scv.pp.filter_and_normalize()` | Basic filtering and normalization | ❌ **No** - Purely deterministic operations |
| `scv.pp.moments()` | Computes first and second moments for smoothing | ❌ **No** - Deterministic given neighbors graph |
| `scv.tl.velocity()` | Computes RNA velocity using stochastic/dynamical model | ❌ **Mostly deterministic** - Depends on moments and model parameters |

### Best Practices for Reproducibility

**Setting Random Seeds**

```python
import numpy as np
import random

# Set global random seeds
np.random.seed(42)
random.seed(42)

# For specific functions, you can also set seeds locally
scv.settings.set_figure_params('scvelo', dpi_save=300, dpi=80, transparent=True, fontsize=14, color_map='viridis')
```

### Workflow Recommendations

1. **Set seeds early** in your analysis pipeline
2. **Document seed values** used in your analysis
3. **Test sensitivity** by running with multiple different seeds
4. **Save intermediate results** after computationally expensive stochastic steps

---

## C) scVelo object structure

### Gene-specific parameters (same for all cells)

* `.var['fit_alpha']` - Transcription rate for each gene ({{< math >}} $\alpha_g$ {{< /math >}})
* `.var['fit_beta']` - Splicing rate for each gene ({{< math >}} $\beta_g$ {{< /math >}})  
* `.var['fit_gamma']` - Degradation rate for each gene ({{< math >}} $\gamma_g$ {{< /math >}})
* `.var['fit_t_']` - Switching time for each gene ({{< math >}} $t_{switch,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{genes},)$ {{< /math >}} - one value per gene

### Cell-and-gene-specific data

* `.layers['fit_t']` - Latent time varies by both cell AND gene ({{< math >}} $t_{c,g}$ {{< /math >}})
* `.layers['velocity']` - Velocity varies by both cell AND gene ({{< math >}} $v_{c,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{cells} \times n_{genes})$ {{< /math >}} - different value for each cell-gene combination

### Mathematical representation

For a given gene {{< math >}} $g$ {{< /math >}} and cell {{< math >}} $c$ {{< /math >}}:

{{< math >}}
$$
\begin{align}
\alpha_g &= \text{transcription rate (gene-specific)} \\
\beta_g &= \text{splicing rate (gene-specific)} \\
\gamma_g &= \text{degradation rate (gene-specific)} \\
t_{c,g} &= \text{latent time (cell- and gene-specific)} \\
v_{c,g} &= \text{velocity (cell- and gene-specific)}
\end{align}
$$
{{< /math >}}

### Why this makes sense biologically

* Each gene has **intrinsic kinetic properties** - the rates {{< math >}} $\alpha_g$, $\beta_g$, $\gamma_g$ {{< /math >}} describe how fast the gene transcribes, splices, and degrades
* But different cells are at **different stages** of the process for each gene - the latent time {{< math >}} $t_{c,g}$ {{< /math >}} varies
* The **rates are gene properties**, the **timing/stage varies by cell**

This means:
- If you want to compare kinetic rates between genes → look at `.var`
- If you want to see where different cells are in the process for a specific gene → look at `.layers['fit_t']`

### Data access patterns

```python
# Same gene has same intrinsic kinetic rates across all cells
gene_alpha = adata.var['fit_alpha']['Ins1']  # α_Ins1 (same for all cells)
gene_beta = adata.var['fit_beta']['Ins1']    # β_Ins1 (same for all cells)  
gene_gamma = adata.var['fit_gamma']['Ins1']  # γ_Ins1 (same for all cells)

# But latent time varies per cell for the same gene
ins1_latent_times = adata.layers['fit_t'][:, gene_idx]  # t_c,Ins1 (different per cell)
# Cell 1 might be at t=0.2, Cell 2 at t=0.8 for the same gene
```

