---
title: "Scvelo Fitting Based On Neighbors "
date: 2025-07-10
draft: True
---
# scVelo Functions: Random Seed Dependencies

This guide outlines which scVelo functions depend on random seeds and which are deterministic.

## Steps that **depend** on random seed

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

## Steps that **do not depend** on random seed

{{% callout note %}}
These functions are deterministic and will produce identical results given the same input data.
{{% /callout %}}

| Function | Description | Seed Dependency |
|----------|-------------|-----------------|
| `scv.pp.filter_and_normalize()` | Basic filtering and normalization | ❌ **No** - Purely deterministic operations |
| `scv.pp.moments()` | Computes first and second moments for smoothing | ❌ **No** - Deterministic given neighbors graph |
| `scv.tl.velocity()` | Computes RNA velocity using stochastic/dynamical model | ❌ **Mostly deterministic** - Depends on moments and model parameters |

## Best Practices for Reproducibility

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

*This documentation covers scVelo version compatibility and function behavior regarding random seed dependencies.*