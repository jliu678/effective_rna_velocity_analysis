---
title: "Scvelo Fitting Based On Neighbors "
date: 2025-07-10
draft: True
---

# Steps in scVelo that Use Neighbors

## 1. scv.pp.moments(adata, ...)

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

The  `use_raw`  parameter of  `scvelo.tl.recover_dynamics`  controls whether to use neighbors-smoothed moments or raw data as shown in:  `use_raw (bool or None (default: None)) – if True, use .layers['spliced'], else use moments from .layers['Ms']` 

## 2. scv.tl.velocity_graph(adata)

**Uses:**  `.obsp['connectivities']`  or recomputes internally if not present (deprecated).

**Purpose:** Computes a transition probability matrix (cosine similarity of velocity vectors) between neighboring cells.

The velocity graph computation involves:

{{< math >}} 
$$V_{ij} = \frac{\vec{v}_i \cdot (\vec{x}_j - \vec{x}_i)}{|\vec{v}_i| \cdot |\vec{x}_j - \vec{x}_i|} \tag{3}$$
{{< /math >}}

where {{< math >}} $V_{ij}$ {{< /math >}} represents the transition probability from cell {{< math >}} $i$ {{< /math >}} to cell {{< math >}} $j$ {{< /math >}}, {{< math >}} $\vec{v}_i$ {{< /math >}} is the velocity vector of cell {{< math >}} $i$ {{< /math >}}, and {{< math >}} $\vec{x}_j - \vec{x}_i$ {{< /math >}} is the displacement vector between cells.

**Depends on:** cell-cell proximity and velocity vectors – neighbor graph ensures computational efficiency and biological relevance.

## 3. scv.tl.latent_time(adata)

**Uses:** Neighbor-smoothed velocity graph to determine pseudotime-like ordering.

**Purpose:** Computes the global latent time based on inferred dynamics and the transition probabilities (again, neighbor-dependent).

The latent time computation can be expressed as:

{{< math >}} 
$$t_{\text{latent}}^{(i)} = \sum_{j} V_{ij} \cdot t_j \tag{4}$$
{{< /math >}}

where {{< math >}} $t_{\text{latent}}^{(i)}$ {{< /math >}} is the latent time for cell {{< math >}} $i$ {{< /math >}}, and the sum is over neighboring cells weighted by their transition probabilities {{< math >}} $V_{ij}$ {{< /math >}}.

## 4. scv.pl.velocity_embedding_stream(...)

**Indirectly uses:**  `velocity_graph` , which depends on the neighbor graph.

**Purpose:** Projects streamlines based on velocity vectors, assuming a local neighborhood flow.

The streamline projection follows the differential equation:

{{< math >}} 
$$\frac{d\vec{x}}{dt} = \vec{v}(\vec{x}) \tag{5}$$
{{< /math >}}

where {{< math >}} $\vec{x}$ {{< /math >}} represents the position in embedding space and {{< math >}} $\vec{v}(\vec{x})$ {{< /math >}} is the interpolated velocity field based on the neighbor-dependent velocity graph.