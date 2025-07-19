---
title: 🧬 Dynamic RNA velocity model-- (3) post hoc velocity graph 
summary: In this third blog on effectively applying the dynamic model of RNA velocity, we look into post hoc computed cosine similarity and the exponential kernel that shape the RNA velocity graph and embedding. This begins our deep dive into scVelo’s post hoc computations that determine visualization and interpretation.
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, cosine similarity, velocity graph
  - Dynamic model
  - scVelo
  - exponential kernel
  - embedding
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Introduction

In this third installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq, we start our deep dive into the post hoc computations in scVelo that shape both the visualization and interpretation of RNA velocity.

Here specifically looks into the two key components that are computed post hoc in scVelo to derive the velocity graph:
- cosine similarity in section A)
- exponential kernel transformation in section B)

And how velocity graph are projected onto an embedding, such as UMAP or t-SNE, which is the common way to visualize the inferred RNA velocity in section C)

Finally, the reconstructability score {{< math >}} $r$ {{< /math >}} in section D), which quantifies how well a subset of genes can recapitulate the overall RNA velocity dynamics inferred from the full set of genes.

## A) Cosine Similarity to compute velocity graph

### **Example Setup**
We have two cells \( i \) and \( j \), with gene expression vectors:

{{< math >}} 
$$ s_i = \begin{bmatrix} 2.0 \\ 3.5 \\ 1.0 \end{bmatrix}, \quad s_j = \begin{bmatrix} 3.5 \\ 2.5 \\ 1.5 \end{bmatrix} $$ 
{{< /math >}}

We also have the velocity vector for cell \( i \):

{{< math >}} 
$$ v_i = \begin{bmatrix} 1.2 \\ -0.5 \\ 0.3 \end{bmatrix} $$ 
{{< /math >}}

First, we compute the difference between the expression values of the two cells:

{{< math >}} 
$$ \delta_{ij} = s_j - s_i = \begin{bmatrix} 3.5 - 2.0 \\ 2.5 - 3.5 \\ 1.5 - 1.0 \end{bmatrix} = \begin{bmatrix} 1.5 \\ -1.0 \\ 0.5 \end{bmatrix} $$ 
{{< /math >}}

---

### **1. Standard Cosine Similarity**
The standard cosine similarity is given by:

{{< math >}} 
$$ \pi_{ij} = \cos(\delta_{ij}, v_i) = \frac{\delta_{ij}^T v_i}{\|\delta_{ij}\| \|v_i\|} $$ 
{{< /math >}}

**Step 1: Compute the dot product**
{{< math >}} 
$$ \delta_{ij}^T v_i = (1.5)(1.2) + (-1.0)(-0.5) + (0.5)(0.3) = 1.8 + 0.5 + 0.15 = 2.45 $$ 
{{< /math >}}

**Step 2: Compute the magnitudes**
{{< math >}} 
$$ \|\delta_{ij}\| = \sqrt{(1.5)^2 + (-1.0)^2 + (0.5)^2} = \sqrt{2.25 + 1.00 + 0.25} = \sqrt{3.5} \approx 1.87 $$ 
{{< /math >}}

{{< math >}} 
$$ \|v_i\| = \sqrt{(1.2)^2 + (-0.5)^2 + (0.3)^2} = \sqrt{1.44 + 0.25 + 0.09} = \sqrt{1.78} \approx 1.34 $$ 
{{< /math >}}

**Step 3: Compute cosine similarity**
{{< math >}} 
$$ \pi_{ij} = \frac{2.45}{(1.87)(1.34)} = \frac{2.45}{2.51} \approx 0.98 $$ 
{{< /math >}}

The vectors are strongly aligned in direction.

---

### **2. Variance-Stabilized Transformation**
Instead of using raw values, we apply:

{{< math >}} 
$$ \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|}, \quad \text{sign}(v_i) \sqrt{|v_i|} $$ 
{{< /math >}}

**Step 1: Transform \( \delta_{ij} \)**

{{< math >}} 
$$ \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|} = \begin{bmatrix} \text{sign}(1.5) \sqrt{1.5} \\ \text{sign}(-1.0) \sqrt{1.0} \\ \text{sign}(0.5) \sqrt{0.5} \end{bmatrix} = \begin{bmatrix} 1.22 \\ -1.00 \\ 0.71 \end{bmatrix} $$ 
{{< /math >}}

**Step 2: Transform \( v_i \)**

{{< math >}} 
$$ \text{sign}(v_i) \sqrt{|v_i|} = \begin{bmatrix} \text{sign}(1.2) \sqrt{1.2} \\ \text{sign}(-0.5) \sqrt{0.5} \\ \text{sign}(0.3) \sqrt{0.3} \end{bmatrix} = \begin{bmatrix} 1.10 \\ -0.71 \\ 0.55 \end{bmatrix} $$ 
{{< /math >}}

**Step 3: Compute transformed cosine similarity**
{{< math >}} 
$$ \pi_{ij}^{\text{stabilized}} = \frac{\left( 1.22, -1.00, 0.71 \right) \cdot \left( 1.10, -0.71, 0.55 \right)}{\| \text{sign}(\delta_{ij}) \sqrt{|\delta_{ij}|}\| \|\text{sign}(v_i) \sqrt{|v_i|}\|} $$ 
{{< /math >}}

**Dot product**
{{< math >}} 
$$ 1.22(1.10) + (-1.00)(-0.71) + (0.71)(0.55) = 1.342 + 0.71 + 0.391 = 2.443 $$ 
{{< /math >}}

**Compute magnitudes**
{{< math >}} 
$$ \sqrt{(1.22)^2 + (-1.00)^2 + (0.71)^2} = \sqrt{1.488 + 1.000 + 0.504} = \sqrt{2.992} \approx 1.73 $$ 
{{< /math >}}

{{< math >}} 
$$ \sqrt{(1.10)^2 + (-0.71)^2 + (0.55)^2} = \sqrt{1.210 + 0.504 + 0.302} = \sqrt{2.016} \approx 1.42 $$ 
{{< /math >}}

**Final computation**
{{< math >}} 
$$ \pi_{ij}^{\text{stabilized}} = \frac{2.443}{(1.73)(1.42)} = \frac{2.443}{2.46} \approx 0.99 $$ 
{{< /math >}}

---

## B) Exponential Kernel Transformation of Cosine Similarity

The raw cosine similarity \( \pi_{ij} \) measures directional alignment between velocity vectors, but to obtain meaningful **transition probabilities**, we apply an **exponential kernel**:

{{< math >}} 
$$ \tilde{\pi}_{ij} = \frac{1}{z_i} \exp \left(\frac{\pi_{ij}}{\sigma^2} \right) $$ 
{{< /math >}}

where:
- {{< math >}} $ \pi_{ij} $ {{< /math >}} is the **cosine similarity** between cell {{< math >}} $ i $ {{< /math >}} and {{< math >}} $ j $ {{< /math >}}.
- {{< math >}} $ \sigma $ {{< /math >}} is the **kernel width parameter**, which controls the spread of probabilities.
- {{< math >}} $ z_i $ {{< /math >}} is a **row normalization factor**, ensuring probabilities sum to 1:

{{< math >}} 
$$ z_i = \sum_j \exp \left(\frac{\pi_{ij}}{\sigma^2} \right) $$ 
{{< /math >}}

---

### **1. Why Exponential Transformation**
The **exponential function** ensures:
- **Amplification of strong correlations**: Large values of {{< math >}} $ \pi_{ij} $ {{< /math >}} get boosted.  
- **Suppression of weak correlations**: Small or negative values decay rapidly, reducing their transition influence.  
- **Nonlinear scaling**: Instead of treating all similarity scores equally, this transformation **sharpens distinctions** between stronger and weaker connections.  

Think of it as a **softmax-like weighting** that enhances **directional flow probability** based on cosine similarity.

---

### **2. Row Normalization Ensures Probabilities Sum to 1**
The **normalization factor** {{< math >}} $ z_i $ {{< /math >}} guarantees valid probability distributions:

{{< math >}} 
$$ \sum_j \tilde{\pi}_{ij} = 1 $$ 
{{< /math >}}

Without this step, raw exponentiated values might grow arbitrarily large, leading to **improper probability distributions**.

---

### **3. Role of Kernel Width σ**
- **Large {{< math >}} $ \sigma $ {{< /math >}}** → Softens probability differences, making transitions smoother.
- **Small {{< math >}} $ \sigma $ {{< /math >}}** → Sharper transitions, favoring **strong** directional alignments.

Choosing an appropriate {{< math >}} $ \sigma $ {{< /math >}} ensures **biological interpretability** in applications like RNA velocity and cell fate prediction.

## C) Project RNA Velocities into an Embedding (e.g., UMAP)

In scVelo, you have:

- **High-dimensional RNA velocities** (in gene expression space)
- **Low-dimensional embeddings** (e.g., UMAP, t-SNE)

The challenge is to project velocities into the embedding so that arrows in UMAP space reflect true dynamics.

### **1. Key Variables**
Let
- {{< math >}} $\vec{s}_i$ {{< /math >}} be the embedding position of cell {{< math >}} $i$ {{< /math >}} (e.g., in 2D UMAP space)

- {{< math >}} $\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|}$ {{< /math >}} be the normalized direction vector from cell {{< math >}} $i$ {{< /math >}} to cell {{< math >}} $j$ {{< /math >}}

- {{< math >}} $\pi_{ij}$ {{< /math >}} be the transition probability from cell {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} in the velocity-derived transition matrix

- {{< math >}} $\vec{\nu}_i$ {{< /math >}} be the velocity vector in the embedding space for cell {{< math >}} $i$ {{< /math >}}

### **2. Formula: Projected Velocity Vector**

The projected velocity for cell {{< math >}} $i$ {{< /math >}} in the embedding space is:

{{< math >}} $$\vec{\nu}_i = E_{\vec{\pi}_i}[\vec{\delta}_{ij}] = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij} - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

Or in compact form:

{{< math >}} $$\vec{\nu}_i = \vec{\delta}_{i \cdot}^{\top} \vec{\pi}_i - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

### **3. Intuition**

The **first term** is the expected movement direction in embedding space, weighted by the transition probabilities derived from RNA velocity.

The **second term**, {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}}, corrects for non-uniform density in the embedding. Without it, areas with more cells would bias the velocity field.

This ensures that the resulting velocity arrows:

- Follow the dynamics inferred from spliced/unspliced RNA
- Are not artifacts of cell density in UMAP
- Represent true directionality in cellular state transitions

### **4. In Practice**

This velocity vector {{< math >}} $\vec{\nu}_i$ {{< /math >}} is:

- Computed for each cell
- Overlaid as arrows on UMAP or t-SNE plots using `scv.pl.velocity_embedding()`

### **5. Summary**

| Term | Meaning |
|------|---------|
| {{< math >}} $\vec{s}_i$ {{< /math >}} | Position of cell {{< math >}} $i$ {{< /math >}} in embedding (e.g., UMAP) |
| {{< math >}} $\vec{\delta}_{ij}$ {{< /math >}} | Normalized direction vector from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\pi_{ij}$ {{< /math >}} | Velocity-based transition probability from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\vec{\nu}_i$ {{< /math >}} | Estimated velocity vector in embedding for cell {{< math >}} $i$ {{< /math >}} |
| {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}} | Density correction term |

### **6. Mathematical Breakdown**

#### Step 1: Compute Direction Vectors
For each cell {{< math >}} $i$ {{< /math >}} and all its neighbors {{< math >}} $j$ {{< /math >}}:

{{< math >}} $$\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|_2}$$ {{< /math >}}

#### Step 2: Get Transition Probabilities
From the velocity graph computed by scVelo:

{{< math >}} $$\pi_{ij} = P(\text{cell } i \rightarrow \text{cell } j | \text{RNA velocity})$$ {{< /math >}}

#### Step 3: Weighted Average Direction
{{< math >}} $$\text{Expected direction} = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij}$$ {{< /math >}}

#### Step 4: Density Correction
{{< math >}} $$\text{Uniform correction} = \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

#### Step 5: Final Velocity
{{< math >}} $$\vec{\nu}_i = \text{Expected direction} - \text{Uniform correction}$$ {{< /math >}}

### **7. Alternative Formulations**

#### Kernel-Based Approach
Some implementations use a kernel-weighted version:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j K(\vec{s}_i, \vec{s}_j) \pi_{ij} \vec{\delta}_{ij}}{\sum_j K(\vec{s}_i, \vec{s}_j)}$$ {{< /math >}}

where {{< math >}} $K(\vec{s}_i, \vec{s}_j)$ {{< /math >}} is a distance-based kernel (e.g., Gaussian).

#### Confidence Weighting
Incorporate velocity confidence scores:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j c_{ij} \pi_{ij} \vec{\delta}_{ij}}{\sum_j c_{ij}}$$ {{< /math >}}

where {{< math >}} $c_{ij}$ {{< /math >}} represents the confidence in the transition from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}}.

### **8. Implementation Notes**

The projected velocities are typically stored in:
- `adata.obsm['velocity_umap']` for UMAP projections
- `adata.obsm['velocity_tsne']` for t-SNE projections

These can be accessed and visualized using scVelo's plotting functions to show the directional flow of cellular development in low-dimensional space.

## D) The Reconstructability Score r

### 1. Mathematical Definition

{{< math >}} $$r = \text{median}_i \, \text{corr}(\pi_i, \pi_i')$$ {{< /math >}}

Where:

- {{< math >}} $\pi_i$ {{< /math >}}: The row vector of the full velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities to all other cells.

- {{< math >}} $\pi_i'$ {{< /math >}}: The row vector of the reduced velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities based only on the selected genes.

- {{< math >}} $\text{corr}(\pi_i, \pi_i')$ {{< /math >}}: This calculates the correlation between the full outgoing transition probabilities for cell {{< math >}} $i$ {{< /math >}} and the outgoing transition probabilities derived from the subset of genes for cell {{< math >}} $i$ {{< /math >}}. A high correlation means that the subset of genes largely reproduces the directional information provided by all genes for that particular cell.

- {{< math >}} $\text{median}_i$ {{< /math >}}: The score is then computed as the median of these per-cell correlations across all cells {{< math >}} $i$ {{< /math >}}. Taking the median makes the score robust to outliers (e.g., a few cells where the subset of genes might perform poorly due to local noise).

### 2. Meaning of Reconstructability Score

The reconstructability score quantifies how well a specific subset of genes can "reconstruct" or recapitulate the overall RNA velocity dynamics inferred from the full set of genes.

#### High Score (close to 1)
If {{< math >}} $r$ {{< /math >}} is high (e.g., {{< math >}} $> 0.8$ {{< /math >}}), it suggests that the selected subset of genes (e.g., top likelihood genes) are indeed the primary drivers of the observed cellular transitions. Their dynamics alone are largely sufficient to explain the overall directionality and speed of differentiation. This validates their "driver" status and indicates that focusing on these genes provides a good summary of the system's dynamics.

#### Low Score
If {{< math >}} $r$ {{< /math >}} is low, it means that the selected subset of genes does not adequately capture the full dynamic picture. The overall cellular transitions are likely influenced by a broader set of genes not included in your selection, or the selected genes might not be truly representative of the overall dynamics.

### 3. Mathematical Breakdown

#### Transition Probability Vectors

For each cell {{< math >}} $i$ {{< /math >}}, the transition probabilities are derived from the velocity graph:

{{< math >}} $$\pi_i = [\pi_{i,1}, \pi_{i,2}, \ldots, \pi_{i,N}]$$ {{< /math >}}

where {{< math >}} $\pi_{i,j}$ {{< /math >}} represents the probability of cell {{< math >}} $i$ {{< /math >}} transitioning to cell {{< math >}} $j$ {{< /math >}}, and {{< math >}} $N$ {{< /math >}} is the total number of cells.

#### Correlation Calculation

The correlation between full and reduced transition probabilities:

{{< math >}} $$\text{corr}(\pi_i, \pi_i') = \frac{\text{cov}(\pi_i, \pi_i')}{\sigma_{\pi_i} \sigma_{\pi_i'}}$$ {{< /math >}}

where:
- {{< math >}} $\text{cov}(\pi_i, \pi_i')$ {{< /math >}} is the covariance
- {{< math >}} $\sigma_{\pi_i}$ {{< /math >}} and {{< math >}} $\sigma_{\pi_i'}$ {{< /math >}} are the standard deviations

#### Median Aggregation

The final score uses median aggregation for robustness:

{{< math >}} $$r = \text{median}\{c_1, c_2, \ldots, c_N\}$$ {{< /math >}}

where {{< math >}} $c_i = \text{corr}(\pi_i, \pi_i')$ {{< /math >}} for each cell {{< math >}} $i$ {{< /math >}}.

### 4. Usage and Biology

#### Validation of Driver Gene Selection
It's commonly used to validate that the "top likelihood genes" (which scVelo identifies as dynamically strong) are indeed the key players shaping the global trajectory.

#### Identifying Essential Gene Sets
You can test specific gene modules or pathways to see if they are collectively responsible for a particular fate decision or developmental progression.

#### Dimensionality Reduction/Summarization
If a small subset of genes yields a high reconstructability score, it suggests that these genes form a powerful, low-dimensional representation of the system's dynamics.

#### Biological Interpretation
High reconstructability from a specific gene set points to their significant biological role in driving the observed cellular changes.