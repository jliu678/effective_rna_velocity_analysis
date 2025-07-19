---
title: "Scvelo Post Fitting Calculation"
date: 2025-07-10
draft: True
---

# D) Cosine Similarity with Variance-Stabilized Transformation

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

## **1. Standard Cosine Similarity**
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

## **2. Variance-Stabilized Transformation**
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

# E) Exponential Kernel Transformation for Transition Probabilities

### **1. Exponential Transformation of Cosine Correlation**
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

### **2. Why Exponential Transformation?**
The **exponential function** ensures:
✅ **Amplification of strong correlations**: Large values of {{< math >}} $ \pi_{ij} $ {{< /math >}} get boosted.  
✅ **Suppression of weak correlations**: Small or negative values decay rapidly, reducing their transition influence.  
✅ **Nonlinear scaling**: Instead of treating all similarity scores equally, this transformation **sharpens distinctions** between stronger and weaker connections.  

Think of it as a **softmax-like weighting** that enhances **directional flow probability** based on cosine similarity.

---

### **3. Row Normalization Ensures Probabilities Sum to 1**
The **normalization factor** {{< math >}} $ z_i $ {{< /math >}} guarantees valid probability distributions:

{{< math >}} 
$$ \sum_j \tilde{\pi}_{ij} = 1 $$ 
{{< /math >}}

Without this step, raw exponentiated values might grow arbitrarily large, leading to **improper probability distributions**.

---

### **4. Role of Kernel Width \( \sigma \)**
- **Large \( \sigma \)** → Softens probability differences, making transitions smoother.
- **Small \( \sigma \)** → Sharper transitions, favoring **strong** directional alignments.

Choosing an appropriate {{< math >}} $ \sigma $ {{< /math >}} ensures **biological interpretability** in applications like RNA velocity and cell fate prediction.

Would you like a **numerical example** demonstrating how this transformation works with actual cosine similarity values? 🚀


# F) How scVelo Identifies the Starting Points of Differentiation Using Markov Chain Analysis

## **1. What Are Root Cells in Differentiation?**
In single-cell RNA velocity analysis, root cells are the **initial states** in a developmental trajectory. These are cells that:  
✅ Have **low latent time** (early progenitors).  
✅ Are **unlikely** to have transitioned from other states.  
✅ Serve as the **starting points** for differentiation pathways.  

Identifying these root cells is crucial for understanding **how differentiation originates**.

---

## **2. Stationary State Equation in RNA Velocity**
The root cells are found by solving the **stationary distribution equation** for the **transition probability matrix** {{< math >}} $ \tilde{\pi} $ {{< /math >}}:

{{< math >}} 
$$ \mu^* = \mu^* \tilde{\pi}^T $$ 
{{< /math >}}

where:
- {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} is the **transpose of the transition probability matrix** that encodes RNA velocity-driven transitions.
- {{< math >}} $ \mu^* $ {{< /math >}} represents the **stationary distribution**, which describes the **long-term probabilities** of cells staying in each state.
- The **left eigenvectors** of {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} corresponding to the **eigenvalue of 1** define these **stationary states**.

---

## **3. Why Does This Identify Root Cells?**

### **a. Probability Evolution in a Markov Process**
In a discrete-time **Markov chain**, the state of the system at time {{< math >}} $ t $ {{< /math >}} is represented by a **probability distribution** {{< math >}} $ \mu_t $ {{< /math >}}, which evolves according to the **transition probability matrix** {{< math >}} $ P $ {{< /math >}}:

{{< math >}} 
$$ \mu_{t+1} = \mu_t P $$ 
{{< /math >}}

where:
- {{< math >}} $ \mu_t $ {{< /math >}} is a row vector of probabilities over all possible states at time {{< math >}} $ t $ {{< /math >}}.
- {{< math >}} $ P $ {{< /math >}} is an {{< math >}} $ n \times n $ {{< /math >}} **transition matrix**, where {{< math >}} $ P_{ij} $ {{< /math >}} represents the probability of transitioning from state {{< math >}} $ i $ {{< /math >}} to state {{< math >}} $ j $ {{< /math >}}.

Each step updates the probability distribution, progressively altering the likelihood of being in each state.

---

### **b. Convergence to a Stable Distribution**
As the Markov chain progresses, repeated multiplications by {{< math >}} $ P $ {{< /math >}} lead the probability distribution toward a **stable equilibrium**:

{{< math >}} 
$$ \mu_{t+k} = \mu_t P^k, \quad \text{as } k \to \infty $$ 
{{< /math >}}

After many steps, the probabilities stabilize, meaning **the system reaches a final long-term behavior** where further transitions do not significantly change the probabilities.

---

### **c. Definition of the Stationary Distribution**
The **stationary distribution** {{< math >}} $ \mu^* $ {{< /math >}} is a probability vector that remains **unchanged** under transitions:

{{< math >}} 
$$ \mu^* = \mu^* P $$ 
{{< /math >}}

This implies that if a system starts in {{< math >}} $ \mu^* $ {{< /math >}}, applying the transition matrix does **not** alter its probabilities—it stays in equilibrium.

---

### **d. Why Is the Stationary Distribution Important?**
✅ **Predicts long-term state occupancy**—the fraction of time spent in each state.  
✅ **Defines steady-state probabilities** in applications like RNA velocity analysis.  
✅ **Identifies stable states (e.g., progenitor/root cells in differentiation)** when applied to biological modeling.


## **e. Example: Computing Stationary States Using Eigenvectors**
### **Step 1: Define the Transition Probability Matrix**
Consider a **3-state system** with transition probabilities:

{{< math >}} 
$$ \tilde{\pi} = \begin{bmatrix} 
0.8 & 0.1 & 0.1 \\ 
0.2 & 0.7 & 0.1 \\ 
0.1 & 0.2 & 0.7 
\end{bmatrix} $$ 
{{< /math >}}

This matrix describes **how probabilities shift** between states.

---

### **Step 2: Find the Left Eigenvector for Eigenvalue 1**
We solve:

{{< math >}} 
$$ v^T \tilde{\pi}^T = v^T $$ 
{{< /math >}}

This means we need the **eigenvector corresponding to eigenvalue \( \lambda = 1 \)**.

Computing the eigenvectors of {{< math >}} $ \tilde{\pi}^T $ {{< /math >}}, we obtain:

{{< math >}} 
$$ v = \begin{bmatrix} 0.45 \\ 0.35 \\ 0.20 \end{bmatrix} $$ 
{{< /math >}}

This **left eigenvector** represents the **stationary distribution**, meaning:
✅ **State 1 holds 45% of the long-term probability**.  
✅ **State 2 holds 35%**.  
✅ **State 3 holds 20%**.  

These **steady-state values** describe **how the system behaves in the long run**, identifying the states **least likely to transition further**.

---

### **4. Biological Interpretation in RNA Velocity**
In **single-cell analysis**:
✅ **Cells corresponding to the stationary distribution** act as **root cells**, initiating differentiation.  
✅ This eigenvector-driven method ensures an **objective, data-driven way** to define **progenitor states**.  
✅ **Eigenvalue 1 captures states that remain stable over differentiation**, meaning they are **starting points** in biological development.



# G) Is p in eq(16) dependent on gene or cell?
In scVelo's global time normalization, we have:

{{< math >}} $$t_{i;o} = Q^p_g \left( t_{i;g} - t_{o;g} \right)$$ {{< /math >}}

Where:
- {{< math >}} $t_{i;o}$ {{< /math >}}: normalized time for cell {{< math >}} $i$ {{< /math >}} relative to root cell {{< math >}} $o$ {{< /math >}}
- {{< math >}} $Q^p_g$ {{< /math >}}: {{< math >}} $p$ {{< /math >}}-quantile computed for gene {{< math >}} $g$ {{< /math >}}
- {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}: time shift between cell {{< math >}} $i$ {{< /math >}} and root cell {{< math >}} $o$ {{< /math >}} for gene {{< math >}} $g$ {{< /math >}}

## Is p Gene-Specific?

**❌ No** — {{< math >}} $p$ {{< /math >}} is **not gene-specific**.

### Explanation

{{< math >}} $Q^p_g$ {{< /math >}} means: for **each gene** {{< math >}} $g$ {{< /math >}}, you compute the **{{< math >}} $p$ {{< /math >}}-quantile** over the time shifts {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}, across all cells {{< math >}} $i$ {{< /math >}}.

But the value of **{{< math >}} $p$ {{< /math >}}** itself is the **same for all genes** — it is **not gene-specific**.

So even though the quantile is applied **per gene** (i.e., each gene gets its own {{< math >}} $Q^p_g$ {{< /math >}}), the **{{< math >}} $p$ {{< /math >}} used in all those calculations is shared**.

### Mathematical Illustration

For a fixed value {{< math >}} $p = 0.5$ {{< /math >}} (median):

{{< math >}} $$Q^{0.5}_{gene1} = \text{median}(\{t_{i;gene1} - t_{o;gene1}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{gene2} = \text{median}(\{t_{i;gene2} - t_{o;gene2}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$\vdots$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{geneM} = \text{median}(\{t_{i;geneM} - t_{o;geneM}\}_{i=1}^N)$$ {{< /math >}}

The **same** {{< math >}} $p = 0.5$ {{< /math >}} is used for all genes, but each gene gets its own quantile value.

## Is p Cell-Specific?

**❌ Also no** — {{< math >}} $p$ {{< /math >}} is not adjusted for each cell. It is only used to **transform the gene-wise temporal shifts** into a consistent time estimate **for each cell**, **relative to the root cell** {{< math >}} $o$ {{< /math >}}.

### Explanation

The only place where indexing over cells happens is in computing the latent time {{< math >}} $t_{i;o}$ {{< /math >}}, but the **value of {{< math >}} $p$ {{< /math >}} stays fixed** for all cells.

### Mathematical Process

1. **Fixed {{< math >}} $p$ {{< /math >}} across all computations**: {{< math >}} $p = p_{\text{fixed}}$ {{< /math >}}

2. **Gene-wise quantile computation**:
   {{< math >}} $$Q^{p_{\text{fixed}}}_g = \text{quantile}_{p_{\text{fixed}}}(\{t_{i;g} - t_{o;g}\}_{i=1}^N)$$ {{< /math >}}

3. **Cell-specific time assignment**:
   {{< math >}} $$t_{1;o} = Q^{p_{\text{fixed}}}_g \left( t_{1;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$t_{2;o} = Q^{p_{\text{fixed}}}_g \left( t_{2;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$\vdots$$ {{< /math >}}
   {{< math >}} $$t_{N;o} = Q^{p_{\text{fixed}}}_g \left( t_{N;g} - t_{o;g} \right)$$ {{< /math >}}

## Is p Root-Specific?

**✅ Yes** — scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness.

### Optimization Process

For each potential root cell {{< math >}} $o$ {{< /math >}}, scVelo optimizes:

{{< math >}} $$p^*_o = \arg\min_p \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

Where the loss function might be:

{{< math >}} $$\mathcal{L}_{\text{smoothness}}(p, o) = \sum_{i,j \in \text{neighbors}} \left| t_{i;o}(p) - t_{j;o}(p) \right|^2$$ {{< /math >}}

### Root Selection

The final root and {{< math >}} $p$ {{< /math >}} are chosen as:

{{< math >}} $$o^*, p^* = \arg\min_{o,p} \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

## Summary Table

| Aspect | Is {{< math >}} $p$ {{< /math >}} dependent? | Explanation |
|--------|-----|-------------|
| **Per gene** {{< math >}} $g$ {{< /math >}} | ❌ No | The quantile {{< math >}} $Q^p$ {{< /math >}} is computed per gene, but **{{< math >}} $p$ {{< /math >}} itself is fixed** |
| **Per cell** {{< math >}} $i$ {{< /math >}} | ❌ No | All cells use the same {{< math >}} $p$ {{< /math >}}; only their computed values differ |
| **Per root** {{< math >}} $o$ {{< /math >}} | ✅ Yes | scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness |


## Optimization Step

Once **scVelo** computes all {{< math >}} $ t_{i,o} $ {{< /math >}} for different candidate root cells {{< math >}} $ o $ {{< /math >}}, it looks for the root and quantile {{< math >}} $ p $ {{< /math >}} that maximize correlation between the **latent time series**  
{{< math >}} 
$$ \{ t_{i,o} \} $$ 
{{< /math >}}  
and its **convolution over a local neighborhood graph** (e.g., KNN graph of cells).

### **Mathematical Formulation**
The optimization step seeks:

{{< math >}} 
$$ \arg\max_{o,p} \text{corr} \left( t_{i,o}, \sum_{j \in N(i)} w_{ij} t_{j,o} \right) $$ 
{{< /math >}}

where:
- {{< math >}} $ N(i) $ {{< /math >}} is the **neighborhood of cell** {{< math >}} $ i $ {{< /math >}}.
- {{< math >}} $ w_{ij} $ {{< /math >}} represents **graph weights**, determining connectivity between cells.



# H) Project RNA Velocities into an Embedding (e.g., UMAP)

In scVelo, you have:

- **High-dimensional RNA velocities** (in gene expression space)
- **Low-dimensional embeddings** (e.g., UMAP, t-SNE)

The challenge is to project velocities into the embedding so that arrows in UMAP space reflect true dynamics.

## 📐 Key Variables

Let {{< math >}} $\vec{s}_i$ {{< /math >}} be the embedding position of cell {{< math >}} $i$ {{< /math >}} (e.g., in 2D UMAP space)

Let {{< math >}} $\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|}$ {{< /math >}} be the normalized direction vector from cell {{< math >}} $i$ {{< /math >}} to cell {{< math >}} $j$ {{< /math >}}

Let {{< math >}} $\pi_{ij}$ {{< /math >}} be the transition probability from cell {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} in the velocity-derived transition matrix

Let {{< math >}} $\vec{\nu}_i$ {{< /math >}} be the velocity vector in the embedding space for cell {{< math >}} $i$ {{< /math >}}

## 📊 Formula: Projected Velocity Vector

The projected velocity for cell {{< math >}} $i$ {{< /math >}} in the embedding space is:

{{< math >}} $$\vec{\nu}_i = E_{\vec{\pi}_i}[\vec{\delta}_{ij}] = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij} - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

Or in compact form:

{{< math >}} $$\vec{\nu}_i = \vec{\delta}_{i \cdot}^{\top} \vec{\pi}_i - \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

## 🧠 Intuition

The **first term** is the expected movement direction in embedding space, weighted by the transition probabilities derived from RNA velocity.

The **second term**, {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}}, corrects for non-uniform density in the embedding. Without it, areas with more cells would bias the velocity field.

This ensures that the resulting velocity arrows:

- Follow the dynamics inferred from spliced/unspliced RNA
- Are not artifacts of cell density in UMAP
- Represent true directionality in cellular state transitions

## 🖼️ In Practice

This velocity vector {{< math >}} $\vec{\nu}_i$ {{< /math >}} is:

- Computed for each cell
- Overlaid as arrows on UMAP or t-SNE plots using `scv.pl.velocity_embedding()`

## ✅ Summary

| Term | Meaning |
|------|---------|
| {{< math >}} $\vec{s}_i$ {{< /math >}} | Position of cell {{< math >}} $i$ {{< /math >}} in embedding (e.g., UMAP) |
| {{< math >}} $\vec{\delta}_{ij}$ {{< /math >}} | Normalized direction vector from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\pi_{ij}$ {{< /math >}} | Velocity-based transition probability from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}} |
| {{< math >}} $\vec{\nu}_i$ {{< /math >}} | Estimated velocity vector in embedding for cell {{< math >}} $i$ {{< /math >}} |
| {{< math >}} $\frac{1}{n} \sum \vec{\delta}_{ij}$ {{< /math >}} | Density correction term |

## Mathematical Breakdown

### Step 1: Compute Direction Vectors
For each cell {{< math >}} $i$ {{< /math >}} and all its neighbors {{< math >}} $j$ {{< /math >}}:

{{< math >}} $$\vec{\delta}_{ij} = \frac{\vec{s}_j - \vec{s}_i}{\|\vec{s}_j - \vec{s}_i\|_2}$$ {{< /math >}}

### Step 2: Get Transition Probabilities
From the velocity graph computed by scVelo:

{{< math >}} $$\pi_{ij} = P(\text{cell } i \rightarrow \text{cell } j | \text{RNA velocity})$$ {{< /math >}}

### Step 3: Weighted Average Direction
{{< math >}} $$\text{Expected direction} = \sum_{j \neq i} \pi_{ij} \vec{\delta}_{ij}$$ {{< /math >}}

### Step 4: Density Correction
{{< math >}} $$\text{Uniform correction} = \frac{1}{n} \sum_{j \neq i} \vec{\delta}_{ij}$$ {{< /math >}}

### Step 5: Final Velocity
{{< math >}} $$\vec{\nu}_i = \text{Expected direction} - \text{Uniform correction}$$ {{< /math >}}

## Alternative Formulations

### Kernel-Based Approach
Some implementations use a kernel-weighted version:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j K(\vec{s}_i, \vec{s}_j) \pi_{ij} \vec{\delta}_{ij}}{\sum_j K(\vec{s}_i, \vec{s}_j)}$$ {{< /math >}}

where {{< math >}} $K(\vec{s}_i, \vec{s}_j)$ {{< /math >}} is a distance-based kernel (e.g., Gaussian).

### Confidence Weighting
Incorporate velocity confidence scores:

{{< math >}} $$\vec{\nu}_i = \frac{\sum_j c_{ij} \pi_{ij} \vec{\delta}_{ij}}{\sum_j c_{ij}}$$ {{< /math >}}

where {{< math >}} $c_{ij}$ {{< /math >}} represents the confidence in the transition from {{< math >}} $i$ {{< /math >}} to {{< math >}} $j$ {{< /math >}}.

## Implementation Notes

The projected velocities are typically stored in:
- `adata.obsm['velocity_umap']` for UMAP projections
- `adata.obsm['velocity_tsne']` for t-SNE projections

These can be accessed and visualized using scVelo's plotting functions to show the directional flow of cellular development in low-dimensional space.

# I) The Reconstructability Score (r)

## Mathematical Definition

{{< math >}} $$r = \text{median}_i \, \text{corr}(\pi_i, \pi_i')$$ {{< /math >}}

Where:

- {{< math >}} $\pi_i$ {{< /math >}}: The row vector of the full velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities to all other cells.

- {{< math >}} $\pi_i'$ {{< /math >}}: The row vector of the reduced velocity graph for cell {{< math >}} $i$ {{< /math >}}, representing its outgoing transition probabilities based only on the selected genes.

- {{< math >}} $\text{corr}(\pi_i, \pi_i')$ {{< /math >}}: This calculates the correlation between the full outgoing transition probabilities for cell {{< math >}} $i$ {{< /math >}} and the outgoing transition probabilities derived from the subset of genes for cell {{< math >}} $i$ {{< /math >}}. A high correlation means that the subset of genes largely reproduces the directional information provided by all genes for that particular cell.

- {{< math >}} $\text{median}_i$ {{< /math >}}: The score is then computed as the median of these per-cell correlations across all cells {{< math >}} $i$ {{< /math >}}. Taking the median makes the score robust to outliers (e.g., a few cells where the subset of genes might perform poorly due to local noise).

## What the Reconstructability Score Tells You

The reconstructability score quantifies how well a specific subset of genes can "reconstruct" or recapitulate the overall RNA velocity dynamics inferred from the full set of genes.

### High Score (close to 1)
If {{< math >}} $r$ {{< /math >}} is high (e.g., {{< math >}} $> 0.8$ {{< /math >}}), it suggests that the selected subset of genes (e.g., top likelihood genes) are indeed the primary drivers of the observed cellular transitions. Their dynamics alone are largely sufficient to explain the overall directionality and speed of differentiation. This validates their "driver" status and indicates that focusing on these genes provides a good summary of the system's dynamics.

### Low Score
If {{< math >}} $r$ {{< /math >}} is low, it means that the selected subset of genes does not adequately capture the full dynamic picture. The overall cellular transitions are likely influenced by a broader set of genes not included in your selection, or the selected genes might not be truly representative of the overall dynamics.

## Mathematical Components Breakdown

### Transition Probability Vectors

For each cell {{< math >}} $i$ {{< /math >}}, the transition probabilities are derived from the velocity graph:

{{< math >}} $$\pi_i = [\pi_{i,1}, \pi_{i,2}, \ldots, \pi_{i,N}]$$ {{< /math >}}

where {{< math >}} $\pi_{i,j}$ {{< /math >}} represents the probability of cell {{< math >}} $i$ {{< /math >}} transitioning to cell {{< math >}} $j$ {{< /math >}}, and {{< math >}} $N$ {{< /math >}} is the total number of cells.

### Correlation Calculation

The correlation between full and reduced transition probabilities:

{{< math >}} $$\text{corr}(\pi_i, \pi_i') = \frac{\text{cov}(\pi_i, \pi_i')}{\sigma_{\pi_i} \sigma_{\pi_i'}}$$ {{< /math >}}

where:
- {{< math >}} $\text{cov}(\pi_i, \pi_i')$ {{< /math >}} is the covariance
- {{< math >}} $\sigma_{\pi_i}$ {{< /math >}} and {{< math >}} $\sigma_{\pi_i'}$ {{< /math >}} are the standard deviations

### Median Aggregation

The final score uses median aggregation for robustness:

{{< math >}} $$r = \text{median}\{c_1, c_2, \ldots, c_N\}$$ {{< /math >}}

where {{< math >}} $c_i = \text{corr}(\pi_i, \pi_i')$ {{< /math >}} for each cell {{< math >}} $i$ {{< /math >}}.

## Use Cases and Significance

### Validation of Driver Gene Selection
It's commonly used to validate that the "top likelihood genes" (which scVelo identifies as dynamically strong) are indeed the key players shaping the global trajectory.

### Identifying Essential Gene Sets
You can test specific gene modules or pathways to see if they are collectively responsible for a particular fate decision or developmental progression.

### Dimensionality Reduction/Summarization
If a small subset of genes yields a high reconstructability score, it suggests that these genes form a powerful, low-dimensional representation of the system's dynamics.

### Biological Interpretation
High reconstructability from a specific gene set points to their significant biological role in driving the observed cellular changes.

# J) Transport Maps and Cell Distribution Analysis

## Mathematical Framework

### Original Notation

{{< math >}} $\mu = (\mu_1, \ldots, \mu_n)$ {{< /math >}}: a distribution over cells

{{< math >}} $\pi_e$ {{< /math >}}: a transport map, e.g., inferred from RNA velocity

{{< math >}} $\tilde{\pi}^T$ {{< /math >}}: the transpose (or reverse) of the transport map

{{< math >}} $\mu_{\text{des}} = \mu \cdot \pi_e$ {{< /math >}}: the descendant distribution

{{< math >}} $\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T$ {{< /math >}}: the ancestor distribution

## Core Concepts

### 1. Forward Transport: "Where are cells going?"

A distribution of cells {{< math >}} $\mu$ {{< /math >}} can be pushed through the transport map to determine future cell states.

Starting with a distribution over some cells {{< math >}} $\mu$ {{< /math >}}, for example:
- All cells in a particular cluster or cell type
- Or a uniform distribution over some selected cells

Pushing forward this distribution through the transport map {{< math >}} $\pi_e$ {{< /math >}} gives the distribution of where these cells go in the future:

{{< math >}}
$$\mu_{\text{des}} = \mu \cdot \pi_e \tag{1}$$
{{< /math >}}

This answers the question: "Given that my cells are here now, where are they going?"

### 2. Reverse Transport: "Where did cells come from?"

Conversely, a distribution {{< math >}} $\mu$ {{< /math >}} can be pulled back to determine cell origins.

You take a set of cells at a later time or later state, and ask: "Where did they come from?"

You pull back the distribution via the transpose of the transport matrix:

{{< math >}}
$$\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T \tag{2}$$
{{< /math >}}

The transpose reverses the direction of flow — it computes how much each earlier state contributed to the current cell distribution.

### 3. Indicator-Based Cell Selection

For a set of cells {{< math >}} $S = \{s_1, \ldots, s_n\}$ {{< /math >}}, the distribution {{< math >}} $\mu$ {{< /math >}} is defined using an indicator vector:

{{< math >}}
$$\mu_i = \mathbf{1}_{s_i \in S} \tag{3}$$
{{< /math >}}

This means:
- {{< math >}} $\mu_i = 1$ {{< /math >}} if cell {{< math >}} $i$ {{< /math >}} is in the set {{< math >}} $S$ {{< /math >}}
- {{< math >}} $\mu_i = 0$ {{< /math >}} otherwise

The input distribution puts mass 1 on each selected cell, and 0 elsewhere.

## Biological Interpretation

### Example Use Cases

**Progenitor Analysis:**
1. Select a set of progenitor cells (set {{< math >}} $S$ {{< /math >}})
2. Compute {{< math >}} $\mu_{\text{des}} = \mu \cdot \pi_e$ {{< /math >}} → get the expected distribution of their descendants

**Lineage Tracing:**
1. Select differentiated cells
2. Compute {{< math >}} $\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T$ {{< /math >}} → get their likely origins

This framework helps answer:
- "Where are these cells going?"
- "Where did these cells come from?"

## Numerical Example

Consider 4 cells ({{< math >}} $n=4$ {{< /math >}}), selecting cells 1 and 3:

{{< math >}}
$$\mu = [1, 0, 1, 0] \tag{4}$$
{{< /math >}}

Assume a learned transport map:

{{< math >}}
$$\pi_e = \begin{bmatrix}
0.6 & 0.2 & 0.1 & 0.1 \\
0.1 & 0.7 & 0.1 & 0.1 \\
0.2 & 0.1 & 0.6 & 0.1 \\
0.3 & 0.2 & 0.1 & 0.4
\end{bmatrix} \tag{5}$$
{{< /math >}}

Then the descendant distribution is:

{{< math >}}
$$\mu_{\text{des}} = \mu \cdot \pi_e = [1, 0, 1, 0] \cdot \pi_e = \text{row 1} + \text{row 3} \tag{6}$$
{{< /math >}}

Computing the result:

{{< math >}}
$$\mu_{\text{des}} = [0.6 + 0.2, 0.2 + 0.1, 0.1 + 0.6, 0.1 + 0.1] = [0.8, 0.3, 0.7, 0.2] \tag{7}$$
{{< /math >}}

This gives a distribution over where the selected cells (1 & 3) are headed — their descendant cell states.

## Reverse Transport Example
### Step 1: Compute the Ancestor Distribution

To find the ancestor distribution, we use the transpose of {{< math >}} $\pi_e$ {{< /math >}}:

{{< math >}} 
$$(\pi_e)^T = \begin{bmatrix}
0.6 & 0.1 & 0.2 & 0.3 \\
0.2 & 0.7 & 0.1 & 0.2 \\
0.1 & 0.1 & 0.6 & 0.1 \\
0.1 & 0.1 & 0.1 & 0.4
\end{bmatrix} \tag{3}$$
{{< /math >}}

### Step 2: Calculate the Ancestor Distribution

The ancestor distribution is computed as:

{{< math >}} 
$$\mu_{\text{anc}} = \mu_{\text{des}} \cdot (\pi_e)^T \tag{4}$$
{{< /math >}}

Multiplying the row vector {{< math >}} $\mu_{\text{des}} = [0.8, 0.3, 0.7, 0.2]$ {{< /math >}} by the matrix {{< math >}} $(\pi_e)^T$ {{< /math >}}:

{{< math >}} 
$$\mu_{\text{anc},1} = 0.8 \times 0.6 + 0.3 \times 0.1 + 0.7 \times 0.2 + 0.2 \times 0.3 = 0.48 + 0.03 + 0.14 + 0.06 = 0.71 $$
{{< /math >}}

{{< math >}} 
$$\mu_{\text{anc},2} = 0.8 \times 0.2 + 0.3 \times 0.7 + 0.7 \times 0.1 + 0.2 \times 0.2 = 0.16 + 0.21 + 0.07 + 0.04 = 0.48 $$
{{< /math >}}

{{< math >}} 
$$\mu_{\text{anc},3} = 0.8 \times 0.1 + 0.3 \times 0.1 + 0.7 \times 0.6 + 0.2 \times 0.1 = 0.08 + 0.03 + 0.42 + 0.02 = 0.55 $$
{{< /math >}}

{{< math >}} 
$$\mu_{\text{anc},4} = 0.8 \times 0.1 + 0.3 \times 0.1 + 0.7 \times 0.1 + 0.2 \times 0.4 = 0.08 + 0.03 + 0.07 + 0.08 = 0.26 $$
{{< /math >}}

### Step 3: Resulting Ancestor Distribution

{{< math >}} 
$$\mu_{\text{anc}} = [0.71, 0.48, 0.55, 0.26] \tag{9}$$
{{< /math >}}

**Interpretation**: Given that the descendant distribution was {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}}, pulling it back through the transpose of the transport map estimates the likely ancestor distribution of these descendants. The weights show which earlier states contributed more to this descendant cell population.

### Normalization

We can normalize {{< math >}} $\mu_{\text{anc}}$ {{< /math >}} to turn it into a proper probability distribution:

{{< math >}} 
$$\text{Sum} = 0.71 + 0.48 + 0.55 + 0.26 = 2.00 $$
{{< /math >}}

{{< math >}} 
$$\mu_{\text{anc,normalized}} = [0.355, 0.24, 0.275, 0.13] $$
{{< /math >}}

### Why Perfect Inversion is Impossible

The descendant distribution {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}} cannot lead exactly back to the original indicator vector {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} when using the transpose of the transport map for several fundamental reasons:

#### 1. Non-invertibility of the Transport Map

The transport matrix {{< math >}} $\pi_e$ {{< /math >}} typically represents a probabilistic, many-to-many transition between cell states. It's a stochastic matrix where rows sum to 1, but it's often not invertible or the inverse is not unique. The operation of pushing forward loses information through mixing and averaging, and pulling back is not a perfect inverse.

#### 2. Smoothing and Mixing

When we multiplied {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} by {{< math >}} $\pi_e$ {{< /math >}}, we obtained a weighted average of rows 1 and 3, resulting in the smoothed vector {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}}. This vector is continuous-valued and represents a distribution over many possible descendant states. The reverse step with the transpose spreads these probabilities back to ancestors but cannot undo the smoothing perfectly.

#### 3. Loss of Sharpness

The original vector {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} is a sharp indicator (only cells 1 and 3 are selected). After forward mapping, the distribution reflects uncertainty or spreading. Pulling back cannot perfectly re-sharpen because the transport matrix mixes information.

#### 4. Biological Analogy

In biological systems, a population of cells can differentiate into multiple descendant states probabilistically. The future cell states reflect a blend of many ancestors. Hence, you cannot simply invert descendants back to a unique set of ancestors — instead, you get a probabilistic ancestor distribution that represents the uncertainty inherent in the biological process.

# K) Calculating the gene shared latent time in scVelo

## Step-by-step explanation

### 1. **Latent time**
* scVelo computes a **latent time** {{< math >}} $t_n$ {{< /math >}} for each cell {{< math >}} $n$ {{< /math >}}, indicating **how far along a differentiation trajectory** it is, starting from the inferred root cells.
* This time is estimated based on how consistent a cell's transcriptional state is with the **direction of gene expression changes (RNA velocity)**.

### 2. **Gene-shared latent time course**
* This is the original latent time vector computed across all cells.
* It's "gene-shared" because it's not just one gene's trajectory — it aggregates information across genes to produce a **common timeline**.

### 3. **Neighborhood convolution**
* scVelo smooths the latent time values using **k-nearest neighbor (KNN)** information (e.g., based on transcriptomic similarity).
* The **convolved time** for a cell is the **average of the latent times of its neighboring cells**.

Think of it as:

{{< math >}} 
$$\tilde{t}_n = \frac{1}{|N(n)|} \sum_{m \in N(n)} t_m \tag{1}$$
{{< /math >}}

where {{< math >}} $N(n)$ {{< /math >}} is the neighborhood of cell {{< math >}} $n$ {{< /math >}}.

### 4. **Regression & consistency check**
* scVelo performs **linear regression**: it tries to predict each cell's latent time using its convolved value.
* If a cell's latent time significantly **deviates from the trend of its neighbors**, it is likely noisy or unreliable.

### 5. **Correction / Replacement**
* For these **outlier cells**, scVelo **replaces their latent time** with the convolved (smoothed) version.
* This acts like a denoising step — ensuring latent time values are **locally smooth and robust** to noise.

**scVelo uses neighborhood convolution** to detect and fix inconsistent latent times using regression.

## 🧪 Setup

Let's say we have **5 cells**, and we already computed latent times for them as:

| Cell | Latent Time {{< math >}} $t_n$ {{< /math >}} |
|------|----------------------------------------------|
| A    | 0.10                                        |
| B    | 0.15                                        |
| C    | 0.85 ❗                                     |
| D    | 0.20                                        |
| E    | 0.25                                        |

Clearly, cell **C** looks suspicious — it's way ahead in time compared to its neighbors.

## 👯 Step 1: Define neighbors

Let's say the **neighborhoods (e.g., via KNN)** are:
* A → neighbors: [B, D]
* B → neighbors: [A, D]
* C → neighbors: [B, D, E]
* D → neighbors: [A, B, E]
* E → neighbors: [C, D]

## 🔁 Step 2: Compute convolved (smoothed) latent times

{{< math >}} 
$$\tilde{t}_n = \text{average of neighbor latent times} \tag{2}$$
{{< /math >}}

| Cell | Neighbors | Neighbor Times    | Convolved Time {{< math >}} $\tilde{t}_n$ {{< /math >}} |
|------|-----------|-------------------|--------------------------------------------------------|
| A    | B, D      | 0.15, 0.20        | 0.175                                                  |
| B    | A, D      | 0.10, 0.20        | 0.15                                                   |
| C    | B, D, E   | 0.15, 0.20, 0.25  | 0.20                                                   |
| D    | A, B, E   | 0.10, 0.15, 0.25  | 0.167                                                  |
| E    | C, D      | 0.85, 0.20        | 0.525 ❗                                              |

## 📉 Step 3: Compare original vs convolved times

| Cell | {{< math >}} $t_n$ {{< /math >}} | {{< math >}} $\tilde{t}_n$ {{< /math >}} | Difference |
|------|-----------------------------------|-------------------------------------------|------------|
| A    | 0.10                             | 0.175                                     | 0.075      |
| B    | 0.15                             | 0.15                                      | 0.000      |
| C    | 0.85 ❗                          | 0.20                                      | 0.65 ❗    |
| D    | 0.20                             | 0.167                                     | 0.033      |
| E    | 0.25                             | 0.525 ❗                                  | 0.275 ❗   |

## 🧽 Step 4: Detect inconsistent cells

* scVelo runs a **regression of** {{< math >}} $t_n \sim \tilde{t}_n$ {{< /math >}} and checks for large residuals.
* Cells **C** and **E** have large mismatches between original and smoothed values → likely noisy or incorrect.

## 🔁 Step 5: Replace inconsistent values

Let's say we define a cutoff: if difference > 0.3, replace with {{< math >}} $\tilde{t}_n$ {{< /math >}}.

| Cell | New {{< math >}} $t_n$ {{< /math >}} |
|------|--------------------------------------|
| A    | 0.10 (unchanged)                    |
| B    | 0.15 (unchanged)                    |
| C    | **0.20** (replaced!)               |
| D    | 0.20 (unchanged)                    |
| E    | **0.525** (replaced!)              |

## ✅ Final corrected latent time:

{{< math >}} 
$$\mathbf{t} = [0.10, 0.15, \mathbf{0.20}, 0.20, \mathbf{0.525}] \tag{3}$$
{{< /math >}}

Cell C is now brought **back in line** with its local trajectory.