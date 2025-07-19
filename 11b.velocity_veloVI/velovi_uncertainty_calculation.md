---
title: "Velovi Uncertainty Calculation"
date: 2025-07-10
draft: True
---

# VeloVI Intrinsic Uncertainty Quantification

This document explains how VeloVI quantifies the intrinsic uncertainty in the predicted velocity for a given cell. This uncertainty reflects how consistent the model's velocity predictions are for a cell, given the inherent variability captured by the variational posterior.

Let's break down the concepts and the formula step-by-step:

## 1. Posterior Predictive Velocity Mean ({{< math >}}$\bar{v}_n${{< /math >}})

**What it is:** {{< math >}}$\bar{v}_n${{< /math >}} is the average or mean velocity vector for cell n, taking into account the model's uncertainty in the underlying latent variables ({{< math >}}$z_n${{< /math >}}).

**How it's computed (conceptual):**

We repeatedly sample {{< math >}}$z_n${{< /math >}} from its approximate posterior {{< math >}}$q_\phi(z_n | u_n, s_n)${{< /math >}}. Let's say we draw L such samples: {{< math >}}$\{z_n^{(1)}, z_n^{(2)}, \ldots, z_n^{(L)}\}${{< /math >}}.

For each sampled {{< math >}}$z_n^{(l)}${{< /math >}}:

- We compute the expected velocity for every gene g in that cell, given {{< math >}}$z_n^{(l)}${{< /math >}}. This involves calculating {{< math >}}$E_{q_\phi(\pi_{ng} | z_n^{(l)})}[v^{(g)}(t(k_{ng})_{ng}, k_{ng})]${{< /math >}} for each gene g. (As explained in our previous discussion, this results in a scalar for each gene g).

- If we put these expected gene velocities together for all G genes, we get a velocity vector sample for cell n: {{< math >}}$v_n^{(l)} = (E_{q_\phi(\pi_{n1} | z_n^{(l)})}[v^{(1)}(\ldots)], \ldots, E_{q_\phi(\pi_{nG} | z_n^{(l)})}[v^{(G)}(\ldots)])${{< /math >}}.

After obtaining L such velocity vector samples, {{< math >}}$\{v_n^{(1)}, v_n^{(2)}, \ldots, v_n^{(L)}\}${{< /math >}}, the posterior predictive velocity mean {{< math >}}$\bar{v}_n${{< /math >}} is simply their average:

{{< math >}}
$$\bar{v}_n = \frac{1}{L} \sum_{l=1}^{L} v_n^{(l)} \quad (1)$$
{{< /math >}}

This {{< math >}}$\bar{v}_n${{< /math >}} is the model's best guess for the overall velocity vector of cell n.

## 2. Intrinsic Uncertainty: The Variance of Cosine Similarity

The goal is to quantify how much the individual velocity samples {{< math >}}$v_n^{(l)}${{< /math >}} vary around this mean {{< math >}}$\bar{v}_n${{< /math >}}. However, simply calculating the variance of the vectors themselves might be misleading, as magnitude differences could dominate. VeloVI is often interested in the direction of cellular change.

**Cosine Similarity** {{< math >}}$c(v_1, v_2)${{< /math >}}: This measures the angular similarity between two vectors, ranging from -1 (opposite directions) to 1 (same direction). It's defined as:

{{< math >}}
$$c(v_1, v_2) = \frac{v_1 \cdot v_2}{\|v_1\| \|v_2\|} \quad (2)$$
{{< /math >}}

where {{< math >}}$\cdot${{< /math >}} is the dot product and {{< math >}}$\|\cdot\|${{< /math >}} is the Euclidean norm (magnitude) of the vector.

**Why Cosine Similarity for Uncertainty?**

- It focuses on the direction of the velocity, which is often biologically more interpretable than its exact magnitude for trajectories.
- It helps quantify how consistently the model predicts the direction of differentiation or change for a cell.

**The Concept:** The intrinsic uncertainty for cell n is computed as the variance of the cosine similarity between each velocity vector sample {{< math >}}$v_n^{(l)}${{< /math >}} and the overall mean velocity vector {{< math >}}$\bar{v}_n${{< /math >}}.

{{< math >}}
$$\text{Intrinsic Uncertainty} = \text{Var}_{q_\phi(v_n | u_n, s_n)}[c(v_n, \bar{v}_n)] \quad (3)$$
{{< /math >}}

This means we're looking at how spread out the cosine similarity values are. If all {{< math >}}$v_n^{(l)}${{< /math >}} point in roughly the same direction as {{< math >}}$\bar{v}_n${{< /math >}}, the cosine similarities will be close to 1, and their variance will be low. If they point in wildly different directions, the cosine similarities will vary a lot, leading to high variance.

## 3. The Formula for {{< math >}}$\hat{\sigma}_n^2${{< /math >}} (Sample Variance)

The provided formula {{< math >}}$\hat{\sigma}_n^2${{< /math >}} is the standard sample variance applied to the set of scalar cosine similarity values.

Let {{< math >}}$C_l = c(v_n^{(l)}, \bar{v}_n) = \frac{v_n^{(l)} \cdot \bar{v}_n}{\|v_n^{(l)}\| \|\bar{v}_n\|}${{< /math >}} be the cosine similarity for the l-th velocity sample.

The formula for {{< math >}}$\hat{\sigma}_n^2${{< /math >}} is:

{{< math >}}
$$\hat{\sigma}_n^2 = \frac{1}{L-1} \sum_{l=1}^{L} \left(C_l - \left(\frac{1}{L} \sum_{j=1}^{L} C_j\right)\right)^2 \quad (4)$$
{{< /math >}}

Let's break down each part of the formula:

- {{< math >}}$C_l = \frac{v_n^{(l)} \cdot \bar{v}_n}{\|v_n^{(l)}\| \|\bar{v}_n\|}${{< /math >}}: This is the cosine similarity between the l-th sampled velocity vector ({{< math >}}$v_n^{(l)}${{< /math >}}) and the overall mean velocity vector ({{< math >}}$\bar{v}_n${{< /math >}}). This value is a scalar between -1 and 1.

- {{< math >}}$\frac{1}{L} \sum_{j=1}^{L} C_j${{< /math >}}: This is the mean of all L cosine similarity values. It represents the average directional alignment of the samples with the mean. Let's call this {{< math >}}$\bar{C}${{< /math >}}.

- {{< math >}}$(C_l - \bar{C})^2${{< /math >}}: This calculates the squared difference between each individual cosine similarity value ({{< math >}}$C_l${{< /math >}}) and the overall mean cosine similarity ({{< math >}}$\bar{C}${{< /math >}}). This measures how much each sample's direction deviates from the average predicted direction.

- {{< math >}}$\sum_{l=1}^{L} (\ldots)^2${{< /math >}}: This sums up all these squared deviations across all L samples.

- {{< math >}}$\frac{1}{L-1}${{< /math >}}: This is the scaling factor for sample variance (using L-1 instead of L to provide an unbiased estimate of the population variance).

## Interpretation of {{< math >}}$\hat{\sigma}_n^2${{< /math >}}:

**High {{< math >}}$\hat{\sigma}_n^2${{< /math >}}:** Indicates high intrinsic uncertainty. The individual velocity samples {{< math >}}$v_n^{(l)}${{< /math >}} (derived from different {{< math >}}$z_n${{< /math >}} samples) are pointing in significantly different directions relative to their average direction ({{< math >}}$\bar{v}_n${{< /math >}}). This means the model is "unsure" about the precise direction of velocity for cell n, even after averaging over gene-specific kinetic states. This might indicate that the cell is at a branching point in its trajectory or that its latent representation is ambiguous.

**Low {{< math >}}$\hat{\sigma}_n^2${{< /math >}}:** Indicates low intrinsic uncertainty. The individual velocity samples {{< math >}}$v_n^{(l)}${{< /math >}} are tightly clustered around the mean direction. This means the model is confident and consistent in its prediction of the velocity direction for cell n.

# VeloVI Future Cell State Prediction and Uncertainty Quantification

This document explains how VeloVI uses the predicted velocities to predict future cell states and, crucially, how it quantifies the uncertainty in these future state predictions. This moves beyond just the uncertainty of a single cell's velocity to the uncertainty of where a cell might move to in the future.

Let's break down the explanation:

## 1. The Cell-Cell Transition Matrix {{< math >}}$T(v_{1:N}, s_{1:N})${{< /math >}}

**Purpose:** This matrix describes the likelihood of a cell transitioning from one state to another (or from cell i to cell j). It's a core component of RNA velocity analysis.

**Inputs:**
- {{< math >}}$v_{1:N}${{< /math >}}: The velocity vectors for all N cells in the dataset. Each {{< math >}}$v_i${{< /math >}} is the velocity vector for cell i.
- {{< math >}}$s_{1:N}${{< /math >}}: The spliced RNA abundance matrix for all N cells. Each {{< math >}}$s_i${{< /math >}} is the spliced abundance vector for cell i.

**How it's computed (Conceptual, "as described previously¹¹"):**

For each cell i, VeloVI identifies its nearest neighbors in the spliced abundance space (using {{< math >}}$s_{1:N}${{< /math >}}). Let's say {{< math >}}$s_j${{< /math >}} is a neighbor of {{< math >}}$s_i${{< /math >}}.

**Displacement Vector** ({{< math >}}$\delta_{ij}${{< /math >}}): This vector represents the difference in spliced RNA abundance between cell j and cell i:

{{< math >}}
$$\delta_{ij} = s_j - s_i \quad (5)$$
{{< /math >}}

It points from cell i to cell j in the gene expression space.

**Cosine Similarity**: The key idea is to compare the direction of cell i's velocity ({{< math >}}$v_i${{< /math >}}) with the direction of the displacement to its neighbor j ({{< math >}}$\delta_{ij}${{< /math >}}):

{{< math >}}
$$\cos(\delta_{ij}, v_i) = \frac{\delta_{ij}^T v_i}{\|\delta_{ij}\| \|v_i\|} \quad (6)$$
{{< /math >}}

- A high positive cosine similarity means that cell i's velocity is pointing towards cell j.
- A low or negative cosine similarity means cell i's velocity is pointing away from cell j (or in a very different direction).

**Transition Probabilities:** These cosine similarities are then typically transformed into non-negative, normalized probabilities. A common approach is to use an exponential kernel (e.g., softmax-like) on positive cosine similarities to determine the probability that cell i will transition towards cell j. This forms the entries of the transition matrix T.

The matrix T will have dimensions N×N, where {{< math >}}$T_{ij}${{< /math >}} represents the probability of transitioning from cell i to cell j. (Or sometimes {{< math >}}$T_{ji}${{< /math >}} depending on convention).

## 2. Predicted Future Cell State: {{< math >}}$T(v_{1:N}, s_{1:N})S${{< /math >}}

**Concept:** Once we have the transition matrix T, we can use it to predict the future state of the entire cellular population.

**Operation:** The matrix multiplication {{< math >}}$T(v_{1:N}, s_{1:N})S${{< /math >}} works as follows:

- S is the N×G matrix of spliced RNA abundances (cells by genes).
- The product TS (where T is N×N) will result in a new N×G matrix.
- The i-th row of this new matrix will represent the predicted future spliced RNA abundance vector for cell i. This is essentially a weighted average of the current spliced abundances of other cells, with the weights determined by the transition probabilities from cell i (i.e., the i-th row of T). If cell i is predicted to move towards cell j, its future state will be influenced more by {{< math >}}$s_j${{< /math >}}.

**Result:** This operation gives us a set of predicted future cell state vectors, one for each cell.

## 3. Quantifying Uncertainty in Predicted Future Cell States

The crucial part here is that the velocities {{< math >}}$v_{1:N}${{< /math >}} that go into computing T are themselves samples from a posterior predictive velocity distribution.

**One "sample of velocity":** The first step in the overall procedure (from the previous explanation) is to:
1. Sample {{< math >}}$z_n${{< /math >}} from {{< math >}}$q_\phi(z_n|u_n,s_n)${{< /math >}}. If we do this for all N cells, we get a specific set of {{< math >}}$z_n${{< /math >}} values {{< math >}}$\{z_1^{(l)}, \ldots, z_N^{(l)}\}${{< /math >}}.

**One set of velocity vectors:** Using these {{< math >}}$z_n^{(l)}${{< /math >}} values, and then applying the second step:
2. Compute {{< math >}}$E_{q_\phi(\pi_{ng}|z_n)}[v^{(g)}(t(k_{ng})_{ng}, k_{ng})]${{< /math >}} for each gene, we obtain a single set of velocity vectors {{< math >}}$v^{(l)} = \{v_1^{(l)}, \ldots, v_N^{(l)}\}${{< /math >}} for the entire dataset. (Each {{< math >}}$v_i^{(l)}${{< /math >}} is the expected velocity vector for cell i given that specific sampled {{< math >}}$z_i^{(l)}${{< /math >}}).

**One transition matrix and one predicted future state:** This single set of velocity vectors {{< math >}}$v^{(l)}${{< /math >}} (along with {{< math >}}$s_{1:N}${{< /math >}}) allows us to compute one sample of the transition matrix, {{< math >}}$T(v^{(l)}, s_{1:N})${{< /math >}}. Multiplying this by S gives one sample of the predicted future cell state matrix, let's call it:

{{< math >}}
$$S_{\text{future}}^{(l)} = T(v^{(l)}, s_{1:N})S \quad (7)$$
{{< /math >}}

**Repeating the Process:** To understand the uncertainty, we repeat this entire procedure L times:

1. Sample {{< math >}}$z_n${{< /math >}} for all cells to get {{< math >}}$v^{(l)} = \{v_1^{(l)}, \ldots, v_N^{(l)}\}${{< /math >}}.
2. Compute {{< math >}}$T^{(l)} = T(v^{(l)}, s_{1:N})${{< /math >}}.
3. Compute {{< math >}}$S_{\text{future}}^{(l)} = T^{(l)}S${{< /math >}}.

This results in a set of L predicted future cell state matrices: {{< math >}}$\{S_{\text{future}}^{(1)}, \ldots, S_{\text{future}}^{(L)}\}${{< /math >}}.

**Variance Computation:** For each cell i, we now have L predicted future cell state vectors: {{< math >}}$\{S_{\text{future},i}^{(1)}, S_{\text{future},i}^{(2)}, \ldots, S_{\text{future},i}^{(L)}\}${{< /math >}}.

Let {{< math >}}$\bar{S}_{\text{future},i}${{< /math >}} be the mean of these L predicted future cell state vectors for cell i:

{{< math >}}
$$\bar{S}_{\text{future},i} = \frac{1}{L} \sum_{l=1}^{L} S_{\text{future},i}^{(l)} \quad (8)$$
{{< /math >}}

The uncertainty in the predicted future cell state for cell i is then computed using the same variance computation procedure as described for the intrinsic uncertainty (variance of the cosine similarity).

Specifically, for cell i, you would compute:

{{< math >}}
$$\text{Uncertainty}(S_{\text{future},i}) = \frac{1}{L-1} \sum_{l=1}^{L} \left(\cos(S_{\text{future},i}^{(l)}, \bar{S}_{\text{future},i}) - \left(\frac{1}{L} \sum_{j=1}^{L} \cos(S_{\text{future},i}^{(j)}, \bar{S}_{\text{future},i})\right)\right)^2 \quad (9)$$
{{< /math >}}

This means:

1. Calculate the mean predicted future state vector for cell i: {{< math >}}$\bar{S}_{\text{future},i} = \frac{1}{L} \sum_{l=1}^{L} S_{\text{future},i}^{(l)}${{< /math >}}.

2. For each sample l, compute the cosine similarity between the l-th predicted future state vector {{< math >}}$S_{\text{future},i}^{(l)}${{< /math >}} and the mean predicted future state vector {{< math >}}$\bar{S}_{\text{future},i}${{< /math >}}.

3. Calculate the variance of these L cosine similarity values.

## In Essence:

This process allows VeloVI to:

**Predict where cells are likely to go** in gene expression space, using the learned velocities and current gene expression.

**Quantify the confidence** (or lack thereof) in these predictions. A high variance in cosine similarity for the predicted future states of a cell indicates that the model's different "stochastic runs" (due to sampling {{< math >}}$z_n${{< /math >}}) lead to future state predictions that point in significantly different directions. This could highlight branching points, unstable states, or regions of high biological variability that the model captures as uncertainty.