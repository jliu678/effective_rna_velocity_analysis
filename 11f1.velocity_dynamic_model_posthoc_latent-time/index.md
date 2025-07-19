---
title: 🧬 Dynamic RNA velocity model-- (4) latent time 
summary: Here reveals the mathematical foundations of **latent time**, which enables in-depth grasp and interpretation of the latent time of RNA velocity. This is the fourth installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq. 
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, latent time, root cell
  - Transport maps
  - Neighborhood convolution
  - scVelo
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Introduction

In this fourth installment of our blog series on effectively applying the dynamic model to infer RNA velocity from single-cell RNA-seq, we reveal the mathematical foundations of **latent time** in the context of RNA velocity analysis, which enables in-depth grasp and interpretation of the latent time of RNA velocity.

Latent time is a crucial interpretation of RNA velocity independent of pre-defined (never perfect) dimentional reduction and it embedding. And latent time provide expressive visulization of the **temporal dynamics of cell differentiation** according to RNA velocity. 

Here specifically focuses on:
- Identification of root cells in section A)
- Discussion and numeric examples of transport maps in section B)
- Gene and cell dependence of latent time in section C)
- Neighborhood-convolved latent time in section D)

## A) Identify Starting Points Using Markov Chain Analysis

### **1. What Are Root Cells in Differentiation**
In single-cell RNA velocity analysis, root cells are the **initial states** in a developmental trajectory. These are cells that:  
✅ Have **low latent time** (early progenitors).  
✅ Are **unlikely** to have transitioned from other states.  
✅ Serve as the **starting points** for differentiation pathways.  

Identifying these root cells is crucial for understanding **how differentiation originates**.

---

### **2. Stationary State Equation in RNA Velocity**
The root cells are found by solving the **stationary distribution equation** for the **transition probability matrix** {{< math >}} $ \tilde{\pi} $ {{< /math >}}:

{{< math >}} 
$$ \mu^* = \mu^* \tilde{\pi}^T $$ 
{{< /math >}}

where:
- {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} is the **transpose of the transition probability matrix** that encodes RNA velocity-driven transitions.
- {{< math >}} $ \mu^* $ {{< /math >}} represents the **stationary distribution**, which describes the **long-term probabilities** of cells staying in each state.
- The **left eigenvectors** of {{< math >}} $ \tilde{\pi}^T $ {{< /math >}} corresponding to the **eigenvalue of 1** define these **stationary states**.

---

### **3. Why Does This Identify Root Cells**

#### **a. Probability Evolution in a Markov Process**
In a discrete-time **Markov chain**, the state of the system at time {{< math >}} $ t $ {{< /math >}} is represented by a **probability distribution** {{< math >}} $ \mu_t $ {{< /math >}}, which evolves according to the **transition probability matrix** {{< math >}} $ P $ {{< /math >}}:

{{< math >}} 
$$ \mu_{t+1} = \mu_t P $$ 
{{< /math >}}

where:
- {{< math >}} $ \mu_t $ {{< /math >}} is a row vector of probabilities over all possible states at time {{< math >}} $ t $ {{< /math >}}.
- {{< math >}} $ P $ {{< /math >}} is an {{< math >}} $ n \times n $ {{< /math >}} **transition matrix**, where {{< math >}} $ P_{ij} $ {{< /math >}} represents the probability of transitioning from state {{< math >}} $ i $ {{< /math >}} to state {{< math >}} $ j $ {{< /math >}}.

Each step updates the probability distribution, progressively altering the likelihood of being in each state.

---

#### **b. Convergence to a Stable Distribution**
As the Markov chain progresses, repeated multiplications by {{< math >}} $ P $ {{< /math >}} lead the probability distribution toward a **stable equilibrium**:

{{< math >}} 
$$ \mu_{t+k} = \mu_t P^k, \quad \text{as } k \to \infty $$ 
{{< /math >}}

After many steps, the probabilities stabilize, meaning **the system reaches a final long-term behavior** where further transitions do not significantly change the probabilities.

---

#### **c. Definition of the Stationary Distribution**
The **stationary distribution** {{< math >}} $ \mu^* $ {{< /math >}} is a probability vector that remains **unchanged** under transitions:

{{< math >}} 
$$ \mu^* = \mu^* P $$ 
{{< /math >}}

This implies that if a system starts in {{< math >}} $ \mu^* $ {{< /math >}}, applying the transition matrix does **not** alter its probabilities--it stays in equilibrium.

---

#### **d. Interpretation of Stationary Distribution**
✅ **Predicts long-term state occupancy**--the fraction of time spent in each state.  
✅ **Defines steady-state probabilities** in applications like RNA velocity analysis.  
✅ **Identifies stable states (e.g., progenitor/root cells in differentiation)** when applied to biological modeling.


### **4. Compute Stationary States Using Eigenvectors**
#### **Step 1: Define the Transition Probability Matrix**
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

#### **Step 2: Find the Left Eigenvector for Eigenvalue 1**
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

#### **Summery & Interpretation**
In **single-cell analysis**:
- **Cells corresponding to the stationary distribution** act as **root cells**, initiating differentiation.  
- This eigenvector-driven method ensures an **objective, data-driven way** to define **progenitor states**.  
- **Eigenvalue 1 captures states that remain stable over differentiation**, meaning they are **starting points** in biological development.

---

## B) Details and Numeric Examples of Transport Maps for Cell Development Analysis

### 1. Mathematical Notation

{{< math >}} $\mu = (\mu_1, \ldots, \mu_n)$ {{< /math >}}: a distribution over cells

{{< math >}} $\pi_e$ {{< /math >}}: a transport map, e.g., inferred from RNA velocity

{{< math >}} $\tilde{\pi}^T$ {{< /math >}}: the transpose (or reverse) of the transport map

{{< math >}} $\mu_{\text{des}} = \mu \cdot \pi_e$ {{< /math >}}: the descendant distribution

{{< math >}} $\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T$ {{< /math >}}: the ancestor distribution

### 2. Core Concepts

#### Forward Transport: "Where are cells going?"

A distribution of cells {{< math >}} $\mu$ {{< /math >}} can be pushed through the transport map to determine future cell states.

Starting with a distribution over some cells {{< math >}} $\mu$ {{< /math >}}, for example:
- All cells in a particular cluster or cell type
- Or a uniform distribution over some selected cells

Pushing forward this distribution through the transport map {{< math >}} $\pi_e$ {{< /math >}} gives the distribution of where these cells go in the future:

{{< math >}}
$$\mu_{\text{des}} = \mu \cdot \pi_e \tag{1}$$
{{< /math >}}

This answers the question: "Given that my cells are here now, where are they going?"

#### Reverse Transport: "Where did cells come from?"

Conversely, a distribution {{< math >}} $\mu$ {{< /math >}} can be pulled back to determine cell origins.

You take a set of cells at a later time or later state, and ask: "Where did they come from?"

You pull back the distribution via the transpose of the transport matrix:

{{< math >}}
$$\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T \tag{2}$$
{{< /math >}}

The transpose reverses the direction of flow -- it computes how much each earlier state contributed to the current cell distribution.

#### Indicator-Based Cell Selection

For a set of cells {{< math >}} $S = \{s_1, \ldots, s_n\}$ {{< /math >}}, the distribution {{< math >}} $\mu$ {{< /math >}} is defined using an indicator vector:

{{< math >}}
$$\mu_i = \mathbf{1}_{s_i \in S} \tag{3}$$
{{< /math >}}

This means:
- {{< math >}} $\mu_i = 1$ {{< /math >}} if cell {{< math >}} $i$ {{< /math >}} is in the set {{< math >}} $S$ {{< /math >}}
- {{< math >}} $\mu_i = 0$ {{< /math >}} otherwise

The input distribution puts mass 1 on each selected cell, and 0 elsewhere.

### 3. Biological Cases

**Progenitor Analysis:**
1. Select a set of progenitor cells (set {{< math >}} $S$ {{< /math >}})
2. Compute {{< math >}} $\mu_{\text{des}} = \mu \cdot \pi_e$ {{< /math >}} → get the expected distribution of their descendants

**Lineage Tracing:**
1. Select differentiated cells
2. Compute {{< math >}} $\mu_{\text{anc}} = \mu \cdot \tilde{\pi}^T$ {{< /math >}} → get their likely origins

This framework helps answer:
- "Where are these cells going?"
- "Where did these cells come from?"

### 4. Numerical Example: Forward Transport

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
$$\mu_{\text{des}} = [0.6 + 0.2, 0.2 + 0.1, 0.1 + 0.6, 0.1 + 0.1] = [0.8, 0.3, 0.7, 0.2] $$
{{< /math >}}

This gives a distribution over where the selected cells (1 & 3) are headed -- their descendant cell states.

### 5. Numerical Example: Reverse Transport
#### Step 1: Compute the Ancestor Distribution

To find the ancestor distribution, we use the transpose of {{< math >}} $\pi_e$ {{< /math >}}:

{{< math >}} 
$$(\pi_e)^T = \begin{bmatrix}
0.6 & 0.1 & 0.2 & 0.3 \\
0.2 & 0.7 & 0.1 & 0.2 \\
0.1 & 0.1 & 0.6 & 0.1 \\
0.1 & 0.1 & 0.1 & 0.4
\end{bmatrix} \tag{3}$$
{{< /math >}}

#### Step 2: Calculate the Ancestor Distribution

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

#### Step 3: Resulting Ancestor Distribution

{{< math >}} 
$$\mu_{\text{anc}} = [0.71, 0.48, 0.55, 0.26] \tag{9}$$
{{< /math >}}

**Interpretation**: Given that the descendant distribution was {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}}, pulling it back through the transpose of the transport map estimates the likely ancestor distribution of these descendants. The weights show which earlier states contributed more to this descendant cell population.

#### Normalization

We can normalize {{< math >}} $\mu_{\text{anc}}$ {{< /math >}} to turn it into a proper probability distribution:

{{< math >}} 
$$\text{Sum} = 0.71 + 0.48 + 0.55 + 0.26 = 2.00 $$
{{< /math >}}

{{< math >}} 
$$\mu_{\text{anc,normalized}} = [0.355, 0.24, 0.275, 0.13] $$
{{< /math >}}

### 6. Why Perfect Inversion is Impossible

The descendant distribution {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}} cannot lead exactly back to the original indicator vector {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} when using the transpose of the transport map for several fundamental reasons:

#### Non-invertibility of the Transport Map

The transport matrix {{< math >}} $\pi_e$ {{< /math >}} typically represents a probabilistic, many-to-many transition between cell states. It's a stochastic matrix where rows sum to 1, but it's often not invertible or the inverse is not unique. The operation of pushing forward loses information through mixing and averaging, and pulling back is not a perfect inverse.

#### Smoothing and Mixing

When we multiplied {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} by {{< math >}} $\pi_e$ {{< /math >}}, we obtained a weighted average of rows 1 and 3, resulting in the smoothed vector {{< math >}} $[0.8, 0.3, 0.7, 0.2]$ {{< /math >}}. This vector is continuous-valued and represents a distribution over many possible descendant states. The reverse step with the transpose spreads these probabilities back to ancestors but cannot undo the smoothing perfectly.

#### Loss of Sharpness

The original vector {{< math >}} $[1, 0, 1, 0]$ {{< /math >}} is a sharp indicator (only cells 1 and 3 are selected). After forward mapping, the distribution reflects uncertainty or spreading. Pulling back cannot perfectly re-sharpen because the transport matrix mixes information.

#### Biological Analogy

In biological systems, a population of cells can differentiate into multiple descendant states probabilistically. The future cell states reflect a blend of many ancestors. Hence, you cannot simply invert descendants back to a unique set of ancestors -- instead, you get a probabilistic ancestor distribution that represents the uncertainty inherent in the biological process.

---
## C) Gene and cell dependence in the calculation of latent time 
scVelo's latent time is calculated as:

{{< math >}} $$t_{i;o} = Q^p_g \left( t_{i;g} - t_{o;g} \right)$$ {{< /math >}}

Where:
- {{< math >}} $t_{i;o}$ {{< /math >}}: latent time for cell {{< math >}} $i$ {{< /math >}} relative to root cell {{< math >}} $o$ {{< /math >}}
- {{< math >}} $Q^p_g$ {{< /math >}}: {{< math >}} $p$ {{< /math >}}-quantile computed for gene {{< math >}} $g$ {{< /math >}}
- {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}: time shift between cell {{< math >}} $i$ {{< /math >}} and root cell {{< math >}} $o$ {{< /math >}} for gene {{< math >}} $g$ {{< /math >}}

Please note {{< math >}} $t_{i;g}$ {{< /math >}} and {{< math >}} $t_{o;g}$ {{< /math >}} have been normalized to a common overall time scale, which worths further discussion in the next blog.

### 1. Is p Gene-Specific?

**❌ No** -- {{< math >}} $p$ {{< /math >}} is **not gene-specific**.

#### Explanation

{{< math >}} $Q^p_g$ {{< /math >}} means: for **each gene** {{< math >}} $g$ {{< /math >}}, you compute the **{{< math >}} $p$ {{< /math >}}-quantile** over the time shifts {{< math >}} $t_{i;g} - t_{o;g}$ {{< /math >}}, across all cells {{< math >}} $i$ {{< /math >}}.

But the value of **{{< math >}} $p$ {{< /math >}}** itself is the **same for all genes** -- it is **not gene-specific**.

So even though the quantile is applied **per gene** (i.e., each gene gets its own {{< math >}} $Q^p_g$ {{< /math >}}), the **{{< math >}} $p$ {{< /math >}} used in all those calculations is shared**.

#### Mathematical Illustration

For a fixed value {{< math >}} $p = 0.5$ {{< /math >}} (median):

{{< math >}} $$Q^{0.5}_{gene1} = \text{median}(\{t_{i;gene1} - t_{o;gene1}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{gene2} = \text{median}(\{t_{i;gene2} - t_{o;gene2}\}_{i=1}^N)$$ {{< /math >}}

{{< math >}} $$\vdots$$ {{< /math >}}

{{< math >}} $$Q^{0.5}_{geneM} = \text{median}(\{t_{i;geneM} - t_{o;geneM}\}_{i=1}^N)$$ {{< /math >}}

The **same** {{< math >}} $p = 0.5$ {{< /math >}} is used for all genes, but each gene gets its own quantile value.

### 2. Is p Cell-Specific?

**❌ Also no** -- {{< math >}} $p$ {{< /math >}} is not adjusted for each cell. It is only used to **transform the gene-wise temporal shifts** into a consistent time estimate **for each cell**, **relative to the root cell** {{< math >}} $o$ {{< /math >}}.

#### Explanation

The only place where indexing over cells happens is in computing the latent time {{< math >}} $t_{i;o}$ {{< /math >}}, but the **value of {{< math >}} $p$ {{< /math >}} stays fixed** for all cells.

#### Mathematical Process

1. **Fixed {{< math >}} $p$ {{< /math >}} across all computations**: {{< math >}} $p = p_{\text{fixed}}$ {{< /math >}}

2. **Gene-wise quantile computation**:
   {{< math >}} $$Q^{p_{\text{fixed}}}_g = \text{quantile}_{p_{\text{fixed}}}(\{t_{i;g} - t_{o;g}\}_{i=1}^N)$$ {{< /math >}}

3. **Cell-specific time assignment**:
   {{< math >}} $$t_{1;o} = Q^{p_{\text{fixed}}}_g \left( t_{1;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$t_{2;o} = Q^{p_{\text{fixed}}}_g \left( t_{2;g} - t_{o;g} \right)$$ {{< /math >}}
   {{< math >}} $$\vdots$$ {{< /math >}}
   {{< math >}} $$t_{N;o} = Q^{p_{\text{fixed}}}_g \left( t_{N;g} - t_{o;g} \right)$$ {{< /math >}}

### 3. Is p Root-Specific?

**✅ Yes** -- scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness.

#### Optimization Process

For each potential root cell {{< math >}} $o$ {{< /math >}}, scVelo optimizes:

{{< math >}} $$p^*_o = \arg\min_p \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

Where the loss function might be:

{{< math >}} $$\mathcal{L}_{\text{smoothness}}(p, o) = \sum_{i,j \in \text{neighbors}} \left| t_{i;o}(p) - t_{j;o}(p) \right|^2$$ {{< /math >}}

#### Root Selection

The final root and {{< math >}} $p$ {{< /math >}} are chosen as:

{{< math >}} $$o^*, p^* = \arg\min_{o,p} \mathcal{L}_{\text{smoothness}}(p, o)$$ {{< /math >}}

### 4. Summary Table

| Aspect | Is {{< math >}} $p$ {{< /math >}} dependent? | Explanation |
|--------|-----|-------------|
| **Per gene** {{< math >}} $g$ {{< /math >}} | ❌ No | The quantile {{< math >}} $Q^p$ {{< /math >}} is computed per gene, but **{{< math >}} $p$ {{< /math >}} itself is fixed** |
| **Per cell** {{< math >}} $i$ {{< /math >}} | ❌ No | All cells use the same {{< math >}} $p$ {{< /math >}}; only their computed values differ |
| **Per root** {{< math >}} $o$ {{< /math >}} | ✅ Yes | scVelo searches for the best {{< math >}} $p$ {{< /math >}} **per root candidate** to optimize time smoothness |


### 5. Optimization Step

Once **scVelo** computes all {{< math >}} $ t_{i,o} $ {{< /math >}} for different candidate root cells {{< math >}} $ o $ {{< /math >}}, it looks for the root and quantile {{< math >}} $ p $ {{< /math >}} that maximize correlation between the **latent time series**  
{{< math >}} 
$$ \{ t_{i,o} \} $$ 
{{< /math >}}  
and its **convolution over a local neighborhood graph** (e.g., KNN graph of cells).

#### **Mathematical Formulation**
The optimization step seeks:

{{< math >}} 
$$ \arg\max_{o,p} \text{corr} \left( t_{i,o}, \sum_{j \in N(i)} w_{ij} t_{j,o} \right) $$ 
{{< /math >}}

where:
- {{< math >}} $ N(i) $ {{< /math >}} is the **neighborhood of cell** {{< math >}} $ i $ {{< /math >}}.
- {{< math >}} $ w_{ij} $ {{< /math >}} represents **graph weights**, determining connectivity between cells.


---
## D) Neighborhood-convolved Latent Time in scVelo

### Step-by-step explanation

#### 1. **Latent time**
* scVelo computes a **latent time** {{< math >}} $t_n$ {{< /math >}} for each cell {{< math >}} $n$ {{< /math >}}, indicating **how far along a differentiation trajectory** it is, starting from the inferred root cells.
* This time is estimated based on how consistent a cell's transcriptional state is with the **direction of gene expression changes (RNA velocity)**.

#### 2. **Gene-shared latent time course**
* This is the original latent time vector computed across all cells.
* It's "gene-shared" because it's not just one gene's trajectory -- it aggregates information across genes to produce a **common timeline**.

#### 3. **Neighborhood convolution**
* scVelo smooths the latent time values using **k-nearest neighbor (KNN)** information (e.g., based on transcriptomic similarity).
* The **convolved time** for a cell is the **average of the latent times of its neighboring cells**.

Think of it as:

{{< math >}} 
$$\tilde{t}_n = \frac{1}{|N(n)|} \sum_{m \in N(n)} t_m \tag{1}$$
{{< /math >}}

where {{< math >}} $N(n)$ {{< /math >}} is the neighborhood of cell {{< math >}} $n$ {{< /math >}}.

#### 4. **Regression & consistency check**
* scVelo performs **linear regression**: it tries to predict each cell's latent time using its convolved value.
* If a cell's latent time significantly **deviates from the trend of its neighbors**, it is likely noisy or unreliable.

#### 5. **Correction / Replacement**
* For these **outlier cells**, scVelo **replaces their latent time** with the convolved (smoothed) version.
* This acts like a denoising step -- ensuring latent time values are **locally smooth and robust** to noise.

### Numeric Explanation of Neighborhood Convolution

Let's say we have **5 cells**, and we already computed latent times for them as:

| Cell | Latent Time {{< math >}} $t_n$ {{< /math >}} |
|------|----------------------------------------------|
| A    | 0.10                                        |
| B    | 0.15                                        |
| C    | 0.85 ❗                                     |
| D    | 0.20                                        |
| E    | 0.25                                        |

Clearly, cell **C** looks suspicious -- it's way ahead in time compared to its neighbors.

#### Step 1: Define neighbors

Let's say the **neighborhoods (e.g., via KNN)** are:
* A → neighbors: [B, D]
* B → neighbors: [A, D]
* C → neighbors: [B, D, E]
* D → neighbors: [A, B, E]
* E → neighbors: [C, D]

#### Step 2: Compute convolved (smoothed) latent times

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

#### Step 3: Compare original vs convolved times

| Cell | {{< math >}} $t_n$ {{< /math >}} | {{< math >}} $\tilde{t}_n$ {{< /math >}} | Difference |
|------|-----------------------------------|-------------------------------------------|------------|
| A    | 0.10                             | 0.175                                     | 0.075      |
| B    | 0.15                             | 0.15                                      | 0.000      |
| C    | 0.85 ❗                          | 0.20                                      | 0.65 ❗    |
| D    | 0.20                             | 0.167                                     | 0.033      |
| E    | 0.25                             | 0.525 ❗                                  | 0.275 ❗   |

#### Step 4: Detect inconsistent cells

* scVelo runs a **regression of** {{< math >}} $t_n \sim \tilde{t}_n$ {{< /math >}} and checks for large residuals.
* Cells **C** and **E** have large mismatches between original and smoothed values → likely noisy or incorrect.

#### Step 5: Replace inconsistent values

Let's say we define a cutoff: if difference > 0.3, replace with {{< math >}} $\tilde{t}_n$ {{< /math >}}.

| Cell | New {{< math >}} $t_n$ {{< /math >}} |
|------|--------------------------------------|
| A    | 0.10 (unchanged)                    |
| B    | 0.15 (unchanged)                    |
| C    | **0.20** (replaced!)               |
| D    | 0.20 (unchanged)                    |
| E    | **0.525** (replaced!)              |

#### Final corrected latent time:

{{< math >}} 
$$\mathbf{t} = [0.10, 0.15, \mathbf{0.20}, 0.20, \mathbf{0.525}] \tag{3}$$
{{< /math >}}

Cell C is now brought **back in line** with its local trajectory.