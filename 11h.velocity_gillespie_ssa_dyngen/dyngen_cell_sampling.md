---
title: "Dyngen Cell Sampling"
date: 2025-07-10
draft: True
---

# A) Cell Sampling Methods

## Backbone-Based Sampling from Multiple Simulations

### 📌 Context
When dealing with simulated developmental trajectories, we need to represent the process compactly using snapshots. The "backbone" serves as a simplified representation, a skeleton path through cell states consisting of discrete states connected by transition edges.

### 🧠 Method Breakdown

#### 1. Backbone Structure
The backbone consists of several states linked together by transition edges with length {{< math >}} $L_i$ {{< /math >}}.

Each edge represents a transition from one cell state to another, where {{< math >}} $L_i$ {{< /math >}} can represent:
- A time duration
- A pseudotime interval  
- Any metric capturing transition "length"

#### 2. Simulation Mapping
Multiple stochastic simulations are run, and their outputs (simulated cell states) are projected onto this backbone. Each simulation step is assigned to a point along one of the edges in the backbone.

#### 3. Proportional Sampling
From each transition, {{< math >}} $$N_i = N \cdot \frac{L_i}{\sum L_i}$$ {{< /math >}} cells are sampled uniformly.

This ensures:
- **Fair allocation**: Samples are distributed proportionally to edge length
- **Longer transitions get more samples**: Reflecting their relative importance
- **Proportional representation**: If edge {{< math >}} $L_i$ {{< /math >}} is 20% of total trajectory length, it gets 20% of samples

#### 4. Integer Adjustment
After calculating {{< math >}} $N_i$ {{< /math >}}, we may get decimals (e.g., 8.4, 3.6). We:
1. Round each to the nearest integer
2. Adjust totals to ensure {{< math >}} $$\sum N_i = N$$ {{< /math >}}

### ✅ Summary
This process ensures:
- All trajectory parts are fairly represented
- Sampling is proportional to transition length
- Final snapshot contains exactly {{< math >}} $N$ {{< /math >}} cells
- Uniform sampling prevents bias toward fast or slow regions

---

## Time Series Snapshots from Molecular Dynamics

### 📌 Key Concepts
For simulations running from time 0 to {{< math >}} $T$ {{< /math >}}, we extract representative snapshots evenly spaced across time with clear temporal separation.

### 🧠 Detailed Method

#### 1. Simulation Duration
The simulation ends at time point {{< math >}} $T$ {{< /math >}} (often between 10 and 20 in dyngen).

#### 2. Interval Division
The interval {{< math >}} $[0, T]$ {{< /math >}} is divided into {{< math >}} $k$ {{< /math >}} equal intervals of width {{< math >}} $w$ {{< /math >}} separated by {{< math >}} $k-1$ {{< /math >}} gaps of width {{< math >}} $g$ {{< /math >}}.

The total time relationship is:
{{< math >}} $$T = k \cdot w + (k-1) \cdot g$$ {{< /math >}}

Solving for {{< math >}} $w$ {{< /math >}}:
{{< math >}} $$w = \frac{T - (k-1)g}{k}$$ {{< /math >}}

These non-overlapping windows ensure snapshots from distinct time periods, providing better coverage of dynamics.

#### 3. Cell Sampling
{{< math >}} $N_i = N/k$ {{< /math >}} cells are sampled uniformly from each interval, where:
- {{< math >}} $N$ {{< /math >}} = total cells to represent the system
- Each time window contributes equally
- Cells are sampled uniformly within their respective intervals

#### 4. Integer Correction
After calculating {{< math >}} $N_i$ {{< /math >}}, round and adjust so:
{{< math >}} $$\sum N_i = N$$ {{< /math >}}

#### 5. Default Parameters
By default: {{< math >}} $k = 8$ {{< /math >}} and {{< math >}} $g = 0.75$ {{< /math >}}

This setup:
- Breaks {{< math >}} $[0, T]$ {{< /math >}} into 8 intervals
- Has 7 gaps between them, each 0.75 units wide
- Total gap time: {{< math >}} $7 \times 0.75 = 5.25$ {{< /math >}}

The usable sampling span becomes:
{{< math >}} $$T = 8w + 5.25 \Rightarrow w = \frac{T - 5.25}{8}$$ {{< /math >}}

**Example**: If {{< math >}} $T = 15$ {{< /math >}}, then {{< math >}} $w = \frac{15 - 5.25}{8} = 1.21875$ {{< /math >}}

#### 6. Parameter Scaling
For usual dyngen simulations: {{< math >}} $10 \leq T \leq 20$ {{< /math >}}

For larger values of {{< math >}} $T$ {{< /math >}}, both {{< math >}} $k$ {{< /math >}} and {{< math >}} $g$ {{< /math >}} should be increased accordingly to:
- Maintain temporal resolution
- Avoid undersampling early and late dynamics

### ✅ Summary
This method ensures:
- Time series snapshots are evenly spaced and non-overlapping
- Each system stage is well represented  
- Parameters can be tuned for longer simulations or finer temporal resolution


# B) Molecular Sampling Methods

This document explains how synthetic (in silico) single-cell RNA-seq data is simulated to resemble real experimental data. Here's a breakdown of what each step does:

## 🔢 Step-by-Step Explanation

### 1. Library Size Sampling per Cell (Cell-level total counts)

"For each in silico cell {{< math >}} $i$ {{< /math >}}, draw its library size {{< math >}} $ls_i$ {{< /math >}} from the distribution of transcript counts per cell in the real dataset."

The library size is the total number of RNA molecules (UMIs or reads) detected in a single cell.

To keep synthetic data realistic, you sample {{< math >}} $ls_i$ {{< /math >}} from the real distribution of observed transcript counts in a real dataset.

**Example:** If real cells range from 1000 to 5000 transcripts, then each synthetic cell's total count is drawn from this distribution.

### 2. Capture Rate Sampling per Molecule Type (Gene-level efficiency)

"The capture rate {{< math >}} $cr_j$ {{< /math >}} of each in silico molecule type {{< math >}} $j$ {{< /math >}} is drawn from {{< math >}} $N(1, 0.05)$ {{< /math >}}."

Capture rate {{< math >}} $cr_j$ {{< /math >}} represents the probability that a molecule of type (gene) {{< math >}} $j$ {{< /math >}} is captured during sequencing.

These are sampled from a normal distribution centered at 1, with 5% standard deviation, reflecting small, random variability in capture efficiency across genes.

Some genes are slightly more or less likely to be captured than others.

### 3. Simulated Read Count Generation (Per cell × gene matrix)

"Finally, for each cell {{< math >}} $i$ {{< /math >}}, draw {{< math >}} $ls_i$ {{< /math >}} molecules from the multinomial distribution with probabilities {{< math >}} $cr_j \times ab_{i,j}$ {{< /math >}}"

Now you simulate the actual observed counts per gene.

For each cell:

- You have the true underlying abundance {{< math >}} $ab_{i,j}$ {{< /math >}} for every molecule {{< math >}} $j$ {{< /math >}} (e.g., precomputed from the model).

- You scale these abundances by the capture rates {{< math >}} $cr_j$ {{< /math >}}, which modifies the effective expression level per gene.

- Normalize the resulting vector to make it a probability distribution.

- Sample {{< math >}} $ls_i$ {{< /math >}} molecules from this multinomial distribution, assigning them to genes.

## 🔁 In Simpler Terms

For cell {{< math >}} $i$ {{< /math >}}, simulate gene expression counts by:

1. Adjusting the "true" gene abundance by technical noise (capture rates),
2. Drawing a total of {{< math >}} $ls_i$ {{< /math >}} reads (like the real cells),
3. Assigning reads to genes probabilistically.

## ✅ Why This Works Well

- It mimics realistic technical noise (capture rates).
- It preserves realistic library size variation across cells.
- It produces synthetic count matrices that closely resemble real scRNA-seq data,useful for benchmarking methods like clustering, trajectory inference, or denoising.

In this context, the multinomial distribution models the process of randomly assigning a fixed number of transcript reads (or molecules) to different genes based on their relative probabilities (i.e. expression levels adjusted by capture rates).

## 🔍 What is a Multinomial Distribution?

The multinomial distribution is a generalization of the binomial distribution.

- **Binomial:** How many times does a single event (e.g. heads or tails) occur over {{< math >}} $n$ {{< /math >}} trials.
- **Multinomial:** How many times each of multiple events (e.g. heads, tails, edge) occurs over {{< math >}} $n$ {{< /math >}} trials.

Formally:

If you do {{< math >}} $n$ {{< /math >}} independent trials, and each trial results in one of {{< math >}} $k$ {{< /math >}} categories with probabilities {{< math >}} $p_1, p_2, \ldots, p_k$ {{< /math >}} (where {{< math >}} $\sum p_j = 1$ {{< /math >}}), then:

{{< math >}} 
$$\text{Multinomial}(n; p_1, p_2, \ldots, p_k)$$
{{< /math >}}

gives the number of times each category (gene) occurs.

## 📌 In Your Case

Each cell {{< math >}} $i$ {{< /math >}} has:

- A total transcript count {{< math >}} $ls_i$ {{< /math >}} (library size).

- A set of probabilities {{< math >}} $p_j = \frac{cr_j \cdot ab_{i,j}}{\sum_{j'} cr_{j'} \cdot ab_{i,j'}}$ {{< /math >}} for each molecule (gene) {{< math >}} $j$ {{< /math >}}.

These reflect how likely a molecule from gene {{< math >}} $j$ {{< /math >}} is captured, based on abundance and capture efficiency.

Then you sample:

{{< math >}} 
$$[x_{i,1}, x_{i,2}, \ldots, x_{i,k}] \sim \text{Multinomial}(ls_i; p_1, p_2, \ldots, p_k)$$
{{< /math >}}

Where:

- {{< math >}} $x_{i,j}$ {{< /math >}} is the number of times gene {{< math >}} $j$ {{< /math >}} is detected in cell {{< math >}} $i$ {{< /math >}}.
- The result is one row of the synthetic count matrix, i.e. the observed counts for all genes in cell {{< math >}} $i$ {{< /math >}}.

## 🔧 Example

Let's say cell {{< math >}} $i$ {{< /math >}} has a library size of 100:

Adjusted (normalized) probabilities:
- Gene A: 0.2
- Gene B: 0.5  
- Gene C: 0.3

Then:

{{< math >}} 
$$[x_A, x_B, x_C] \sim \text{Multinomial}(100; 0.2, 0.5, 0.3)$$
{{< /math >}}

Could result in something like: [19, 52, 29], one random realization of the process.