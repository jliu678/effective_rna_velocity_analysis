---
title: "Unraveled Brust Cme Simulation Fourier Computation"
date: 2025-07-10
draft: True
---

# A) Step-by-Step Summary of the Simulation Pipeline for Modeling Transcriptional Bursting Using the Chemical Master Equation (CME)

## 🔁 Overall Goal

Simulate steady-state mRNA/protein distributions for multiple genes across single cells, under a bursting gene expression model, using CME-based probability distributions.

## 1. Function: `simulate_burst_model(...)`

### 📌 Inputs:
- {{< math >}} $\text{nCells}$ {{< /math >}}: Total number of simulated cells
- {{< math >}} $\text{nGenes}$ {{< /math >}}: Number of genes
- {{< math >}} $T$ {{< /math >}}: (Unused) placeholder
- {{< math >}} $\text{n\_cell\_types}$ {{< /math >}}: Number of distinct cell types
- {{< math >}} $\text{seed}$ {{< /math >}}: For reproducibility

### 🔧 Step-by-Step Simulation Logic:

#### a. Initialize Variables:
- {{< math >}} $X$ {{< /math >}}: shape {{< math >}} $(2, \text{nCells}, \text{nGenes})$ {{< /math >}} – stores molecule counts per gene per cell
- {{< math >}} $g_{\text{true}}[j,i]$ {{< /math >}}: degradation rate for gene {{< math >}} $j$ {{< /math >}} in cell type {{< math >}} $i$ {{< /math >}}
- {{< math >}} $b_{\text{true}}[j]$ {{< /math >}}: burst size for gene {{< math >}} $j$ {{< /math >}}
- {{< math >}} $K[j,i]$ {{< /math >}}: burst frequency (initiation rate) for gene {{< math >}} $j$ {{< /math >}} in cell type {{< math >}} $i$ {{< /math >}}

#### b. Loop Over Genes:
For each gene {{< math >}} $j$ {{< /math >}}:

1. Sample a burst size {{< math >}} $b$ {{< /math >}} (log-normal)
2. Sample a mean degradation rate {{< math >}} $\bar{\gamma}$ {{< /math >}}

🔁 Then loop over each cell type {{< math >}} $i$ {{< /math >}}:

1. Sample:
   - Burst frequency {{< math >}} $k_{\text{ini}}$ {{< /math >}}
   - Degradation rate {{< math >}} $\gamma$ {{< /math >}} around {{< math >}} $\bar{\gamma}$ {{< /math >}}

2. Compute mean vector:

{{< math >}} 
$$\mu_1 = \frac{k_{\text{ini}} b}{\beta} \tag{1}$$
{{< /math >}}

{{< math >}} 
$$\mu_2 = \frac{k_{\text{ini}} b}{\gamma} \tag{2}$$
{{< /math >}}

3. Compute variances and set a support limit {{< math >}} $l_m$ {{< /math >}} based on a safety margin

4. Call `cme_integrator(...)` to compute steady-state distribution over molecule counts

5. Use `rv_discrete` to sample from CME-distribution

6. Store samples in {{< math >}} $X[:, \text{IND}, j]$ {{< /math >}} for all cells in this type

## 2. Function: `cme_integrator(p, lm, ...)`

### 🧮 Goal:
Compute steady-state joint distribution of two molecule species (e.g., nascent and mature RNA), using inverse Fourier transform of generating function.

### 🔧 Key Steps:

#### a. Define Generating Function Grid:
- Set up {{< math >}} $u$ {{< /math >}} for each dimension (Fourier basis)
- Compute evaluation mesh {{< math >}} $g = [g_0, g_1]$ {{< /math >}}

#### b. Adjust Generating Function:
- Compute weights from burst model
- Transform the mesh to represent the generating function of the CME

#### c. Numerical Integration:
Use fixed Gaussian quadrature (or optionally `quad_vec`) to integrate over time.

Evaluate:

{{< math >}} 
$$\text{gf} = \int_0^T \frac{e^{-\beta t} g_0 + e^{-\gamma t} g_1}{1 - (e^{-\beta t} g_0 + e^{-\gamma t} g_1)} dt \tag{3}$$
{{< /math >}}

#### d. Back-Transform to Probability:
- Apply exponential and inverse FFT: `irfft2(...)`
- Normalize and clip small values
- Returns steady-state probability mass function over molecule count pairs

## 3. Function: `INTFUN(x, g, β, γ)`

Core function for integration in the generating function:

{{< math >}} 
$$\frac{e^{-\beta x} g_0 + e^{-\gamma x} g_1}{1 - (e^{-\beta x} g_0 + e^{-\gamma x} g_1)} \tag{4}$$
{{< /math >}}

This function captures the input to the integral representing the generating function's logarithm over time.

## ✅ Final Outputs:
- {{< math >}} $X$ {{< /math >}}: Simulated count data (2 species × cells × genes)
- `cell_types`: Cell type labels per cell
- {{< math >}} $K$ {{< /math >}}: Burst frequencies
- {{< math >}} $g_{\text{true}}$ {{< /math >}}: Degradation rates
- {{< math >}} $b_{\text{true}}$ {{< /math >}}: Burst sizes

## 📌 Summary Diagram (Conceptual):

```
Gene j + Cell type i
    ↓ (sample parameters: b, γ, k_ini)
CME → Generating Function (via INTFUN)
    ↓
Inverse FFT → P(n₁, n₂)
    ↓
Sampling → X[:, cells, j]
```
# B) Biological RNA Model: Transcription and Splicing Dynamics

## 🧬 Biological Interpretation

The model assumes that:

- Transcription initiation occurs at rate {{< math >}} $k_{\text{ini}}$ {{< /math >}}, leading to bursts of RNA production.

- Each burst produces a random number of nascent transcripts (pre-mRNA), with mean burst size {{< math >}} $b$ {{< /math >}}.

- Pre-mRNA is processed/degraded at rate {{< math >}} $\gamma$ {{< /math >}} (e.g., splicing, decay).

- Processed mRNA is degraded at rate {{< math >}} $\beta$ {{< /math >}}.

This describes two RNA species per gene:

- **Unspliced** (pre-mRNA)
- **Spliced** (mature mRNA)

## ⚙️ Parameter Definitions

In your model:

- {{< math >}} $k_{\text{ini}}$ {{< /math >}}: transcription initiation rate (bursts per unit time)
- {{< math >}} $b$ {{< /math >}}: mean burst size (molecules per burst)
- {{< math >}} $\beta$ {{< /math >}}: decay rate of spliced mRNA
- {{< math >}} $\gamma$ {{< /math >}}: decay/splicing rate of unspliced (pre-mRNA)

So, the parameter vector is:

{{< math >}} 
$$\mathbf{p} = [k_{\text{ini}}, b, \beta, \gamma] \tag{1}$$
{{< /math >}}

## 📐 The Model Equations (Stochastic View)

Let {{< math >}} $u$ {{< /math >}} be the unspliced (pre-mRNA), and {{< math >}} $s$ {{< /math >}} be the spliced mRNA.

Then, in a stochastic formulation:

**Transcription:**
{{< math >}} 
$$\emptyset \xrightarrow{k_{\text{ini}}} b \cdot u \tag{2}$$
{{< /math >}}
(i.e., a burst of size ~Geom or ~ExpPoisson of {{< math >}} $b$ {{< /math >}})

**Splicing:**
{{< math >}} 
$$u \xrightarrow{\gamma} s \tag{3}$$
{{< /math >}}

**Degradation:**
{{< math >}} 
$$u \xrightarrow{\gamma} \emptyset, \quad s \xrightarrow{\beta} \emptyset \tag{4}$$
{{< /math >}}

The joint distribution {{< math >}} $P(u,s)$ {{< /math >}} is solved via the generating function approach.

## 🎯 Key Modeling Assumptions

1. **No feedback or switching**: All parameters are static per gene and per cell type.

2. **Independent genes**: No gene-gene regulation.

3. **Fast bursts**: Bursts are instantaneous; durations are neglected.

4. **Quasi–steady-state**: You're solving for the stationary distribution.

## 🔄 Generating Function Method

Instead of simulating trajectories (e.g., Gillespie SSA), your model uses Fourier transform of the generating function to directly compute the joint distribution {{< math >}} $P(u,s)$ {{< /math >}}.

The generating function approach provides:

{{< math >}} 
$$G(z_u, z_s) = \sum_{u=0}^{\infty} \sum_{s=0}^{\infty} P(u,s) z_u^u z_s^s \tag{5}$$
{{< /math >}}

This is efficient and allows accurate sampling from the true steady-state distribution under Chemical Master Equation (CME) dynamics.

## Model Summary

The complete stochastic system can be represented as:

{{< math >}} 
$$\frac{\partial P(u,s,t)}{\partial t} = k_{\text{ini}} b P(u-b,s,t) - k_{\text{ini}} P(u,s,t) + \gamma (u+1) P(u+1,s-1,t) - \gamma u P(u,s,t) + \beta (s+1) P(u,s+1,t) - \beta s P(u,s,t) \tag{6}$$
{{< /math >}}

where the steady-state solution {{< math >}} $P(u,s) = \lim_{t \to \infty} P(u,s,t)$ {{< /math >}} is obtained through the generating function method.

# C) Solving Generating Function via Method of Characteristics

## Step 1: Model Setup and Generating Function Definition

We consider a bursty transcription model with two molecular species:

- {{< math >}} $U$ {{< /math >}}: unspliced mRNA count
- {{< math >}} $S$ {{< /math >}}: spliced mRNA count

with parameters:

- {{< math >}} $k_{\text{ini}}$ {{< /math >}}: burst initiation rate (Poisson rate of bursts)
- {{< math >}} $b$ {{< /math >}}: burst size (number of transcripts per burst, modeled as geometric or fixed burst)
- {{< math >}} $\gamma$ {{< /math >}}: splicing/degradation rate of unspliced mRNA
- {{< math >}} $\beta$ {{< /math >}}: degradation rate of spliced mRNA

The joint generating function is:

{{< math >}} 
$$G(z_1, z_2, t) = \sum_{u=0}^{\infty} \sum_{s=0}^{\infty} P(u,s,t) z_1^u z_2^s \tag{1}$$
{{< /math >}}

where {{< math >}} $P(u,s,t)$ {{< /math >}} is the probability of having {{< math >}} $u$ {{< /math >}} unspliced and {{< math >}} $s$ {{< /math >}} spliced mRNA molecules at time {{< math >}} $t$ {{< /math >}}.

## Step 2: Write the CME in Generating Function Form

The CME describes the time evolution of {{< math >}} $P(u,s,t)$ {{< /math >}}. Using standard approaches, the PDE for the generating function is:

{{< math >}} 
$$\frac{\partial G}{\partial t} = k_{\text{ini}}(F(z_1, z_2) - 1)G + \gamma(z_2 - z_1)\frac{\partial G}{\partial z_1} + \beta(1 - z_2)\frac{\partial G}{\partial z_2} \tag{2}$$
{{< /math >}}

where {{< math >}} $F(z_1, z_2)$ {{< /math >}} is the generating function for the burst size distribution. For a geometric burst:

{{< math >}} 
$$F(z_1, z_2) = \frac{1}{1 - b(z_1 - 1)} \tag{3}$$
{{< /math >}}

In your model, bursts produce unspliced RNA only, so:

{{< math >}} 
$$F(z_1, z_2) = \frac{1}{1 - b(z_1 - 1)} \tag{4}$$
{{< /math >}}

## Step 3: Method of Characteristics Solution

### Given PDE:

{{< math >}} 
$$\frac{\partial G}{\partial t} = k_{\text{ini}}\left(\frac{1}{1 - b(z_1 - 1)} - 1\right)G + \gamma(z_2 - z_1)\frac{\partial G}{\partial z_1} + \beta(1 - z_2)\frac{\partial G}{\partial z_2} \tag{5}$$
{{< /math >}}

### Step 3.1: Identify Characteristic System

We write the PDE in the form:

{{< math >}} 
$$\frac{\partial G}{\partial t} + a_1(z_1, z_2)\frac{\partial G}{\partial z_1} + a_2(z_1, z_2)\frac{\partial G}{\partial z_2} = \phi(z_1, z_2)G \tag{6}$$
{{< /math >}}

This gives:

{{< math >}} 
$$a_1 = -\gamma(z_2 - z_1) \tag{7}$$
{{< /math >}}

{{< math >}} 
$$a_2 = -\beta(1 - z_2) \tag{8}$$
{{< /math >}}

{{< math >}} 
$$\phi(z_1, z_2) = k_{\text{ini}}\left(\frac{1}{1 - b(z_1 - 1)} - 1\right) \tag{9}$$
{{< /math >}}

### Step 3.2: Characteristic Equations

The method of characteristics reduces this PDE to a system of ODEs along curves where the PDE becomes an ODE. The system is:

{{< math >}} 
$$\frac{dz_1}{dt} = \gamma(z_2 - z_1) \tag{10}$$
{{< /math >}}

{{< math >}} 
$$\frac{dz_2}{dt} = \beta(1 - z_2) \tag{11}$$
{{< /math >}}

{{< math >}} 
$$\frac{dG}{dt} = k_{\text{ini}}\left(\frac{1}{1 - b(z_1 - 1)} - 1\right)G \tag{12}$$
{{< /math >}}

We now solve this system step-by-step.

### Step 3.3: Solve for {{< math >}} $z_2(t)$ {{< /math >}}

{{< math >}} 
$$\frac{dz_2}{dt} = \beta(1 - z_2) \Rightarrow \frac{dz_2}{1 - z_2} = \beta dt \tag{13}$$
{{< /math >}}

Integrate both sides:

{{< math >}} 
$$-\ln|1 - z_2| = \beta t + C_2 \Rightarrow 1 - z_2 = Ce^{-\beta t} \Rightarrow z_2(t) = 1 - Ce^{-\beta t} \tag{14}$$
{{< /math >}}

Let {{< math >}} $z_2(0) = z_{2,0}$ {{< /math >}}, then:

{{< math >}} 
$$C = 1 - z_{2,0} \Rightarrow z_2(t) = 1 - (1 - z_{2,0})e^{-\beta t} \tag{15}$$
{{< /math >}}

### Step 3.4: Solve for {{< math >}} $z_1(t)$ {{< /math >}}

{{< math >}} 
$$\frac{dz_1}{dt} = \gamma(z_2(t) - z_1) \Rightarrow \frac{dz_1}{dt} + \gamma z_1 = \gamma z_2(t) \tag{16}$$
{{< /math >}}

This is a linear ODE. Use the integrating factor:

{{< math >}} 
$$\mu(t) = e^{\gamma t} \Rightarrow \frac{d}{dt}[z_1 e^{\gamma t}] = \gamma z_2(t) e^{\gamma t} \tag{17}$$
{{< /math >}}

Substitute {{< math >}} $z_2(t)$ {{< /math >}}:

{{< math >}} 
$$\frac{d}{dt}[z_1 e^{\gamma t}] = \gamma(1 - (1 - z_{2,0})e^{-\beta t})e^{\gamma t} \tag{18}$$
{{< /math >}}

Break it into two terms:

{{< math >}} 
$$\gamma e^{\gamma t} - \gamma(1 - z_{2,0})e^{(\gamma - \beta)t} \tag{19}$$
{{< /math >}}

Integrate both:

First: {{< math >}} $\int \gamma e^{\gamma t} dt = e^{\gamma t}$ {{< /math >}}

Second: {{< math >}} $\int \gamma(1 - z_{2,0})e^{(\gamma - \beta)t} dt = \frac{\gamma(1 - z_{2,0})}{\gamma - \beta}e^{(\gamma - \beta)t}$ {{< /math >}}

So:

{{< math >}} 
$$z_1(t) e^{\gamma t} = e^{\gamma t} - \frac{\gamma(1 - z_{2,0})}{\gamma - \beta}e^{(\gamma - \beta)t} + C_1 \tag{20}$$
{{< /math >}}

Solve for {{< math >}} $z_1(t)$ {{< /math >}}:

{{< math >}} 
$$z_1(t) = 1 - \frac{\gamma(1 - z_{2,0})}{\gamma - \beta}e^{-\beta t} + C_1 e^{-\gamma t} \tag{21}$$
{{< /math >}}

Find {{< math >}} $C_1$ {{< /math >}} from initial value {{< math >}} $z_1(0) = z_{1,0}$ {{< /math >}}:

{{< math >}} 
$$z_{1,0} = 1 - \frac{\gamma(1 - z_{2,0})}{\gamma - \beta} + C_1 \Rightarrow C_1 = z_{1,0} - 1 + \frac{\gamma(1 - z_{2,0})}{\gamma - \beta} \tag{22}$$
{{< /math >}}

Now we have both {{< math >}} $z_1(t)$ {{< /math >}} and {{< math >}} $z_2(t)$ {{< /math >}}.

### Step 3.5: Solve for {{< math >}} $G(t)$ {{< /math >}}

{{< math >}} 
$$\frac{dG}{dt} = k_{\text{ini}}\left(\frac{1}{1 - b(z_1(t) - 1)} - 1\right)G \tag{23}$$
{{< /math >}}

Let's simplify the coefficient:

{{< math >}} 
$$\phi(t) = k_{\text{ini}}\left(\frac{1}{1 - b(z_1(t) - 1)} - 1\right) \tag{24}$$
{{< /math >}}

Then:

{{< math >}} 
$$\frac{dG}{dt} = \phi(t)G \Rightarrow \ln G = \int \phi(t) dt + \ln G_0 \tag{25}$$
{{< /math >}}

{{< math >}} 
$$\Rightarrow G(t) = G_0 \cdot \exp\left(\int_0^t \phi(s) ds\right) \tag{26}$$
{{< /math >}}

This integral is evaluated numerically in implementation via numerical integration methods.

## Summary

The complete solution consists of:

1. **Characteristic curves**: {{< math >}} $z_1(t)$ {{< /math >}} and {{< math >}} $z_2(t)$ {{< /math >}} given by equations (15) and (21)
2. **Generating function evolution**: {{< math >}} $G(t)$ {{< /math >}} given by equation (26)
3. **Numerical integration**: The exponential integral in equation (26) requires numerical evaluation

# D) Discrete Fourier Transform (DFT) then Inverse Fast Fourier Transform (IFFT) for Probability Distribution Extraction

For this class of models, the generating function at steady state can be expressed as:

{{< math >}}
$$G(z_1, z_2) = \exp\left(k_{\text{ini}} \int_0^{\infty} \frac{b(e^{-\beta t}(z_1 - 1) + e^{-\gamma t}(z_2 - 1))}{1 - b(e^{-\beta t}(z_1 - 1) + e^{-\gamma t}(z_2 - 1))} dt\right) \tag{1}$$
{{< /math >}}

This integral expression comes from solving the CME with generating functions and using the fact that bursts are Poisson events with exponentially distributed waiting times.

The numerator and denominator inside the integral come from the transition probabilities of unspliced and spliced RNAs decaying/exiting.

The exponential decay terms {{< math >}}$e^{-\beta t}${{< /math >}}, {{< math >}}$e^{-\gamma t}${{< /math >}} reflect lifetimes of the RNA species.

The integral sums contributions over all times {{< math >}}$t${{< /math >}}, weighted by the dynamics.

## Correspondence with Your Code

In your code:

{{< math >}}
$$g[0] = z_1 - 1 \tag{2}$$
{{< /math >}}

{{< math >}}
$$g[1] = z_2 - 1 \tag{3}$$
{{< /math >}}

where {{< math >}}$z_1, z_2${{< /math >}} are points on the complex unit circle:

{{< math >}}
$$z_i = e^{-2\pi i l/L_i} \tag{4}$$
{{< /math >}}

for discrete indices {{< math >}}$l${{< /math >}} and {{< math >}}$L_i${{< /math >}} total points.

The function

```python
INTFUN(t, g, bet, gam) = (e^{-β t} g[0] + e^{-γ t} g[1]) / (1 - (e^{-β t} g[0] + e^{-γ t} g[1]))
```

matches the integrand in the integral expression of {{< math >}}$G${{< /math >}}.

## Step 0: You Start with a Continuous Generating Function

You are given the generating function of a stochastic process (such as transcription with bursting):

{{< math >}}
$$G(z_1, z_2) = \exp\left(k_{\text{ini}} \int_0^{\infty} \frac{b(e^{-\beta t}(z_1 - 1) + e^{-\gamma t}(z_2 - 1))}{1 - b(e^{-\beta t}(z_1 - 1) + e^{-\gamma t}(z_2 - 1))} dt\right) \tag{5}$$
{{< /math >}}

Where:

- {{< math >}}$z_1, z_2 \in \mathbb{C}${{< /math >}} (usually on the unit circle)
- {{< math >}}$k_{\text{ini}}${{< /math >}} is the transcription initiation rate
- {{< math >}}$\beta, \gamma${{< /math >}} are degradation/splicing rates
- {{< math >}}$b(x)${{< /math >}} is a function (e.g., {{< math >}}$b(x) = 1 + x${{< /math >}} for geometric bursting)

This function encodes the joint distribution {{< math >}}$P(u,s)${{< /math >}} over unspliced and spliced molecule counts in its power series:

{{< math >}}
$$G(z_1, z_2) = \sum_{u=0}^{\infty} \sum_{s=0}^{\infty} P(u,s) z_1^u z_2^s \tag{6}$$
{{< /math >}}

## Step 1: Truncate and Discretize the State Space

Since you can't compute infinitely many coefficients, you truncate the domain to a finite grid:

- Let {{< math >}}$u = 0, \ldots, N_1 - 1${{< /math >}}
- Let {{< math >}}$s = 0, \ldots, N_2 - 1${{< /math >}}

This gives a discrete matrix {{< math >}}$P(u,s)${{< /math >}} of size {{< math >}}$N_1 \times N_2${{< /math >}}, which we aim to compute.

## Step 2: Sample G(z₁,z₂) on the Complex Unit Torus

The discrete Fourier transform (DFT) interprets function evaluations on the unit circle in the complex plane as encoding frequency components.

Define a grid:

{{< math >}}
$$z_1^{(k)} = e^{2\pi i k/N_1}, \quad z_2^{(m)} = e^{2\pi i m/N_2} \tag{7}$$
{{< /math >}}

for {{< math >}}$k = 0, \ldots, N_1 - 1${{< /math >}} and {{< math >}}$m = 0, \ldots, N_2 - 1${{< /math >}}. These points lie evenly spaced on the unit circle.

You then compute:

{{< math >}}
$$G_{k,m} = G(z_1^{(k)}, z_2^{(m)}) = G(e^{2\pi i k/N_1}, e^{2\pi i m/N_2}) \tag{8}$$
{{< /math >}}

This yields a 2D matrix of size {{< math >}}$N_1 \times N_2${{< /math >}}, which is exactly the 2D Discrete Fourier Transform (DFT) of {{< math >}}$P(u,s)${{< /math >}}.

### How?

Because the generating function has the series:

{{< math >}}
$$G(z_1, z_2) = \sum_{u=0}^{N_1-1} \sum_{s=0}^{N_2-1} P(u,s) z_1^u z_2^s \tag{9}$$
{{< /math >}}

Plugging in {{< math >}}$z_1 = e^{2\pi i k/N_1}, z_2 = e^{2\pi i m/N_2}${{< /math >}}:

{{< math >}}
$$G_{k,m} = \sum_{u=0}^{N_1-1} \sum_{s=0}^{N_2-1} P(u,s) e^{2\pi i ku/N_1} e^{2\pi i ms/N_2} \tag{10}$$
{{< /math >}}

This is exactly the 2D DFT.

## Step 3: Use 2D Inverse FFT to Recover P(u,s)

The inverse of the 2D DFT is:

{{< math >}}
$$P(u,s) = \frac{1}{N_1 N_2} \sum_{k=0}^{N_1-1} \sum_{m=0}^{N_2-1} G_{k,m} e^{-2\pi i ku/N_1} e^{-2\pi i ms/N_2} \tag{11}$$
{{< /math >}}

This is computed efficiently using:

```python
P = np.fft.ifft2(G_vals).real
```

### Post-processing:

1. Take `.real` since numerical FFT introduces tiny imaginary parts
2. Apply `np.maximum(P, 0)` to correct negative numerical noise
3. Normalize: `P /= np.sum(P)` to ensure it sums to 1

## Why This Works — Intuition

- You treat {{< math >}}$G(z_1, z_2)${{< /math >}} as the Fourier transform of {{< math >}}$P(u,s)${{< /math >}}
- By sampling it on a grid of points {{< math >}}$z_i = e^{2\pi i k/N}${{< /math >}}, you extract the Fourier coefficients
- The inverse transform reconstructs {{< math >}}$P(u,s)${{< /math >}}, the distribution of interest
- This is numerical coefficient extraction from a bivariate generating function

## Why Not Use Analytic Coefficient Extraction?

For most generating functions of complex stochastic systems (like the one you have), it's impossible to find a closed-form expression for the coefficients {{< math >}}$P(u,s)${{< /math >}}. So we approximate the inverse using the IFFT, which is fast and accurate.


# E) `cme_integrator()`: Fourier Basis Computation

## Code Section Analysis

In the `cme_integrator` function, the section:

```python
# initialize the generating function evaluation points
mx[-1] = mx[-1] // 2 + 1
for i in range(len(mx)):
    l = np.arange(mx[i])
    u_ = np.exp(-2j * np.pi * l / lm[i]) - 1
    u.append(u_)
```

## Purpose

It sets up Fourier evaluation points {{< math >}} $u$ {{< /math >}} for each dimension of the system. These {{< math >}} $u$ {{< /math >}} points are required to evaluate the probability generating function (PGF) using inverse FFT (`irfft2`) later.

## Step-by-step Explanation

### 1. Efficient Real FFT Storage

```python
mx[-1] = mx[-1] // 2 + 1
```

This is for efficient storage of real FFTs in 2D. The real FFT (`irfft2`) only needs about half of the frequency bins in the last axis. It reduces the number of frequency components in the last dimension accordingly.

### 2. Dimension Iteration

```python
for i in range(len(mx))
```

Iterate over all state dimensions (e.g., mRNA, protein, etc.)

### 3. Index Vector Generation

```python
l = np.arange(mx[i])
```

Get the index vector {{< math >}} $l = [0, 1, \ldots, mx[i]-1]$ {{< /math >}} for that dimension.

### 4. Fourier Basis Computation

```python
u_ = np.exp(-2j * np.pi * l / lm[i]) - 1
```

This defines the Fourier basis values {{< math >}} $u\_$ {{< /math >}} for that dimension:

{{< math >}} 
$$u = e^{-2\pi i l/L} - 1 \tag{1}$$
{{< /math >}}

Where:
- {{< math >}} $e^{-2\pi i l/L}$ {{< /math >}} generates points on the unit circle (complex exponentials)
- Subtracting {{< math >}} $1$ {{< /math >}} gives shifted evaluation points used in the PGF formula

In generating function theory, the variable {{< math >}} $u$ {{< /math >}} is often taken as {{< math >}} $z-1$ {{< /math >}} where {{< math >}} $z \in \mathbb{C}$ {{< /math >}}, the complex unit circle. This is why the shift {{< math >}} $-1$ {{< /math >}} is done.

### 5. Storage

```python
u.append(u_)
```

Store these basis values in the list {{< math >}} $u$ {{< /math >}}, one per dimension (e.g., {{< math >}} $u[0]$ {{< /math >}} for mRNA, {{< math >}} $u[1]$ {{< /math >}} for protein).

## Why This Matters

These {{< math >}} $u$ {{< /math >}} values encode how the probability distribution will be transformed (via PGF). The function `g = np.meshgrid(...)` then builds all combinations {{< math >}} $(u_1, u_2)$ {{< /math >}}, used in the integrand.

## Concrete Numeric Example

### Setup
Assume:
- You are modeling 2 species (e.g., mRNA and protein)
- The maximum number of states in each species = {{< math >}} $lm = [5, 6]$ {{< /math >}}

Then {{< math >}} $mx = [5, 3]$ {{< /math >}} because of:

```python
mx[-1] = lm[-1] // 2 + 1  # mx[1] = 6//2 + 1 = 4
```

but in the code it is then redefined to 3 (could be manually changed to reduce cost).

### Step-by-step Computation (for i = 0, the mRNA dimension)

```python
l = np.arange(5)  # [0, 1, 2, 3, 4]
lm[0] = 5

u_0 = np.exp(-2j * np.pi * l / 5) - 1
```

Using Euler's formula:

{{< math >}} 
$$e^{-2\pi i l/5} = \cos(2\pi l/5) - i\sin(2\pi l/5) \tag{2}$$
{{< /math >}}

| {{< math >}} $l$ {{< /math >}} | {{< math >}} $e^{-2\pi i l/5}$ {{< /math >}} | {{< math >}} $u = e^{-2\pi i l/5} - 1$ {{< /math >}} |
|---|---|---|
| 0 | {{< math >}} $1 + 0i$ {{< /math >}} | {{< math >}} $0 + 0i$ {{< /math >}} |
| 1 | {{< math >}} $0.309 - 0.951i$ {{< /math >}} | {{< math >}} $-0.691 - 0.951i$ {{< /math >}} |
| 2 | {{< math >}} $-0.809 - 0.588i$ {{< /math >}} | {{< math >}} $-1.809 - 0.588i$ {{< /math >}} |
| 3 | {{< math >}} $-0.809 + 0.588i$ {{< /math >}} | {{< math >}} $-1.809 + 0.588i$ {{< /math >}} |
| 4 | {{< math >}} $0.309 + 0.951i$ {{< /math >}} | {{< math >}} $-0.691 + 0.951i$ {{< /math >}} |

So:

```python
u[0] = np.array([0 + 0j,
                 -0.691 - 0.951j,
                 -1.809 - 0.588j,
                 -1.809 + 0.588j,
                 -0.691 + 0.951j])
```

## Geometric View (Diagram)

```
Unit circle in complex plane (|z|=1)

             Im
              ↑
        *     |     *
              |  
  *-----------o-----------* Re
              |  
        *     |     *
              ↓

Legend:
- Points on the circle = exp(-2πil/L) for l = 0 to L-1
- Arrows to each point = u = exp(-2πil/L) - 1 (shifted from 1 to that point)
```

Each {{< math >}} $u\_$ {{< /math >}} is a vector from {{< math >}} $1$ {{< /math >}} (origin in generating function space) to a point on the unit circle.

## Use in CME

These complex-valued {{< math >}} $u$ {{< /math >}} points are used to evaluate the characteristic function or generating function (in `INTFUN`). The process involves:

{{< math >}} 
$$\text{Integrate over time} \rightarrow \text{Apply inverse FFT} \rightarrow \text{Convert to PMF} \tag{3}$$
{{< /math >}}

The inverse FFT converts from the Fourier domain back to the probability mass function (PMF) on integer state space.

# F) Integration Bounds T

```python
T = quad_vec_T * (1/bet + 1/gam + 1/kini)
```

This formula sets the upper limit of integration based on the natural timescales of the system:

| Parameter | Meaning | Timescale |
|-----------|---------|-----------|
| {{< math >}} $\beta$ {{< /math >}} | degradation rate of spliced mRNA | average lifetime ≈ {{< math >}} $1/\beta$ {{< /math >}} |
| {{< math >}} $\gamma$ {{< /math >}} | degradation/splicing rate of pre-mRNA | average lifetime ≈ {{< math >}} $1/\gamma$ {{< /math >}} |
| {{< math >}} $k_{\text{ini}}$ {{< /math >}} | transcription initiation rate | mean time between bursts ≈ {{< math >}} $1/k_{\text{ini}}$ {{< /math >}} |

Adding these together:

{{< math >}} 
$$\text{characteristic time} = \left(\frac{1}{\beta} + \frac{1}{\gamma} + \frac{1}{k_{\text{ini}}}\right) \tag{1}$$
{{< /math >}}

So you're multiplying this by a scalar {{< math >}} $\text{quad\_vec\_T}$ {{< /math >}} (default ∞, or 10 in fixed_quad mode) to make sure you integrate long enough for the system to approach steady state.

This is a heuristic upper bound — not exact, but it ensures you capture the full contribution to the generating function.

## Then What Happens?

```python
gf = scipy.integrate.quad_vec(fun, 0, T)[0]
```

- {{< math >}} $\text{fun}$ {{< /math >}} is a vector-valued function built to compute terms in the generating function.
- {{< math >}} $\text{quad\_vec}$ {{< /math >}} performs vectorized quadrature (numerical integration) over {{< math >}} $[0,T]$ {{< /math >}}.
- The {{< math >}} $[0]$ {{< /math >}} extracts the integral result (not the error estimate).

Afterward:

```python
gf = np.exp(kini * gf)
```

That gives you the generating function {{< math >}} $G(z_1, z_2)$ {{< /math >}} over the frequency grid.

## Summary

- {{< math >}} $T$ {{< /math >}} is a time cutoff that captures all relevant molecular events (transcription, splicing, degradation).
- {{< math >}} $\text{quad\_vec}(\text{fun}, 0, T)$ {{< /math >}} integrates the generating function terms over time.
- This integral is part of an analytical solution to the CME using generating functions.

## Ambiguous Code Needs Clarification

```python
gf = np.exp(kini * gf)
```

The initial {{< math >}} $\text{gf}$ {{< /math >}} contains values of the definite integral:

{{< math >}} 
$$\text{gf}[i,j] = \int_0^t \left(\frac{1}{1-b(z_1(s;i,j)-1)} - 1\right) ds \tag{2}$$
{{< /math >}}

for each Fourier mode {{< math >}} $(i, j)$ {{< /math >}} (i.e., each frequency pair).

Then, the line applies:

{{< math >}} 
$$G[i,j] = \exp(k_{\text{ini}} \cdot \text{that integral}) \tag{3}$$
{{< /math >}}

So it's implementing:

{{< math >}} 
$$G(z_1, z_2, t) = \exp\left(k_{\text{ini}} \int_0^t \left(\frac{1}{1-b(z_1(s)-1)} - 1\right) ds\right) \tag{4}$$
{{< /math >}}

**Note:** It assumes the initial condition is {{< math >}} $G_0 = 1$ {{< /math >}} (the generating function starts from 1 for all modes).

# G) Generating Functions and the Chemical Master Equation (CME)

Generating functions corresponding to solutions of the Chemical Master Equation (CME) can correspond to discrete nonnegative integer-valued distributions — but only under certain conditions.

## Comply Under Proper Conditions

When the CME is solved correctly and represents a proper stochastic process, the generating function does comply:

The CME describes the probability distribution over molecular counts (nonnegative integers).

If {{< math >}} $P(n,t)$ {{< /math >}} is the probability of having {{< math >}} $n$ {{< /math >}} molecules at time {{< math >}} $t$ {{< /math >}}, then the generating function is:

{{< math >}} 
$$G(z,t) = \sum_{n=0}^{\infty} P(n,t) z^n \tag{1}$$
{{< /math >}}

This is a valid probability generating function (PGF):

- Coefficients {{< math >}} $P(n,t) \geq 0$ {{< /math >}}
- {{< math >}} $\sum_n P(n,t) = 1$ {{< /math >}} (equation 2)
- {{< math >}} $n \in \mathbb{N}_0$ {{< /math >}}

✔ **So**: the inverse FFT can be used to recover {{< math >}} $P(n,t)$ {{< /math >}} from evaluations of {{< math >}} $G(z,t)$ {{< /math >}} on the unit circle.

## Not Always in Practice

In numerical or approximate solutions of the CME (e.g., moment closures, perturbative expansions, symbolic approximations), the generating function:

- Might have non-physical coefficients (e.g., negatives)
- Might not converge properly  
- Might not normalize to 1

So while the CME itself ensures that {{< math >}} $G(z,t)$ {{< /math >}} should correspond to a valid PMF, approximations or derivation mistakes may violate this.

## Also: Generating Function Types Matter

**Probability generating function (PGF)**: Used for CME, discrete + integer-valued

**Moment generating function (MGF) or Laplace transform**: Also used, but not invertible by FFT in the same way

## ✔ Summary

| Case | Does Generating Function Correspond to a PMF? | FFT applicable? |
|------|---------------------------------------------|-----------------|
| Exact solution to CME | ✅ Yes | ✅ Yes |
| Approximate CME (e.g., moment closure) | ⚠️ Not guaranteed | ❌ Maybe not |
| Symbolic/misused generating function | ❌ No | ❌ No |