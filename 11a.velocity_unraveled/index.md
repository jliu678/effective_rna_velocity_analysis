---
title: 🧬 Math Derivation of CME-defined Stochastic Model of RNA Velocity
summary: Stochastic model of RNA velocity defined by the Chemical Master Equation (CME) is superior to the deterministic ODE model in capturing the stochastic nature of single-cell RNA sequencing data. I derive the key equations from the paper [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492) and explain their biological significance.  
date: 2025-04-04
authors:
  - admin
tags:
  - scRNAseq,Stochastic Model,Differential Equation
  - probability generating function (PGF)
  - Chemical Master Equation (CME)
  - Method of Characteristics
  - Partial Differential Equation (PDE)
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Overview of Eq(18) from [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492)
The stochastic model defined by the Chemical Master Equation (CME) outperforms deterministic ODE models in capturing the inherent stochasticity of single-cell RNA sequencing (scRNA-seq) data. It is actively developed to provide a more accurate representation of feature counts and their underlying biological processes. And it has also enabled the generation of simulated data to evaluate deterministic ODE models and associated data processing methods commonly used in scRNA-seq analysis. Thus I derive the key equations from the paper [velocity unraveled](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492), a pivotal paper demonstrating the transformative potential of stochastic approaches. 

### 🧬 1. Probability generating function (PGF)

We define the state of a cell as:

{{< math >}} 
$$x = (x_u, x_s)$$ 
{{< /math >}}

Where:
- {{< math >}} $x_u$ {{< /math >}}: unspliced mRNA count
- {{< math >}} $x_s$ {{< /math >}}: spliced mRNA count

We now write the generating function of their joint distribution:

{{< math >}} 
$$G(u_u, u_s, t) = \sum_x P(x,t)(u_u + 1)^{x_u}(u_s + 1)^{x_s}$$ 
{{< /math >}}

This is a modified bivariate probability generating function, where the "+1" shift is standard in certain moment-generating setups. It lets you cleanly extract moments via derivatives of {{< math >}} $G$ {{< /math >}}.

### ⚙️ 2. Characteristic of ODEs derived from CME

{{< math >}} $$U_1(u_u, u_s, s) = \frac{u_s \beta}{\beta - \gamma} e^{-\gamma s} + \left(u_u - \frac{u_s \beta}{\beta - \gamma}\right) e^{-\beta s}$$ {{< /math >}}

This expression arises from [solving a linear system of ODEs for the chemical master equation (CME) via generating functions](#Derive-G), and it essentially encodes how the state propagates in time. It's derived from how unspliced → spliced reactions occur over time.

### 🧠 3. Log-generating function $f$<a id="log-generating-function-first-introduced"></a>


{{< math >}} 
$$f(u_u, u_s, t) := \ln G(u_u, u_s, t) = \int_0^t \alpha(t-s) U_1(u_u, u_s, s) ds$$ 
{{< /math >}}
This integral form is derived in [later section](#Derive-log-generating-function) and it tells us how the log of the generating function evolves, driven by transcription rate {{< math >}} $\alpha(t-s)$ {{< /math >}} and the system dynamics encoded in {{< math >}} $U_1$ {{< /math >}}. Essentially, it's the cumulative effect of production and conversion over time.

Then the log-GF can be written in a linear form in {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}}:

{{< math >}} 
$$f(u_u, u_s, t) = \mu_u(t) u_u + \mu_s(t) u_s$$ 
{{< /math >}}

This suggests that the process is governed by Poisson distributions, since the log-GF is linear in the arguments.

### 📊 4. Explicit distribution — product of Poissons

Given the above, we can now recover the joint distribution {{< math >}} $P(x,t)$ {{< /math >}}:

{{< math >}} 
$$P(x,t) = \frac{\mu_u(t)^{x_u} e^{-\mu_u(t)}}{x_u!} \cdot \frac{\mu_s(t)^{x_s} e^{-\mu_s(t)}}{x_s!}$$ 
{{< /math >}}

So:
- {{< math >}} $x_u \sim \text{Poisson}(\mu_u(t))$ {{< /math >}}
- {{< math >}} $x_s \sim \text{Poisson}(\mu_s(t))$ {{< /math >}}
- Jointly independent

This model assumes that given time {{< math >}} $t$ {{< /math >}}, the spliced and unspliced counts are independent Poisson-distributed variables, whose rates {{< math >}} $\mu_u(t), \mu_s(t)$ {{< /math >}} evolve in time according to the underlying biochemistry.

### 🧮 5. Time-averaged distribution

Finally, since single-cell sequencing samples cells asynchronously in time, the observed distribution over counts is not {{< math >}} $P(x,t)$ {{< /math >}} at a fixed {{< math >}} $t$ {{< /math >}}, but a time-averaged version:

{{< math >}} 
$$P(x) = \frac{1}{T} \int_0^T P(x,t) dt$$ 
{{< /math >}}

This is a mixture of Poissons over time, reflecting the asynchrony of cells in scRNA-seq snapshots. This averaging introduces overdispersion, which is critical to explain the variance observed in real data — greater than what a single Poisson can model.


## A) PGF Introduction

### 1. Multivariate Probability Generating Function (PGF)

This function encodes the entire joint distribution of the random vector {{< math >}} $x = (x_u, x_s)$ {{< /math >}}, i.e., the number of unspliced ({{< math >}} $x_u$ {{< /math >}}) and spliced ({{< math >}} $x_s$ {{< /math >}}) transcripts.

It's a bivariate generating function, meaning it's a function of two complex variables {{< math >}} $u_u$ {{< /math >}} and {{< math >}} $u_s$ {{< /math >}}.

The inclusion of {{< math >}} $+1$ {{< /math >}} makes it a shifted PGF, often done for technical convenience (especially when converting to moment-generating functions).

### 2. Moment Extraction

From the properties of generating functions:

**First Moments:**

{{< math >}} 
$$ \frac{\partial G}{\partial u_u}\bigg|_{u_u = u_s = 0} = E[x_u], \quad \frac{\partial G}{\partial u_s}\bigg|_{u_u = u_s = 0} = E[x_s] $$ 
{{< /math >}}

**Second Moments / Covariances:**

{{< math >}} 
$$ \frac{\partial^2 G}{\partial u_u^2}\bigg|_{u_u = u_s = 0} = E[x_u(x_u - 1)], \quad \frac{\partial^2 G}{\partial u_u \partial u_s}\bigg|_{u_u = u_s = 0} = E[x_u x_s] $$ 
{{< /math >}}

These derivatives allow us to compute variances, covariances, and higher-order statistics of transcript counts.

### 3. Connection to Biological Reactions

In a linear RNA kinetic model:

- {{< math >}} $x_u$ {{< /math >}}: produced at rate {{< math >}} $\alpha$ {{< /math >}}, converted to {{< math >}} $x_s$ {{< /math >}} at rate {{< math >}} $\beta$ {{< /math >}}
- {{< math >}} $x_s$ {{< /math >}}: degraded at rate {{< math >}} $\gamma$ {{< /math >}}

The evolution of {{< math >}} $G(u_u, u_s, t)$ {{< /math >}} over time follows a partial differential equation that arises from the Chemical Master Equation (CME), which governs the time evolution of probability distributions in chemical kinetics.

### 4. Time-Dependence

{{< math >}} $G$ {{< /math >}} is explicitly time-dependent, evolving as the distribution {{< math >}} $P(x,t)$ {{< /math >}} changes.

In some derivations, the log-generating function {{< math >}} $f = \log G$ {{< /math >}} is linear in {{< math >}} $u_u, u_s$ {{< /math >}}, which implies that {{< math >}} $x_u, x_s \sim \text{Poisson}(\mu_u(t)), \text{Poisson}(\mu_s(t))$ {{< /math >}} and are independent.

### 5. Decay and Stationarity

As {{< math >}} $t \to \infty$ {{< /math >}}:

- The means {{< math >}} $\mu_u(t), \mu_s(t)$ {{< /math >}} stabilize.
- So does the generating function, converging to that of a product of Poisson distributions (one for each transcript species).

### Summary: Why It Matters

- **Encodes all statistical information** about {{< math >}} $x_u, x_s$ {{< /math >}}
- **Enables exact computation** of moments, cumulants, and correlations
- **Links stochastic biochemical kinetics** with observed scRNA-seq distributions
- **Supports modeling** of noise, burstiness, and cell-to-cell heterogeneity


## B) Derive $u_u(s)$ via method of characteristics <a id="Derive-G"></a>

{{< math >}} $U_1(u_u, u_s, s)$ {{< /math >}} is claimed to be {{< math >}} $u_u(s)$ {{< /math >}}, which is the solution to the ODE derived from the Chemical Master Equation (CME) for the two-species birth-death process representing unspliced (u) and spliced (s) mRNA dynamics.

{{< math >}} $$U_1(u_u, u_s, s) := u_u(s) = \frac{u_s \beta}{\beta - \gamma} e^{-\gamma s} + \left(u_u - \frac{u_s \beta}{\beta - \gamma}\right) e^{-\beta s}$$ {{< /math >}}

This derivation uses the method of characteristics applied to the generating function of a stochastic process governed by the Chemical Master Equation (CME). We consider a two-species birth-death process representing unspliced (u) and spliced (s) mRNA dynamics.

### Step 1: Define the Generating Function

Let the joint probability distribution of unspliced and spliced mRNA at time {{< math >}} $t$ {{< /math >}} be {{< math >}} $P(x_u, x_s, t)$ {{< /math >}}. The generating function is:

{{< math >}} $$G(z_u, z_s, t) = \sum_{x_u, x_s} P(x_u, x_s, t) z_u^{x_u} z_s^{x_s}$$ {{< /math >}}

We define new variables:

{{< math >}} $$z_u = u_u + 1, \quad z_s = u_s + 1$$ {{< /math >}}

so the generating function becomes:

{{< math >}} $$G(u_u, u_s, t) = \sum_{x_u, x_s} P(x_u, x_s, t) (u_u + 1)^{x_u} (u_s + 1)^{x_s}$$ {{< /math >}}

### Step 2: CME and Corresponding PDE

The CME for this system is governed by the reactions:

- Transcription (birth of unspliced): rate {{< math >}} $\alpha$ {{< /math >}}
- Splicing: {{< math >}} $u \xrightarrow{\beta} s$ {{< /math >}}
- Degradation: {{< math >}} $s \xrightarrow{\gamma} \emptyset$ {{< /math >}}

From the CME, the PDE for {{< math >}} $G$ {{< /math >}} is below as derived in [another section](#Derive-step-3-of-B)):

{{< math >}} $$\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s}$$ {{< /math >}}

### Step 3: Method of Characteristics (turn PDE into ODEs)

We now transform this PDE into ODEs using the method of characteristics. Let {{< math >}} $u_u(s), u_s(s), G(s)$ {{< /math >}} be functions of characteristic time {{< math >}} $s$ {{< /math >}} such that along these paths:

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s - u_u), \quad \frac{du_s}{ds} = -\gamma u_s, \quad \frac{dG}{ds} = \alpha u_u G$$ {{< /math >}}

Let's solve these:

### Step 4: Solve for $u_s(s)$

{{< math >}} $$\frac{du_s}{ds} = -\gamma u_s \Rightarrow u_s(s) = u_s(0) e^{-\gamma s}$$ {{< /math >}}

### Step 5: Solve for $u_u(s)$

Use integrating factor method:

{{< math >}} $$\frac{du_u}{ds} + \beta u_u = \beta u_s(s) = \beta u_s(0) e^{-\gamma s}$$ {{< /math >}}

Multiply both sides by {{< math >}} $e^{\beta s}$ {{< /math >}}:

{{< math >}} $$\frac{d}{ds}(u_u e^{\beta s}) = \beta u_s(0) e^{(\beta - \gamma)s}$$ {{< /math >}}

Integrate:

{{< math >}} $$u_u(s) e^{\beta s} = u_u(0) + \frac{\beta u_s(0)}{\beta - \gamma}(e^{(\beta - \gamma)s} - 1)$$ {{< /math >}}

Solve for {{< math >}} $u_u(s)$ {{< /math >}}:

{{< math >}} $$u_u(s) = u_u(0) e^{-\beta s} + \frac{\beta u_s(0)}{\beta - \gamma}(e^{-\gamma s} - e^{-\beta s})$$ {{< /math >}}

### Step 6: Define $U_1(u_u, u_s, s)$

We identify:

{{< math >}} $$U_1(u_u, u_s, s) := u_u(s) = \frac{u_s \beta}{\beta - \gamma} e^{-\gamma s} + \left(u_u - \frac{u_s \beta}{\beta - \gamma}\right) e^{-\beta s}$$ {{< /math >}}

## B.1) Derive PDE from CME (step 2 of B)

### Step 1: Write the CME explicitly
<a id="Derive-step-3-of-B"></a>
The CME for the joint distribution {{< math >}} $P(x_u, x_s, t)$ {{< /math >}} of the unspliced {{< math >}} $x_u$ {{< /math >}} and spliced {{< math >}} $x_s$ {{< /math >}} RNA is:

{{< math >}}
$$
\frac{d}{dt}P(x_u, x_s, t) = \alpha[P(x_u - 1, x_s, t) - P(x_u, x_s, t)] + \beta[(x_u + 1)P(x_u + 1, x_s - 1, t) - x_u P(x_u, x_s, t)] + \gamma[(x_s + 1)P(x_u, x_s + 1, t) - x_s P(x_u, x_s, t)] \qquad (9)
$$
{{< /math >}}

where:
- {{< math >}} $\alpha$ {{< /math >}}: production rate of unspliced RNA
- {{< math >}} $\beta$ {{< /math >}}: splicing rate (unspliced → spliced)
- {{< math >}} $\gamma$ {{< /math >}}: degradation rate of spliced RNA

### Step 2: Define the generating function $G(u_u, u_s, t)$

{{< math >}}
$$
G(u_u, u_s, t) := \sum_{x_u=0}^{\infty} \sum_{x_s=0}^{\infty} P(x_u, x_s, t)(u_u + 1)^{x_u}(u_s + 1)^{x_s} \qquad (10)
$$
{{< /math >}}

**Note:** Using {{< math >}} $(u_u + 1)^{x_u}$ {{< /math >}} instead of {{< math >}} $u_u^{x_u}$ {{< /math >}} is a common shift to simplify derivatives later, but you can equivalently define with {{< math >}} $u_u^{x_u}$ {{< /math >}}.

### Step 3: Take time derivative of $G$

Using linearity of sums and derivatives:

{{< math >}}
$$
\frac{\partial G}{\partial t} = \sum_{x_u=0}^{\infty} \sum_{x_s=0}^{\infty} \frac{\partial P(x_u, x_s, t)}{\partial t} (u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Substitute the CME expression:

{{< math >}}
$$
= \sum_{x_u,x_s} [\alpha(P(x_u - 1, x_s, t) - P(x_u, x_s, t)) + \beta((x_u + 1)P(x_u + 1, x_s - 1, t) - x_u P(x_u, x_s, t)) + \gamma((x_s + 1)P(x_u, x_s + 1, t) - x_s P(x_u, x_s, t))](u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

### Step 4: Evaluate each term separately

#### Term 1: Transcription

{{< math >}}
$$
\sum_{x_u,x_s} \alpha(P(x_u - 1, x_s) - P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Rewrite the first sum by shifting {{< math >}} $x_u \to x_u + 1$ {{< /math >}} in the first part:

{{< math >}}
$$
\sum_{x_u=0}^{\infty} P(x_u - 1, x_s)(u_u + 1)^{x_u} = \sum_{x_u'=-1}^{\infty} P(x_u', x_s)(u_u + 1)^{x_u' + 1}
$$
{{< /math >}}

Since {{< math >}} $P(x_u', x_s) = 0$ {{< /math >}} for {{< math >}} $x_u' < 0$ {{< /math >}}, this becomes:

{{< math >}}
$$
(u_u + 1) \sum_{x_u'=0}^{\infty} P(x_u', x_s)(u_u + 1)^{x_u'}
$$
{{< /math >}}

Therefore,

{{< math >}}
$$
\sum_{x_u,x_s} \alpha(P(x_u - 1, x_s) - P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s} = \alpha((u_u + 1)G - G) = \alpha u_u G
$$
{{< /math >}}

#### Term 2: Splicing

{{< math >}}
$$
\sum_{x_u,x_s} \beta((x_u + 1)P(x_u + 1, x_s - 1) - x_u P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Split into two sums:

{{< math >}}
$$
S_1 = \beta \sum_{x_u,x_s} (x_u + 1)P(x_u + 1, x_s - 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

{{< math >}}
$$
S_2 = -\beta \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Change indices in {{< math >}} $S_1$ {{< /math >}}:

Let {{< math >}} $x_u' = x_u + 1$ {{< /math >}}, {{< math >}} $x_s' = x_s - 1$ {{< /math >}}, so {{< math >}} $x_u = x_u' - 1$ {{< /math >}}, {{< math >}} $x_s = x_s' + 1$ {{< /math >}}

Then,

{{< math >}}
$$
S_1 = \beta \sum_{x_u'=1}^{\infty} \sum_{x_s'=0}^{\infty} x_u' P(x_u', x_s')(u_u + 1)^{x_u' - 1}(u_s + 1)^{x_s' + 1}
$$
{{< /math >}}

Rearranged:<a id="why_sum_is_equivalence"></a>
 (see [note](#b2-why-double-sum-equivalence-holds-in-term-2-of-step4-in-b1))
 
{{< math >}}
$$
= \beta(u_s + 1) \sum_{x_u',x_s'} x_u' P(x_u', x_s')(u_u + 1)^{x_u' - 1}(u_s + 1)^{x_s'}
$$
{{< /math >}}

{{< math >}} $S_2$ {{< /math >}} is:

{{< math >}}
$$
S_2 = -\beta \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Now recognize:

{{< math >}}
$$
\frac{\partial G}{\partial u_u} = \sum_{x_u,x_s} x_u P(x_u, x_s)(u_u + 1)^{x_u - 1}(u_s + 1)^{x_s}
$$
{{< /math >}}

So:

{{< math >}}
$$
S_1 = \beta(u_s + 1) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

and

{{< math >}}
$$
S_2 = -\beta(u_u + 1) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

Putting together:

{{< math >}}
$$
\text{Splicing term} = \beta((u_s + 1) - (u_u + 1)) \frac{\partial G}{\partial u_u} = \beta(u_s - u_u) \frac{\partial G}{\partial u_u}
$$
{{< /math >}}

#### Term 3: Degradation

{{< math >}}
$$
\sum_{x_u,x_s} \gamma((x_s + 1)P(x_u, x_s + 1) - x_s P(x_u, x_s))(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

Split into:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u,x_s} (x_s + 1)P(x_u, x_s + 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

{{< math >}}
$$
S_4 = -\gamma \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

We have:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u,x_s} (x_s + 1)P(x_u, x_s + 1)(u_u + 1)^{x_u}(u_s + 1)^{x_s}
$$
{{< /math >}}

With the substitution {{< math >}} $x_s' = x_s + 1$ {{< /math >}}, we get {{< math >}} $x_s = x_s' - 1$ {{< /math >}} and:

{{< math >}}
$$
S_3 = \gamma \sum_{x_u=0}^{\infty} \sum_{x_s'=1}^{\infty} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s' - 1}
$$
{{< /math >}}

{{< math >}}
$$
= \frac{\gamma}{u_s + 1} \sum_{x_u,x_s'} x_s' P(x_u, x_s')(u_u + 1)^{x_u}(u_s + 1)^{x_s'}
$$
{{< /math >}}

Recognize

{{< math >}}
$$
\frac{\partial G}{\partial u_s} = \sum_{x_u,x_s} x_s P(x_u, x_s)(u_u + 1)^{x_u}(u_s + 1)^{x_s - 1}
$$
{{< /math >}}
{{< math >}}
$$
= \frac{\gamma}{u_s + 1} \cdot (u_s + 1) \frac{\partial G}{\partial u_s} = \gamma \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

And:

{{< math >}}
$$
S_4 = -\gamma(u_s + 1) \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

So the degradation term is:

{{< math >}}
$$
\text{Degradation term} = \gamma \frac{\partial G}{\partial u_s} - \gamma(u_s + 1) \frac{\partial G}{\partial u_s} = -\gamma u_s \frac{\partial G}{\partial u_s}
$$
{{< /math >}}

### Final Result

Combining all terms:

{{< math >}}
$$
\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u) \frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s} \qquad (11)
$$
{{< /math >}}

## B.2) Why Double Sum Equivalence Holds [in Term 2 of Step4 in B.1)](#why_sum_is_equivalence)

This equivalence depends on what values the indices range over and how the function being summed behaves.

### 1. Notation

When we write:

{{< math >}} $$\sum_{x_u', x_s'} f(x_u', x_s')$$ {{< /math >}}

This is shorthand for:

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty} f(x_u', x_s')$$ {{< /math >}}

That is, summing over all nonnegative integer pairs {{< math >}} $(x_u', x_s') \in \mathbb{N}_0 \times \mathbb{N}_0$ {{< /math >}}.

### 2. 'Suspect' in the derivation

In our original sum:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') (u_u + 1)^{x_u' - 1} (u_s + 1)^{x_s'}$$ {{< /math >}}

The lower bound of {{< math >}} $x_u'$ {{< /math >}} is 1 because we performed a change of variables from {{< math >}} $x_u = x_u' - 1$ {{< /math >}}, and in the original sum, {{< math >}} $x_u \geq 0$ {{< /math >}}, which implies {{< math >}} $x_u' \geq 1$ {{< /math >}}.

So the double sum with bounds:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

is not exactly the same as:

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

But if we define:

{{< math >}} $$\sum_{x_u', x_s'} := \sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

then in your derivation, the support of {{< math >}} $P(x_u', x_s')$ {{< /math >}} makes this safe because:

{{< math >}} $$x_u' P(x_u', x_s') = 0 \text{ when } x_u' = 0$$ {{< /math >}}

since the factor is 0.

So extending the lower limit to 0 adds no contribution to the sum.

### ✅ Conclusion

The expressions:

{{< math >}} $$\sum_{x_u' = 1}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') \cdots$$ {{< /math >}}

and

{{< math >}} $$\sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty} x_u' P(x_u', x_s') \cdots$$ {{< /math >}}

are equal because when {{< math >}} $x_u' = 0$ {{< /math >}}, the term is 0.

Hence, we can write:

{{< math >}} $$\sum_{x_u', x_s'} := \sum_{x_u' = 0}^{\infty} \sum_{x_s' = 0}^{\infty}$$ {{< /math >}}

without affecting the value of the sum.

## B.3) Solve PDE via Method of Characteristics

We start with the PDE derived from the chemical master equation (CME) for a stochastic model of unspliced (u) and spliced (s) RNA:

{{< math >}} $$\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s}$$ {{< /math >}}

This is a first-order linear PDE in 3 variables: {{< math >}} $u_u, u_s, t$ {{< /math >}}.

To solve this, we apply the method of characteristics, which reduces a PDE to a system of ODEs along special curves (characteristics) in the domain {{< math >}} $(u_u, u_s, t)$ {{< /math >}}. The idea is to track how {{< math >}} $G$ {{< /math >}} changes along these curves as we change a parameter {{< math >}} $s$ {{< /math >}} (which can be thought of like an artificial time).

### Step 1: Define Characteristic Curves

We introduce {{< math >}} $s$ {{< /math >}} as a parameter along a characteristic curve and define:

{{< math >}} $$u_u = u_u(s)$$ {{< /math >}}

{{< math >}} $$u_s = u_s(s)$$ {{< /math >}}

{{< math >}} $$t = t(s)$$ {{< /math >}}

{{< math >}} $$G = G(u_u(s), u_s(s), t(s))$$ {{< /math >}}

Then the total derivative of {{< math >}} $G$ {{< /math >}} along the curve is:

{{< math >}} $$\frac{dG}{ds} = \frac{\partial G}{\partial u_u}\frac{du_u}{ds} + \frac{\partial G}{\partial u_s}\frac{du_s}{ds} + \frac{\partial G}{\partial t}\frac{dt}{ds}$$ {{< /math >}}

Now, we substitute the PDE into this expression. From the PDE:

{{< math >}} $$\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s}$$ {{< /math >}}

Plugging this into the total derivative:

{{< math >}} $$\frac{dG}{ds} = \frac{\partial G}{\partial u_u}\frac{du_u}{ds} + \frac{\partial G}{\partial u_s}\frac{du_s}{ds} + \left[\alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s}\right]\frac{dt}{ds}$$ {{< /math >}}

Now, **choose {{< math >}} $\frac{dt}{ds} = 1$ {{< /math >}}**. This simplifies the expression because now {{< math >}} $t = s$ {{< /math >}}, and we can reduce the 3-variable PDE into a system of ODEs in {{< math >}} $s$ {{< /math >}}.

### Step 2: Match Terms

To make the right-hand side cancel cleanly, we group terms:

**Coefficient of** {{< math >}} $\frac{\partial G}{\partial u_u}$ {{< /math >}}:

{{< math >}} $$\frac{du_u}{ds} + \beta(u_s - u_u)$$ {{< /math >}}

**Coefficient of** {{< math >}} $\frac{\partial G}{\partial u_s}$ {{< /math >}}:

{{< math >}} $$\frac{du_s}{ds} - \gamma u_s$$ {{< /math >}}

To cancel the dependence on {{< math >}} $\frac{\partial G}{\partial u_u}$ {{< /math >}} and {{< math >}} $\frac{\partial G}{\partial u_s}$ {{< /math >}}, we set these to zero, yielding:

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s - u_u), \quad \frac{du_s}{ds} = -\gamma u_s$$ {{< /math >}}

Then the remaining term becomes:

{{< math >}} $$\frac{dG}{ds} = \alpha u_u G$$ {{< /math >}}

Now we've reduced the PDE into this system of ODEs:

{{< math >}} $$\begin{align}
\frac{du_s}{ds} &= -\gamma u_s \\
\frac{du_u}{ds} &= \beta(u_s - u_u) \\
\frac{dG}{ds} &= \alpha u_u G
\end{align}$$ {{< /math >}}

These are much easier to solve analytically or numerically.

## B.4) Legitimacy to choose $\frac{dt}{ds} = 1$ in B.3)

We choose {{< math >}} $\frac{dt}{ds} = 1$ {{< /math >}} in the method of characteristics because it simplifies the partial differential equation (PDE) to a more tractable set of ordinary differential equations (ODEs) — and this choice is completely valid and standard in this method. Let me explain why:

### What Does $\frac{dt}{ds} = 1$ Mean

In the method of characteristics, you introduce a new parameter {{< math >}} $s$ {{< /math >}} that traces out a path (or characteristic curve) in the space of independent variables — here, {{< math >}} $(u_u, u_s, t)$ {{< /math >}}. Along this path:

{{< math >}} $$\frac{dG}{ds} = \frac{\partial G}{\partial u_u} \frac{du_u}{ds} + \frac{\partial G}{\partial u_s} \frac{du_s}{ds} + \frac{\partial G}{\partial t} \frac{dt}{ds}$$ {{< /math >}}

So {{< math >}} $\frac{dt}{ds}$ {{< /math >}} tells you how "fast" you're moving in the time direction along the characteristic curve.

Now, if you choose:

{{< math >}} $$\frac{dt}{ds} = 1 \Rightarrow t = s$$ {{< /math >}}

you're saying: "Let the parameter along the path simply equal time." This simplifies the math without changing the problem.

### Why Is This Legitimate

Because in the method of characteristics, {{< math >}} $s$ {{< /math >}} is a dummy variable. You're free to choose how it relates to the original coordinates, as long as it parametrizes a valid path. The PDE solution is determined by the behavior along these characteristic curves, and the parameterization does not affect the final solution.

Choosing {{< math >}} $\frac{dt}{ds} = 1$ {{< /math >}}:

- simplifies the system (fewer variables)
- allows you to think of {{< math >}} $s$ {{< /math >}} as time
- turns the original PDE into a solvable system of ODEs
- does not change the physics or solution of the underlying system

It's not a constraint imposed by the model — it's a strategic mathematical choice to simplify the derivation.

### Alternative Choices

Yes — you could, in principle, choose something else (e.g. {{< math >}} $\frac{dt}{ds} = \gamma$ {{< /math >}}, or make {{< math >}} $s = -t$ {{< /math >}}), but then the resulting ODEs are messier. Since the method is agnostic to how you parametrize the path, you're free to choose the one that leads to the simplest math.

### Summary

Choosing {{< math >}} $\frac{dt}{ds} = 1$ {{< /math >}} is valid because:

1. The method of characteristics allows any parameterization.
2. This choice makes the equations easier to solve.
3. The physics or stochastic model (e.g. CME) remains unchanged.

## C) Derive log-generating function <a id="Derive-log-generating-function"></a>
Below derive the log-generating function first introduced [here](#log-generating-function-first-introduced).

### 1. Setup: Use the Method of Characteristics

We rewrite the PDE:

{{< math >}} $$\frac{\partial G}{\partial t} = \alpha u_u G + \beta(u_s - u_u)\frac{\partial G}{\partial u_u} - \gamma u_s \frac{\partial G}{\partial u_s}$$ {{< /math >}}

This is a linear PDE in {{< math >}} $G(u_u, u_s, t)$ {{< /math >}}, and we apply the method of characteristics.

Let {{< math >}} $s$ {{< /math >}} be the parameter along characteristic curves. Then we solve the system:

{{< math >}} $$\frac{dt}{ds} = 1$$ {{< /math >}}

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s - u_u)$$ {{< /math >}}

{{< math >}} $$\frac{du_s}{ds} = -\gamma u_s$$ {{< /math >}}

{{< math >}} $$\frac{dG}{ds} = \alpha u_u G$$ {{< /math >}}

### 2. Solve the ODE for $u_s(s)$

{{< math >}} $$\frac{du_s}{ds} = -\gamma u_s \Rightarrow u_s(s) = u_s(0)e^{-\gamma s}$$ {{< /math >}}

### 3. Plug into the ODE for $u_u(s)$

{{< math >}} $$\frac{du_u}{ds} = \beta(u_s(s) - u_u(s)) = \beta(u_s(0)e^{-\gamma s} - u_u(s))$$ {{< /math >}}

This is a linear non-homogeneous ODE:

Let's solve using integrating factor:

Integrating factor: {{< math >}} $\mu(s) = e^{\beta s}$ {{< /math >}}

Multiply both sides:

{{< math >}} $$e^{\beta s}\frac{du_u}{ds} + \beta e^{\beta s}u_u(s) = \beta u_s(0)e^{(\beta - \gamma)s}$$ {{< /math >}}

{{< math >}} $$\frac{d}{ds}(e^{\beta s}u_u(s)) = \beta u_s(0)e^{(\beta - \gamma)s}$$ {{< /math >}}

Integrate both sides:

{{< math >}} $$e^{\beta s}u_u(s) = \frac{\beta u_s(0)}{\beta - \gamma}e^{(\beta - \gamma)s} + C$$ {{< /math >}}

Now divide both sides:

{{< math >}} $$u_u(s) = \frac{\beta u_s(0)}{\beta - \gamma}e^{-\gamma s} + Ce^{-\beta s}$$ {{< /math >}}

Apply initial condition {{< math >}} $u_u(0)$ {{< /math >}} to solve for {{< math >}} $C$ {{< /math >}}:

{{< math >}} $$u_u(0) = \frac{\beta u_s(0)}{\beta - \gamma} + C \Rightarrow C = u_u(0) - \frac{\beta u_s(0)}{\beta - \gamma}$$ {{< /math >}}

So we now have:

{{< math >}} $$u_u(s) = \frac{\beta u_s(0)}{\beta - \gamma}e^{-\gamma s} + \left(u_u(0) - \frac{\beta u_s(0)}{\beta - \gamma}\right)e^{-\beta s}$$ {{< /math >}}

### 4. Define $U_1(u_u, u_s, s)$

Recall that in the characteristic solution for {{< math >}} $G$ {{< /math >}}, we solve:

{{< math >}} $$\frac{dG}{ds} = \alpha u_u(s)G \Rightarrow G(s) = \exp\left(\int_0^t \alpha(t-s)u_u(s)ds\right)$$ {{< /math >}}

Define {{< math >}} $U_1(u_u, u_s, s) = u_u(s)$ {{< /math >}}. So:

{{< math >}} $$U_1(u_u, u_s, s) = \frac{\beta u_s}{\beta - \gamma}e^{-\gamma s} + \left(u_u - \frac{\beta u_s}{\beta - \gamma}\right)e^{-\beta s}$$ {{< /math >}}

### 5. Compute $f(u_u, u_s, t) = \ln G$

Now integrate:

{{< math >}} $$\ln G(u_u, u_s, t) = \int_0^t \alpha(t-s)U_1(u_u, u_s, s) \, ds$$ {{< /math >}}