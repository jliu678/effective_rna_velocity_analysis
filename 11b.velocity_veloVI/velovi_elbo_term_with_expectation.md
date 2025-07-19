---
title: "Velovi Elbo Term With Expectation"
date: 2025-07-10
draft: True
---

# VeloVI ELBO Third Term with Expectation

## The Expression

{{< math >}} $$E_{q_\phi(z_1|u_1,s_1)} \left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right]$$ {{< /math >}}

This is the third term of the ELBO for a single cell n (in your example, n=1), as seen in Equation (26) of the VeloVI model description.

Let's dissect each part:

### {{< math >}} $z_1$ {{< /math >}}
This represents the low-dimensional latent variable for Cell 1. It's a vector that summarizes the latent state of Cell 1 (e.g., its overall position in the developmental trajectory or cell cycle).

### {{< math >}} $u_1, s_1$ {{< /math >}}
These are the observed unspliced and spliced transcript abundances for Cell 1, across all genes. So, {{< math >}} $u_1$ {{< /math >}} is a vector {{< math >}} $(u_{1A}, u_{1B}, u_{1C})$ {{< /math >}} and {{< math >}} $s_1$ {{< /math >}} is {{< math >}} $(s_{1A}, s_{1B}, s_{1C})$ {{< /math >}} in your three-gene example.

### {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}
This is the approximate posterior distribution over the latent variable {{< math >}} $z_1$ {{< /math >}}, conditioned on the observed data {{< math >}} $(u_1,s_1)$ {{< /math >}} for Cell 1.

- It's part of the variational posterior {{< math >}} $q_\phi(z,\pi|u,s)$ {{< /math >}} (Eq. 22).
- The subscript {{< math >}} $\phi$ {{< /math >}} indicates that this distribution is parameterized by the neural network parameters {{< math >}} $\phi$ {{< /math >}}. Specifically, an "encoder" neural network takes {{< math >}} $(u_1,s_1)$ {{< /math >}} as input and outputs the parameters (e.g., mean and variance if it's a Gaussian) for this distribution.

### {{< math >}} $E_{q_\phi(z_1|u_1,s_1)}[\ldots]$ {{< /math >}}
This denotes an expectation with respect to the distribution {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}.

- In variational inference, we often can't compute this expectation in closed form.
- Instead, it's typically approximated using Monte Carlo sampling. We sample one or more {{< math >}} $z_1^*$ {{< /math >}} values from {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}} and then compute the term inside the expectation for each sample, averaging the results.
- The "reparameterization trick" is often used to make this sampling process differentiable, allowing gradients to flow back through the sampling process to update {{< math >}} $\phi$ {{< /math >}}.

### {{< math >}} $\sum_g$ {{< /math >}}
This simply means a summation over all genes g (Gene A, Gene B, Gene C in your example).

### {{< math >}} $\pi_{1g}$ {{< /math >}}
This is the probability distribution over the four transcriptional states for gene g in cell 1 (i.e., {{< math >}} $[\pi_{1g,k=1}, \pi_{1g,k=2}, \pi_{1g,k=3}, \pi_{1g,k=4}]$ {{< /math >}}).

### {{< math >}} $q_\phi(\pi_{1g}|z_1)$ {{< /math >}}
This is the approximate posterior distribution over the state assignment probabilities {{< math >}} $\pi_{1g}$ {{< /math >}}, conditioned on the latent variable {{< math >}} $z_1$ {{< /math >}}.

- This is also part of the variational posterior (Eq. 22).
- Another neural network (part of {{< math >}} $\phi$ {{< /math >}}) takes {{< math >}} $z_1$ {{< /math >}} as input and outputs the parameters (e.g., concentration parameters for a Dirichlet distribution) for this distribution.

### {{< math >}} $p(\pi_{1g})$ {{< /math >}}
This is the prior distribution over the state assignment probabilities {{< math >}} $\pi_{1g}$ {{< /math >}}.

As specified in Eq. 14, it's a uniform Dirichlet distribution: {{< math >}} $p(\pi_{1g}) = \text{Dirichlet}(0.25, 0.25, 0.25, 0.25)$ {{< /math >}}. This means, a priori, all four states are considered equally likely for any gene in any cell.

### {{< math >}} $\text{KL}(Q \| P)$ {{< /math >}}
This is the Kullback-Leibler (KL) divergence between two probability distributions Q and P. It's a measure of how one probability distribution Q diverges from a second, expected probability distribution P.

- {{< math >}} $\text{KL}(Q \| P) \geq 0$ {{< /math >}}. It is 0 if and only if {{< math >}} $Q = P$ {{< /math >}}.
- In this context, it measures how much our inferred posterior {{< math >}} $q_\phi(\pi_{1g}|z_1)$ {{< /math >}} deviates from the prior {{< math >}} $p(\pi_{1g})$ {{< /math >}}.

## Putting it all together

The term {{< math >}} $E_{q_\phi(z_1|u_1,s_1)} \left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right]$ {{< /math >}} represents the expected sum of KL divergences between the approximate posterior distribution over state assignments ({{< math >}} $\pi_{1g}$ {{< /math >}}) and its prior, where the expectation is taken over the approximate posterior distribution of the cell's latent variable ({{< math >}} $z_1$ {{< /math >}}).

---

# Understanding Expectations in Variational Inference

The phrase "expectation with respect to the distribution {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}" refers to the average value of a function or random variable, where the averaging is performed according to the probabilities defined by the distribution {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}.

Let's break it down in the context of VeloVI:

## {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}

- This is the approximate posterior distribution of the latent variable {{< math >}} $z_1$ {{< /math >}} (for Cell 1), given the observed unspliced ({{< math >}} $u_1$ {{< /math >}}) and spliced ({{< math >}} $s_1$ {{< /math >}}) expression data for that cell.
- It's an "approximate" posterior because, in variational inference, we cannot directly calculate the true posterior {{< math >}} $p(z_1|u_1,s_1)$ {{< /math >}}. Instead, we use a simpler, parameterized distribution (in this case, often a Gaussian) to approximate it.
- The subscript {{< math >}} $\phi$ {{< /math >}} indicates that the parameters of this distribution (e.g., its mean and variance if it's a Gaussian) are determined by a neural network (an "encoder") that takes the observed data {{< math >}} $(u_1,s_1)$ {{< /math >}} as input.

## Expectation (E)

In probability theory, the expectation of a random variable X with respect to a probability density function {{< math >}} $f(x)$ {{< /math >}} is denoted {{< math >}} $E_{f(x)}[X]$ {{< /math >}} or {{< math >}} $E[X]$ {{< /math >}}.

- For a continuous random variable X with PDF {{< math >}} $f(x)$ {{< /math >}}, the expectation of a function {{< math >}} $h(X)$ {{< /math >}} is calculated as an integral:

{{< math >}} $$E_{f(x)}[h(X)] = \int h(x)f(x)dx$$ {{< /math >}}

- For a discrete random variable X with PMF {{< math >}} $P(x)$ {{< /math >}}, it's a sum:

{{< math >}} $$E_{P(x)}[h(X)] = \sum_x h(x)P(x)$$ {{< /math >}}

## What is being averaged?

In the ELBO term you referenced, {{< math >}} $E_{q_\phi(z_1|u_1,s_1)} \left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right]$ {{< /math >}}, the function being averaged is {{< math >}} $\left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right]$ {{< /math >}}.

- Notice that this function depends on {{< math >}} $z_1$ {{< /math >}}. The KL divergence calculation for {{< math >}} $\pi_{1g}$ {{< /math >}} requires {{< math >}} $z_1$ {{< /math >}} as an input to {{< math >}} $q_\phi(\pi_{1g}|z_1)$ {{< /math >}}.
- Since {{< math >}} $z_1$ {{< /math >}} itself is a random variable drawn from {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}, the entire expression {{< math >}} $\left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right]$ {{< /math >}} is also a random variable whose value depends on the specific {{< math >}} $z_1$ {{< /math >}} sampled. The expectation averages this quantity over all possible {{< math >}} $z_1$ {{< /math >}} values, weighted by their probabilities according to {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}.

## Why is it necessary?

The ELBO is defined as an expectation because the true posterior distribution (which is what we ideally want to work with) is intractable. By taking an expectation with respect to the approximate posterior, we are essentially trying to make our approximation {{< math >}} $q_\phi$ {{< /math >}} as close as possible to the true posterior, on average, across all possible values of {{< math >}} $z_1$ {{< /math >}} that are deemed probable by {{< math >}} $q_\phi$ {{< /math >}}.

## How is it calculated in practice (Monte Carlo Approximation)

Since the integral (for continuous {{< math >}} $z_1$ {{< /math >}}) is usually intractable, this expectation is approximated using Monte Carlo sampling. The process is as follows:

### 1. Sample from {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}
Draw one or more samples, say {{< math >}} $z_1^{(m)}$ {{< /math >}} for {{< math >}} $m=1,\ldots,M$ {{< /math >}}, from the approximate posterior distribution {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}.

- If {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}} is a Gaussian, you sample from a Gaussian.
- The "reparameterization trick" is often used here. Instead of sampling directly, you sample from a simple base distribution (e.g., standard normal) and then transform this sample using the mean and variance produced by the neural network, making the sampling process differentiable.

### 2. Evaluate the function for each sample
For each sampled {{< math >}} $z_1^{(m)}$ {{< /math >}}, calculate the value of the function:

{{< math >}} $$h(z_1^{(m)}) = \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1^{(m)}) \| p(\pi_{1g}))$$ {{< /math >}}

### 3. Average the results
The expectation is then approximated by the average of these function evaluations:

{{< math >}} $$E_{q_\phi(z_1|u_1,s_1)} \left[ \sum_g \text{KL}(q_\phi(\pi_{1g}|z_1) \| p(\pi_{1g})) \right] \approx \frac{1}{M} \sum_{m=1}^M h(z_1^{(m)})$$ {{< /math >}}

In VeloVI, typically {{< math >}} $M=1$ {{< /math >}} (a single sample per element of the mini-batch) is used for efficiency in deep learning models, particularly during training. While this introduces more variance into the gradient estimates, it's often sufficient for practical optimization.

---

# Why Simple Averaging Works in Monte Carlo

The expression {{< math >}} $\frac{1}{M} \sum_{m=1}^M h(z_1^{(m)})$ {{< /math >}} does not explicitly show a "weighting by their probabilities according to {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}" in the way one might typically see in a discrete weighted average (e.g., {{< math >}} $\sum x_i P(X=x_i)$ {{< /math >}}).

However, this simple average is the correct Monte Carlo approximation for an expectation when the samples {{< math >}} $z_1^{(m)}$ {{< /math >}} are drawn directly from the distribution {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}.

Let's clarify why this is the case:

## The Definition of Expectation

For a function {{< math >}} $h(z)$ {{< /math >}} and a probability distribution {{< math >}} $q(z)$ {{< /math >}}, the expectation {{< math >}} $E_{q(z)}[h(z)]$ {{< /math >}} is defined as:

- **For continuous z**: {{< math >}} $\int h(z)q(z)dz$ {{< /math >}}
- **For discrete z**: {{< math >}} $\sum_z h(z)q(z)$ {{< /math >}}

## Monte Carlo Approximation (Simple Sampling)

The fundamental idea of Monte Carlo approximation for an expectation is based on the Law of Large Numbers. If you draw M independent and identically distributed (i.i.d.) samples {{< math >}} $z^{(1)}, z^{(2)}, \ldots, z^{(M)}$ {{< /math >}} from the distribution {{< math >}} $q(z)$ {{< /math >}}, then the sample mean of the function values {{< math >}} $h(z^{(m)})$ {{< /math >}} will converge to the true expectation as {{< math >}} $M \to \infty$ {{< /math >}}:

{{< math >}} $$E_{q(z)}[h(z)] \approx \frac{1}{M} \sum_{m=1}^M h(z^{(m)})$$ {{< /math >}}

## Why this works (implicit weighting)

The "weighting" is implicitly handled by the sampling process itself. When you draw samples {{< math >}} $z^{(m)}$ {{< /math >}} directly from {{< math >}} $q_\phi(z_1|u_1,s_1)$ {{< /math >}}:

- Values of {{< math >}} $z_1$ {{< /math >}} that have a higher probability (higher density) under {{< math >}} $q_\phi$ {{< /math >}} will naturally be sampled more frequently.
- Values of {{< math >}} $z_1$ {{< /math >}} that have a lower probability (lower density) under {{< math >}} $q_\phi$ {{< /math >}} will be sampled less frequently.

So, when you take the simple average of {{