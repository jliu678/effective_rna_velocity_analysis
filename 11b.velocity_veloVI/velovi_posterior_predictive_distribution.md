---
title: "Velovi Posterior Predictive Distribution"
date: 2025-07-10
draft: True
---

# VeloVI Posterior Predictive Distribution: Theory and Numerical Example

## Introduction

Equation (28) describes the concept of a posterior predictive distribution for the unspliced abundance ({{< math >}}$u_n^*${{< /math >}}) of a cell {{< math >}}$n${{< /math >}} in the VeloVI model. It's a powerful idea that allows us to understand what new, unobserved data would look like, given our current model and the observed data.

## Understanding the Components

### Unobserved Random Variables

{{< math >}}$u_n^*${{< /math >}} and {{< math >}}$s_n^*${{< /math >}} are unobserved random variables representing posterior predictive values of unspliced and spliced abundances for cell {{< math >}}$n${{< /math >}}.

This means we're not talking about the actual observed {{< math >}}$u_n, s_n${{< /math >}} that we fed into the model. Instead, we're thinking about hypothetical, new (unseen) unspliced and spliced RNA counts for this same cell {{< math >}}$n${{< /math >}}, assuming its underlying biological state (represented by {{< math >}}$z_n${{< /math >}} and {{< math >}}$\pi_n${{< /math >}}) is consistent with what we observed.

"Posterior predictive" means we're predicting after having observed data ({{< math >}}$u_n, s_n${{< /math >}}) and after having inferred the likely latent states ({{< math >}}$z_n, \pi_n${{< /math >}}) for that data.

### The Posterior Predictive Distribution

The main equation we want to calculate is:

{{< math >}}
$$
p(u_n^* \mid u_n, s_n) = \int p_\theta(u_n^* \mid z_n, \pi_n) q_\phi(z_n, \pi_n \mid u_n, s_n) d\pi_n dz_n \tag{28}
$$
{{< /math >}}

This is what we want to calculate. It's the probability of observing a new unspliced count {{< math >}}$u_n^*${{< /math >}}, given the original observed counts {{< math >}}$u_n${{< /math >}} and {{< math >}}$s_n${{< /math >}} for cell {{< math >}}$n${{< /math >}}.

It essentially answers: "If I were to measure cell {{< math >}}$n${{< /math >}} again, what kind of unspliced count would I expect to see, given what I already know about it?"

### The Integral (Marginalization)

The integral sign {{< math >}}$\int \ldots d\pi_n dz_n${{< /math >}} means we are marginalizing out (summing over all possible values of) the latent variables {{< math >}}$z_n${{< /math >}} and {{< math >}}$\pi_n${{< /math >}}.

Why? Because {{< math >}}$u_n^*${{< /math >}} doesn't directly depend on {{< math >}}$u_n, s_n${{< /math >}}. Instead, {{< math >}}$u_n^*${{< /math >}} depends on the cell's underlying biological state ({{< math >}}$z_n, \pi_n${{< /math >}}), which we inferred from {{< math >}}$u_n, s_n${{< /math >}}. To get the total probability of {{< math >}}$u_n^*${{< /math >}}, we need to consider all possible values of {{< math >}}$z_n${{< /math >}} and {{< math >}}$\pi_n${{< /math >}} that could have led to the observed {{< math >}}$u_n, s_n${{< /math >}}, weighted by how likely they are.

### The Generative Model / Decoder

{{< math >}}$p_\theta(u_n^* \mid z_n, \pi_n)${{< /math >}} (The Generative Model / Decoder):

This term is part of the decoder or generative model of VeloVI. It describes the probability of generating an unspliced count {{< math >}}$u_n^*${{< /math >}}, given specific values of the latent variable {{< math >}}$z_n${{< /math >}} and the kinetic state probabilities {{< math >}}$\pi_n${{< /math >}}, and using the biophysical parameters {{< math >}}$\theta${{< /math >}} (like {{< math >}}$\alpha, \beta, \gamma${{< /math >}}) that govern the RNA dynamics.

In VeloVI, {{< math >}}$z_n${{< /math >}} helps determine the latent time, and {{< math >}}$\pi_n${{< /math >}} determines the kinetic state (e.g., induction, repression). These then feed into the kinetic equations (which use {{< math >}}$\theta${{< /math >}}) to predict the mean expected abundance {{< math >}}$\bar{u}^{(g)}${{< /math >}}. {{< math >}}$p_\theta(u_n^* \mid z_n, \pi_n)${{< /math >}} would then be a distribution (like a Normal distribution) centered around that predicted mean, with variance {{< math >}}$\sigma_{gu}${{< /math >}}.

### The Approximate Posterior / Encoder

{{< math >}}$q_\phi(z_n, \pi_n \mid u_n, s_n)${{< /math >}} (The Approximate Posterior / Encoder):

This is the approximate posterior distribution provided by the encoder network. It represents VeloVI's learned understanding of the probability of a cell having latent state ({{< math >}}$z_n, \pi_n${{< /math >}}), given its observed unspliced and spliced counts ({{< math >}}$u_n, s_n${{< /math >}}).

This is the distribution from which we draw samples ({{< math >}}$z_n^*, \pi_n^*${{< /math >}}) during training. It's parameterized by the neural network parameters {{< math >}}$\phi${{< /math >}}.

## Putting it All Together: The Intuition

Equation (28) essentially says:

"To predict a new, unobserved unspliced count ({{< math >}}$u_n^*${{< /math >}}) for a cell ({{< math >}}$n${{< /math >}}) that we've already measured ({{< math >}}$u_n, s_n${{< /math >}}):

1. Consider all possible underlying biological states ({{< math >}}$z_n, \pi_n${{< /math >}}) that could explain the observed data ({{< math >}}$u_n, s_n${{< /math >}}). The encoder ({{< math >}}$q_\phi${{< /math >}}) tells us how likely each of these states is.

2. For each of these possible states, use the generative model ({{< math >}}$p_\theta${{< /math >}}) to predict what unspliced count ({{< math >}}$u_n^*${{< /math >}}) would be produced.

3. Sum up these predictions, weighted by how likely each underlying state is."

## Why is this Important?

- **Model Evaluation**: It allows you to check if your model makes reasonable predictions for new data. If the predicted {{< math >}}$u_n^*${{< /math >}} looks very different from the original {{< math >}}$u_n${{< /math >}}, it suggests a problem with your model's fit or assumptions.

- **Missing Data Imputation (conceptual)**: While not explicitly used for imputation in VeloVI's core, posterior predictive distributions are the foundation for imputing missing data or generating synthetic data that resembles the observed data.

- **Uncertainty Quantification**: The posterior predictive distribution gives you not just a single predicted value, but a full distribution, allowing you to quantify the uncertainty in your predictions.

## Monte Carlo Approximation

In practice, this integral is usually intractable (cannot be solved analytically). So, in VAEs, it's often approximated using Monte Carlo sampling:

1. Draw many samples of ({{< math >}}$z_n, \pi_n${{< /math >}}) from the approximate posterior {{< math >}}$q_\phi(z_n, \pi_n \mid u_n, s_n)${{< /math >}}.

2. For each sample, draw a new {{< math >}}$u_n^*${{< /math >}} from {{< math >}}$p_\theta(u_n^* \mid z_n, \pi_n)${{< /math >}}.

3. The collection of these {{< math >}}$u_n^*${{< /math >}} samples forms an empirical approximation of the posterior predictive distribution.

## Numerical Example

Let's calculate (or more accurately, approximate) Equation (28) using numeric examples. We'll use our scenario of Cell 1 and Gene A for the numeric example.

### Our Goal
Generate samples of {{< math >}}$u_{1A}^*${{< /math >}} to approximate the distribution {{< math >}}$p(u_{1A}^* \mid u_{1A}, s_{1A})${{< /math >}}.

### What We Assume We Have
After VeloVI has been trained and its parameters {{< math >}}$\theta${{< /math >}} and {{< math >}}$\phi${{< /math >}} have converged:

#### Observed Data for Cell 1, Gene A:
- {{< math >}}$u_{1A} = 0.8${{< /math >}}
- {{< math >}}$s_{1A} = 0.2${{< /math >}} (These are the inputs to our encoder)

#### Learned Biophysical Parameters ({{< math >}}$\theta${{< /math >}}) for Gene A:
(These are fixed values after training)

- Transcription Rate ({{< math >}}$\alpha_{A1}${{< /math >}}) = 1.2
- Splicing Rate ({{< math >}}$\beta_A${{< /math >}}) = 0.2
- Degradation Rate ({{< math >}}$\gamma_A${{< /math >}}) = 0.1
- Switching Time ({{< math >}}$t_A^s${{< /math >}}) = 8.0
- Unspliced Std. Dev. ({{< math >}}$\sigma_{Au}${{< /math >}}) = 0.15
- Spliced Std. Dev. ({{< math >}}$\sigma_{As}${{< /math >}}) = 0.15

(For simplicity, let's assume all state scaling factors {{< math >}}$c_k = 1${{< /math >}} for this example, though VeloVI uses a specific {{< math >}}$c_4 = 0.1${{< /math >}})

#### Learned Neural Network Parameters ({{< math >}}$\phi${{< /math >}}) for the Encoder and Decoder:
(These are the weights and biases of the networks. We don't see the numbers directly, but we know their effect on the outputs.)

**Encoder's Output for {{< math >}}$z_1${{< /math >}} (given {{< math >}}$u_1, s_1${{< /math >}})**: This defines the approximate posterior {{< math >}}$q_\phi(z_1 \mid u_1, s_1)${{< /math >}}, which is a Gaussian.

- Mean ({{< math >}}$\mu_{z_1}${{< /math >}}) = (0.5, -0.3)
- Log-Variance ({{< math >}}$\log\Sigma_{z_1}${{< /math >}}) = (-0.5, -0.5)

This implies {{< math >}}$\Sigma_{z_1} = \text{diag}(e^{-0.5}, e^{-0.5}) \approx \text{diag}(0.6065, 0.6065)${{< /math >}} (Let's assume {{< math >}}$z_n${{< /math >}} is 2-dimensional for this example)

**Encoder's Output for {{< math >}}$\pi_{1A}${{< /math >}} (given {{< math >}}$z_1${{< /math >}})**: This defines {{< math >}}$q_\phi(\pi_{1A} \mid z_1)${{< /math >}}, which is a Dirichlet.

Example concentration parameters (output for a specific {{< math >}}$z_1${{< /math >}}): {{< math >}}$\text{conc}_{1A} = (5.0, 0.1, 0.2, 0.1)${{< /math >}} (This strongly favors State 1: Induction)

**Decoder's Output for {{< math >}}$\rho_{1A}^{(k)}${{< /math >}} (given {{< math >}}$z_1${{< /math >}})**: These networks take {{< math >}}$z_1${{< /math >}} and output scaled latent times.

- Example output for a specific {{< math >}}$z_1${{< /math >}}: {{< math >}}$\rho_{1A}^{(1)} = 0.5${{< /math >}} (for state 1)
- Example output for a specific {{< math >}}$z_1${{< /math >}}: {{< math >}}$\rho_{1A}^{(3)} = 0.8${{< /math >}} (for state 3)

### The Monte Carlo Approximation Process (Step-by-Step Numerical Example)

We will repeat the following steps many times (say, {{< math >}}$M = 1000${{< /math >}} times) to get enough samples of {{< math >}}$u_{1A}^*${{< /math >}}. For this explanation, we'll trace just one single sample to illustrate the process.

#### Iteration 1 of M (Sample #1):

**Step 1**: Sample latent variables ({{< math >}}$z_1^*, \pi_{1A}^*${{< /math >}}, and determine {{< math >}}$k_{1A}${{< /math >}}) from the approximate posterior {{< math >}}$q_\phi(z_1, \pi_{1A} \mid u_{1A}, s_{1A})${{< /math >}}.

**1a. Sample {{< math >}}$z_1^*${{< /math >}} from {{< math >}}$q_\phi(z_1 \mid u_{1A}, s_{1A})${{< /math >}} using the reparameterization trick:**

Recall {{< math >}}$q_\phi(z_1 \mid u_{1A}, s_{1A})${{< /math >}} is {{< math >}}$\mathcal{N}(\mu_{z_1}, \Sigma_{z_1})${{< /math >}}, where {{< math >}}$\mu_{z_1} = (0.5, -0.3)${{< /math >}} and {{< math >}}$\Sigma_{z_1} = \text{diag}(0.6065, 0.6065)${{< /math >}}.

- Draw a random standard normal noise vector {{< math >}}$\epsilon_z${{< /math >}}: Let {{< math >}}$\epsilon_z = (0.2, -0.1)${{< /math >}}.
- Calculate {{< math >}}$z_1^* = \mu_{z_1} + \sqrt{\Sigma_{z_1}} \odot \epsilon_z${{< /math >}} (where {{< math >}}$\odot${{< /math >}} is element-wise multiplication):

{{< math >}}
$$
z_1^* = (0.5, -0.3) + (\sqrt{0.6065}, \sqrt{0.6065}) \odot (0.2, -0.1) \tag{29}
$$
{{< /math >}}

{{< math >}}
$$
z_1^* = (0.5, -0.3) + (0.7788, 0.7788) \odot (0.2, -0.1) \tag{30}
$$
{{< /math >}}

{{< math >}}
$$
z_1^* = (0.5, -0.3) + (0.1558, -0.0779) \tag{31}
$$
{{< /math >}}

{{< math >}}
$$
z_1^* = (0.6558, -0.3779) \tag{32}
$$
{{< /math >}}

**1b. Given {{< math >}}$z_1^*${{< /math >}}, sample {{< math >}}$\pi_{1A}^*${{< /math >}} from {{< math >}}$q_\phi(\pi_{1A} \mid z_1^*)${{< /math >}}:**

- The encoder (a specific part of it) takes this sampled {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}} as input and outputs the concentration parameters for {{< math >}}$\pi_{1A}${{< /math >}}.
- Let's assume for {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}}, the encoder outputs: {{< math >}}$\text{conc}_{1A} = (5.0, 0.1, 0.2, 0.1)${{< /math >}}.
- Sample {{< math >}}$\pi_{1A}^*${{< /math >}} from {{< math >}}$\text{Dirichlet}(5.0, 0.1, 0.2, 0.1)${{< /math >}}. This is a stochastic step.
- Let's say we get the sample: {{< math >}}$\pi_{1A}^* = (0.95, 0.02, 0.02, 0.01)${{< /math >}}.

**1c. Determine the kinetic state {{< math >}}$k_{1A}${{< /math >}} from {{< math >}}$\pi_{1A}^*${{< /math >}}:**

- We choose the state with the highest probability. Here, 0.95 is the highest, corresponding to the first state.
- So, {{< math >}}$k_{1A} = 1${{< /math >}} (Induction).

**Step 2**: Generate a new unspliced abundance {{< math >}}$u_{1A}^*${{< /math >}} from the generative model {{< math >}}$p_\theta(u_{1A}^* \mid z_1^*, \pi_{1A}^*)${{< /math >}}.

**2a. Decode {{< math >}}$z_1^*${{< /math >}} to latent time ({{< math >}}$t_{1A}^{(k)}${{< /math >}}):**

- The decoder network (specifically, {{< math >}}$h_{\text{ind}}${{< /math >}} for Induction states) takes {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}} and outputs a scaled time {{< math >}}$\rho_{1A}^{(1)}${{< /math >}}.
- Let's assume the decoder outputs {{< math >}}$\rho_{1A}^{(1)} = 0.5${{< /math >}}.
- Calculate the actual latent time: {{< math >}}$t_{1A}^{(1)} = \rho_{1A}^{(1)} \times t_A^s = 0.5 \times 8.0 = 4.0${{< /math >}}.

**2b. Calculate the predicted mean unspliced abundance {{< math >}}$\bar{u}^{(A)}${{< /math >}} using the kinetic equations and {{< math >}}$\theta${{< /math >}}:**

Since {{< math >}}$k_{1A} = 1${{< /math >}} (Induction), we use the kinetic equation for the induction phase.

{{< math >}}
$$
\bar{u}^{(A)}(t=4.0, k=1) = \frac{\alpha_{A1}}{\beta_A}(1 - e^{-\beta_A t}) \tag{33}
$$
{{< /math >}}

{{< math >}}
$$
\bar{u}^{(A)}(4.0, 1) = \frac{1.2}{0.2}(1 - e^{-0.2 \times 4.0}) \tag{34}
$$
{{< /math >}}

{{< math >}}
$$
\bar{u}^{(A)}(4.0, 1) = 6(1 - e^{-0.8}) \tag{35}
$$
{{< /math >}}

{{< math >}}
$$
\bar{u}^{(A)}(4.0, 1) = 6(1 - 0.4493) \tag{36}
$$
{{< /math >}}

{{< math >}}
$$
\bar{u}^{(A)}(4.0, 1) = 6 \times 0.5507 = 3.3042 \tag{37}
$$
{{< /math >}}

**2c. Sample {{< math >}}$u_{1A}^*${{< /math >}} from {{< math >}}$p_\theta(u_{1A}^* \mid z_1^*, \pi_{1A}^*)${{< /math >}}:**

- This distribution is a Normal distribution: {{< math >}}$\mathcal{N}(\text{mean} = \bar{u}^{(A)}, \text{variance} = \sigma_{Au}^2)${{< /math >}}.
- So, {{< math >}}$u_{1A}^* \sim \mathcal{N}(3.3042, 0.15^2)${{< /math >}}.
- Draw one sample from this Normal distribution (this is again a stochastic step).
- Let's say we sample: {{< math >}}$u_{1A}^* = 3.25${{< /math >}}. This is our first sample for the posterior predictive distribution.

### Repeat for M Samples (Conceptual)

You would repeat these two main steps (Step 1 and Step 2) many times (e.g., {{< math >}}$M = 1000${{< /math >}}). Each repetition would involve new random draws for {{< math >}}$\epsilon_z${{< /math >}}, new Dirichlet samples for {{< math >}}$\pi_{1A}^*${{< /math >}}, and new Normal samples for {{< math >}}$u_{1A}^*${{< /math >}}.

For instance, if we ran a second iteration, we might get:

**Sample #2**: {{< math >}}$\epsilon_z = (-0.5, 0.3) \rightarrow z_1^* = (0.0106, -0.0663) \rightarrow \text{conc}_{1A} = (4.0, 0.3, 0.1, 0.1) \rightarrow \pi_{1A}^* = (0.90, 0.07, 0.02, 0.01) \rightarrow k_{1A} = 1 \rightarrow \rho_{1A}^{(1)} = 0.4 \rightarrow t_{1A}^{(1)} = 3.2 \rightarrow \bar{u}^{(A)} = 2.75 \rightarrow u_{1A}^* = 2.80${{< /math >}}

## Result

After {{< math >}}$M${{< /math >}} iterations, you would have a collection of {{< math >}}$M${{< /math >}} sampled {{< math >}}$u_{1A}^*${{< /math >}} values: {{< math >}}$\{u_{1A}^{*(1)}, u_{1A}^{*(2)}, \ldots, u_{1A}^{*(M)}\}${{< /math >}}.

This collection of samples empirically approximates the posterior predictive distribution {{< math >}}$p(u_{1A}^* \mid u_{1A}, s_{1A})${{< /math >}}. You can then:

1. Plot a histogram of these samples to visualize the distribution.
2. Calculate the mean of these samples (e.g., {{< math >}}$(3.25 + 2.80 + \ldots)/M${{< /math >}}) as your model's best prediction for a new unspliced measurement.
3. Calculate the standard deviation or percentiles to understand the uncertainty in the prediction.

This entire process allows VeloVI to "imagine" what the unspliced counts should be, given its understanding of the cell's underlying biological state and the learned kinetic parameters, thereby providing a powerful way to interpret its learned representations.