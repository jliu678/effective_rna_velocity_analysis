---
title: "Velovi Model N Inference"
date: 2025-07-10
draft: True
---

# VeloVI Model Description

Let {{< math >}} $\alpha_{gk}$ {{< /math >}} be the gene-state-specific reaction rate of transcription.
Let {{< math >}} $\beta_g$ {{< /math >}} be the gene-specific splicing rate constant and let {{< math >}} $\gamma_g$ {{< /math >}} be the
gene-specific degradation rate constant. Each gene has a switching
time {{< math >}} $t^s_g$ {{< /math >}} when the system switches from induction phase to repression
phase.

Given the solution to the ordinary differential equations{{< math >}} $^{11}$ {{< /math >}}, the
unspliced transcript abundance at time {{< math >}} $t_{ng}$ {{< /math >}} for cell {{< math >}} $n$ {{< /math >}} and gene {{< math >}} $g$ {{< /math >}} is
defined as

{{< math >}} $$\bar{u}^{(g)}(t_{ng}, k) := u^0_{gk} e^{-\beta_g(t_{ng} - t^0_{gk})} + \frac{\alpha_{gk}}{\beta_g}(1 - e^{-\beta_g(t_{ng} - t^0_{gk})}) \tag{1}$$ {{< /math >}}

where {{< math >}} $t^0_{gk}$ {{< /math >}} is the initial time of the system in state {{< math >}} $k$ {{< /math >}}. The spliced transcript
abundance is defined as

{{< math >}} $$\bar{s}^{(g)}(t_{ng}, k) := s^0_{gk} e^{-\gamma_g \tau} + \frac{\alpha_{gk}}{\gamma_g} (1 - e^{-\gamma_g(t_{ng} - t^0_{gk})}) + \frac{\alpha_{gk} - \beta_g u^0_{gk}}{\gamma_g - \beta_g} (e^{-\gamma_g(t_{ng} - t^0_{gk})} - e^{-\beta_g(t_{ng} - t^0_{gk})}) \tag{2}$$ {{< /math >}}

## Induction state

For the induction state, {{< math >}} $k = 1$ {{< /math >}}, we have {{< math >}} $u^0_{g1} = 0$ {{< /math >}}, {{< math >}} $s^0_{g1} = 0$ {{< /math >}},
{{< math >}} $\alpha_{g1} > 0$ {{< /math >}} and {{< math >}} $t^0_{g1} = 0$ {{< /math >}}. Thus, the unspliced transcript abundance can then
be expressed as

{{< math >}} $$\bar{u}^{(g)}(t_{ng}, k=1) := \frac{\alpha_{g1}}{\beta_g}(1 - e^{-\beta_g t_{ng}}) \tag{3}$$ {{< /math >}}

Likewise, the spliced transcript abundance can be simplified to

{{< math >}} $$\bar{s}^{(g)}(t_{ng}, k = 1) := \frac{\alpha_{g1}}{\gamma_g}(1 - e^{-\gamma_g t_{ng}}) + \frac{\alpha_{g1}}{\gamma_g - \beta_g} (e^{-\gamma_g t_{ng}} - e^{-\beta_g t_{ng}}) \tag{4}$$ {{< /math >}}

## Induction steady state

For the induction steady state, {{< math >}} $k = 2$ {{< /math >}}, the
unspliced and spliced transcript abundances are defined as limits of
the system:

{{< math >}} $$\bar{u}^{(g)}(t_{ng}, k=1) := \lim_{t_{ng} \to \infty} \bar{u}^{(g)}(t_{ng}, k=1) = \frac{\alpha_{g1}}{\beta_g} \tag{5}$$ {{< /math >}}

{{< math >}} $$\bar{s}^{(g)}(t_{ng}, k = 2) := \lim_{t_{ng} \to \infty} \bar{s}^{(g)}(t^s_g, k = 1) = \frac{\alpha_{g1}}{\gamma_g} \tag{6}$$ {{< /math >}}

## Repression state

For the repression state, {{< math >}} $k = 3$ {{< /math >}}, we have {{< math >}} $\alpha_{g3} = 0$ {{< /math >}} and
{{< math >}} $t^0_{g3} = t^s_g$ {{< /math >}}. Thus, the number of unspliced transcripts can then be
expressed as

{{< math >}} $$\bar{u}^{(g)}(t_{ng}, k=3) := u^0_{g3} e^{-\beta_g(t_{ng} - t^0_{g3})} \tag{7}$$ {{< /math >}}

Likewise, the number of spliced transcripts can be simplified to

{{< math >}} $$\bar{s}^{(g)}(t_{ng}, k = 3) := s^0_{g3} e^{-\gamma_g(t_{ng} - t^0_{g3})} - \frac{\beta_g u^0_{g3}}{\gamma_g - \beta_g} (e^{-\gamma_g \tau} - e^{-\beta_g(t_{ng} - t^0_{g3})}) \tag{8}$$ {{< /math >}}

The initial conditions, {{< math >}} $u^0_{g3}$ {{< /math >}} and {{< math >}} $s^0_{g3}$ {{< /math >}} are defined by the induction
model at the switching time {{< math >}} $t^s_g$ {{< /math >}}, such that

{{< math >}} $$u^0_{g3} = \bar{u}^{(g)}(t^s_g, k = 2) \tag{9}$$ {{< /math >}}

{{< math >}} $$s^0_{g3} = \bar{s}^{(g)}(t^s_g, k = 2) \tag{10}$$ {{< /math >}}

## Repression steady state

For the repression steady state, the limit
upon which {{< math >}} $t_{ng} \to \infty$ {{< /math >}}, there is no expression, so we have

{{< math >}} $$\bar{u}^{(g)}(t_{ng}, k=4) := 0 \tag{11}$$ {{< /math >}}

{{< math >}} $$\bar{s}^{(g)}(t_{ng}, k = 4) := 0 \tag{12}$$ {{< /math >}}

## Model assumptions

As in ref. 11, this model assumes that for one gene,
at the initial time of the system, cells are first in induction phase in which
both spliced and unspliced expression increases. Then cells potentially
reach a steady state of this induction state. Next at some future time {{< math >}} $t^s_g$ {{< /math >}}
the system switches to repression state. Finally, the repression reaches
a steady state in which there is no expression. Further assumptions are
necessary to identify the dynamical model parameters{{< math >}} $^{44}$ {{< /math >}}; thus, we
assume that each gene is on the same time scale (precisely each gene
has a maximum time of {{< math >}} $t = 20$ {{< /math >}} as shown previously{{< math >}} $^{11}$ {{< /math >}}).

# veloVI generative process

We posit a generative process that takes into account the underlying dynamics of the system. Compared to Bergen et al.{{< math >}} $^{11}$ {{< /math >}}, the model
here does not treat each gene independently; instead, the latent time
and states for each (cell and gene) pair are tied together via a local
low-dimensional latent variable.

For each cell we draw a low-dimensional ({{< math >}} $d = 10$ {{< /math >}} dimensions
throughout this manuscript) latent variable

{{< math >}} $$z_n \sim \text{Normal}(0, I_d) \tag{13}$$ {{< /math >}}

that summarizes the latent state of each cell. Next, for each gene {{< math >}} $g$ {{< /math >}} in
cell {{< math >}} $n$ {{< /math >}} we draw the distribution over the state assignments as well as the
state assignment itself

{{< math >}} $$\pi_{ng} \sim \text{Dirichlet}(0.25, 0.25, 0.25, 0.25) \tag{14}$$ {{< /math >}}

{{< math >}} $$k_{ng} \sim \text{Categorical}(\pi_{ng}) \tag{15}$$ {{< /math >}}

Here {{< math >}} $\pi_{ng}$ {{< /math >}} is sampled from a Dirichlet distribution, which has the
support of the probability simplex. In other words, the Dirichlet provides a distribution over discrete probability distributions. If {{< math >}} $k_{ng} = 1$ {{< /math >}}
(induction), then the time is a function of {{< math >}} $z_n$ {{< /math >}},

{{< math >}} $$\rho^{(1)}_{ng} = [h_{\text{ind}}(z_n)]_g \tag{16}$$ {{< /math >}}

{{< math >}} $$(t^1)_{ng} = \rho^{(1)}_{ng} t^s_g \tag{17}$$ {{< /math >}}

where {{< math >}} $h_{\text{ind}} : \mathbb{R}^d \to (0, 1)^G$ {{< /math >}} is parameterized as a fully connected neural
network. Notably, this parameterization results in an induction-specific
time that is constrained to be less than the switching time.

Else, if {{< math >}} $k_{ng} = 3$ {{< /math >}} (repression),

{{< math >}} $$\rho^{(3)}_{ng} = [h_{\text{rep}}(z_n)]_g \tag{18}$$ {{< /math >}}

{{< math >}} $$(t^3)_{ng} = (t_{\max} - t^s_g)\rho^{(3)}_{ng} + t^s_g \tag{19}$$ {{< /math >}}

where {{< math >}} $t_{\max} := 20$ {{< /math >}} is used to fix the time scale across genes and identify
the rate parameters of the model. Similarly to the previously defined
function, {{< math >}} $h_{\text{rep}} : \mathbb{R}^d \to (0, 1)^G$ {{< /math >}} and is also a neural network.

We also consider two potential steady states. If {{< math >}} $k_{ng} = 2$ {{< /math >}} (induction
steady state) or if {{< math >}} $k_{ng} = 4$ {{< /math >}} (repression steady state), we consider the
limit as time approaches {{< math >}} $\infty$ {{< /math >}}, which is described in the previous section.

Finally, the observed data are sampled from normal distributions as

{{< math >}} $$u_{ng} \sim \text{Normal} \left(\bar{u}^{(g)}\left(t_{ng}^{(k_{ng})}, k_{ng}\right), (c_k\sigma_{gu})^2\right) \tag{20} $${{< /math >}}

{{< math >}} $$s_{ng} \sim \text{Normal} \left(\bar{s}^{(g)}\left(t_{ng}^{(k_{ng})}, k_{ng}\right), (c_k\sigma_{gs})^2\right) \tag{21} $$ {{< /math >}}

For veloVI, we consider the observed data {{< math >}} $\{(s_n, u_n)\}_{n=1}^N$ {{< /math >}} to be the nearest-neighbor smoothed expression data that is also used as input to scVelo as well as velocyto. In addition, we assume the data have been preprocessed such that for each gene, the smoothed spliced and unspliced abundances are independently min-max scaled into [0, 1]. By using the normal distribution, we assume that the smoothed expression (which represents an average of random variables) has a sampling distribution centered on some mean value and that this sampling distribution is approximately normal; however, the flexibility of this modeling framework will enable extensions that consider the discrete nature of unique molecular identifiers used in standard scRNA-seq assays.

We include a state-dependent scaling factor on the variance. For all experiments in this manuscript, we used {{< math >}} $c_k = 1$ {{< /math >}} except for the repression steady state in which {{< math >}} $c_4 = 0.1$ {{< /math >}}. This hyperparameter choice forces the variance of abundance in the repression steady state to be less than that of other transcriptional states, which reflects the notion that the repression steady state corresponds to zero transcriptional activity. Despite the assumption of zero transcriptional activity, the normal distribution here captures noise that arises during the experimental process (ambient transcripts) as well as during preprocessing (for example, KNN smoothing). Finally, in the following, let {{< math >}} $\theta$ {{< /math >}} be the set of parameters of the generative process ({{< math >}} $\alpha, \beta, \gamma, t_s$ {{< /math >}} and neural network parameters).

# Multivariate Normal Distribution Described in Equation (13)

The expression:

{{< math >}} 
$$z_n \sim \text{Normal}(0, I_d) \tag{13}$$
{{< /math >}}

means that the random variable {{< math >}} $z_n$ {{< /math >}} is drawn from a multivariate normal (Gaussian) distribution with:

- **Mean vector** {{< math >}} $0$ {{< /math >}}: a vector of all zeros, of dimension {{< math >}} $d$ {{< /math >}}
- **Covariance matrix** {{< math >}} $I_d$ {{< /math >}}: the {{< math >}} $d \times d$ {{< /math >}} identity matrix

## 🔍 What this means:

{{< math >}} $z_n \in \mathbb{R}^d$ {{< /math >}}, so it's a {{< math >}} $d$ {{< /math >}}-dimensional vector.

Each component of {{< math >}} $z_n$ {{< /math >}} is independent and follows a standard normal distribution:

{{< math >}} 
$$z_n^{(i)} \sim N(0,1) \quad \text{for } i=1,\ldots,d \tag{13.1}$$
{{< /math >}}

The identity matrix {{< math >}} $I_d$ {{< /math >}} as covariance means there's no correlation between dimensions.

## ✅ Common context:

This type of distribution is often used in Variational Autoencoders (VAEs) or other latent variable models, where:

- {{< math >}} $z_n$ {{< /math >}} is a latent (hidden) variable representing compressed features.
- We assume a standard normal prior for simplicity and tractability.

## 📌 Simple Numeric Example

Let's take:

**Dimension** {{< math >}} $d = 3$ {{< /math >}}

So, {{< math >}} $z \sim N(0, I_3)$ {{< /math >}}, where

{{< math >}} 
$$I_3 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \tag{13.2}$$
{{< /math >}}

This means {{< math >}} $z = [z_1, z_2, z_3]$ {{< /math >}}, where each {{< math >}} $z_i \sim N(0,1)$ {{< /math >}}, independently.

## 🎲 Sample values (drawn randomly):

Suppose we draw one sample and get:

{{< math >}} 
$$z = \begin{bmatrix} 0.15 \\ -1.23 \\ 0.84 \end{bmatrix} \tag{13.3}$$
{{< /math >}}

That means:
- {{< math >}} $z_1 = 0.15$ {{< /math >}}
- {{< math >}} $z_2 = -1.23$ {{< /math >}}
- {{< math >}} $z_3 = 0.84$ {{< /math >}}

and each of these values came from a standard normal distribution.

## 🔧 In Python (if you're curious):

```python
import numpy as np

z = np.random.normal(0, 1, size=3)
print(z)
```

This will give a different 3D sample from {{< math >}} $N(0, I_3)$ {{< /math >}} each time you run it.

# veloVI inference procedure

We seek the following: (1) point estimates of the transcription rate, degradation and splicing rate constants and the switching time point; (2) point estimates of the parameters of the neural networks; and (3), a posterior distribution over the latent variables, which in this case includes {{< math >}} $z$ {{< /math >}} and {{< math >}} $\pi$ {{< /math >}}. Noting that the model evidence {{< math >}} $p_\theta(u, s)$ {{< /math >}} cannot be computed in closed form, we use variational inference to approximate the posterior distribution as well as accomplish the other tasks. Following inference, velocity can be calculated as a functional of the variational posterior distribution.

### Variational posterior

We posit the following factorization on the approximate posterior distribution

{{< math >}} $$q_\phi(z, \pi | u, s) := \prod_{n}^N q_\phi(z_n | u_n, s_n) \prod_{g}^G q_\phi(\pi_{ng} | z_n), \tag{22} $$ {{< /math >}}

in which dependencies are specified using neural networks with parameter set {{< math >}} $\phi$ {{< /math >}}. Here {{< math >}} $z$ {{< /math >}} factorizes over all {{< math >}} $n$ {{< /math >}} cells and {{< math >}} $\pi_{ng}$ {{< /math >}} over all {{< math >}} $n$ {{< /math >}} cells and {{< math >}} $g$ {{< /math >}} genes.

For the likelihoods, we integrate over the choice of transcriptional state {{< math >}} $k_{ng}$ {{< /math >}}, such that the likelihoods for unspliced and spliced transcript abundances,

{{< math >}} $$p_\theta(u_{ng} | z_n, \pi_n) = \sum_{k_{ng} \in \{1,2,3,4\}} \pi_{ng}^{k_{ng}} \text{Normal} \left(\bar{u}^{(g)}\left(t_{ng}^{(k_{ng})}, k_{ng}\right), (c_k\sigma_{gu})^2\right) \tag{23} $$ {{< /math >}}

{{< math >}} $$p_\theta(s_{ng} | z_n, \pi_n) = \sum_{k_{ng} \in \{1,2,3,4\}} \pi_{ng}^{k_{ng}} \text{Normal} \left(\bar{s}^{(g)}\left(t_{ng}^{(k_{ng})}, k_{ng}\right), (c_k\sigma_{gs})^2\right)\tag{24} $$ {{< /math >}}

are mixtures of normal distributions.

### Objective

The objective that is minimized during inference is composed of two terms

{{< math >}} $$\mathcal{L}_{\text{velo}}(\theta, \phi; u, s) = \mathcal{L}_{\text{elbo}}(\theta, \phi; u, s) + \lambda \mathcal{L}_{\text{switch}}(\theta; u, s), \tag{25} $$ {{< /math >}}

where {{< math >}} $\mathcal{L}_{\text{elbo}}$ {{< /math >}} is the negative evidence lower bound of {{< math >}} $\log p_\theta(u, s)$ {{< /math >}} and {{< math >}} $\mathcal{L}_{\text{switch}}$ {{< /math >}} is an additional penalty that regularizes the location of the transcriptional switch in the phase portrait. In more detail,

{{< math >}} $$\mathcal{L}_{\text{elbo}}(\theta, \phi; u, s) = \sum_n - \mathbb{E}_{q_\phi(z_n, \pi_n | u_n, s_n)} [\log p_\theta(u_n, s_n | z_n, \pi_n)]$$ {{< /math >}}
{{< math >}} $$+ \text{KL} (q_\phi(z_n | u_n, s_n) \| p(z))$$ {{< /math >}}
{{< math >}} $$+ \mathbb{E}_{q_\phi(z_n | u_n, s_n)} \left[\sum_g \text{KL} (q_\phi(\pi_{ng} | z_n) \| p(\pi_{ng}))\right], \tag{26}$$ {{< /math >}}

{{< math >}} $$\phantom{\mathcal{L}_{\text{elbo}}(\theta, \phi; u, s) = \sum_n - \mathbb{E}_{q_\phi(z_n, \pi_n | u_n, s_n)} [\log p_\theta(u_n, s_n | z_n, \pi_n)]}$$ {{< /math >}}

which can be estimated using minibatches of data. In particular, we use randomly sampled minibatches of 256 cells for inference. For the penalty term {{< math >}} $\mathcal{L}_{\text{switch}}$ {{< /math >}}, we start by only considering cells that are above the 99th percentile of unspliced abundance for each gene. Using these cells we compute the median unspliced and spliced abundance for each gene separately. Let {{< math >}} $u^*$ {{< /math >}} and {{< math >}} $s^*$ {{< /math >}} be the outcome of this procedure, then

{{< math >}} $$\mathcal{L}_{\text{switch}}(\theta; u, s) = \sum_g (u_{0g3} - u_g^*)^2 + (s_{0g3} - s_g^*)^2, \tag{27}$$ {{< /math >}}

where {{< math >}} $u_{0g3}$ {{< /math >}} and {{< math >}} $s_{0g3}$ {{< /math >}} were defined as the initial conditions of the repression phase at the switch time {{< math >}} $t_{sg}$ {{< /math >}}.

### Initialization

We initialize {{< math >}} $\alpha_{g1}$ {{< /math >}} to be equal to the median unspliced abundance for the cells above the 99th percentile for each gene. The other global parameters, including the splicing, degradation and switch time are initialized to a constant value shared by all genes. All neural network initialization uses the default implementation in PyTorch.

### Optimization

To optimize {{< math >}} $\mathcal{L}_{\text{velo}}$ {{< /math >}} we use stochastic gradients along with the Adam optimizer with weight decay as implemented in PyTorch. For all experiments we use {{< math >}} $\lambda = 0.2$ {{< /math >}} for scaling the regularization term in the loss. As a result of minibatching, veloVI's memory usage is constant throughout training. Unless otherwise specified, all neural networks are fully connected feedforward networks that use standard activation functions such as ReLU for hidden layers and softplus or exponential for parameterizing non-negative distributional parameters.

### Architecture 

An overview of the veloVI architecture is shown in Supplementary Fig. 9.


