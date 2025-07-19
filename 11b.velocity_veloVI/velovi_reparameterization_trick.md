---
title: "Velovi Reparameterization Trick"
date: 2025-07-10
draft: True
---

# The Reparameterization Trick

The "reparameterization trick" is fundamental to how Variational Autoencoders (VAEs) and models like VeloVI (which uses VAE principles) can be trained using gradient-based optimization methods like backpropagation.

## The Problem

**Sampling is not differentiable**: If you have a random variable {{< math >}} $z$ {{< /math >}} drawn directly from a distribution {{< math >}} $q(z|\mu,\Sigma)$ {{< /math >}}, you cannot calculate {{< math >}} $\frac{\partial z}{\partial \mu}$ {{< /math >}} or {{< math >}} $\frac{\partial z}{\partial \Sigma}$ {{< /math >}} in a meaningful way for backpropagation. The sampling operation itself is stochastic and doesn't have a well-defined gradient in the traditional sense. This would break the computational graph.

**Need for gradient flow**: However, the loss function (ELBO) depends on the sampled {{< math >}} $z^*$ {{< /math >}}, and we need to update {{< math >}} $\mu$ {{< /math >}} and {{< math >}} $\Sigma$ {{< /math >}} (which are outputs of the encoder neural network) to minimize this loss.

The reparameterization trick solves this by separating the stochasticity from the parameters.

## How the Reparameterization Trick Works for a Gaussian Distribution

For a Gaussian distribution {{< math >}} $z \sim \mathcal{N}(\mu, \Sigma)$ {{< /math >}}, we can reparameterize the sampling process as follows:

### 1. Sample from a fixed, simple base distribution

First, we sample a "noise" variable {{< math >}} $\epsilon$ {{< /math >}} from a simple, fixed distribution, usually a standard normal distribution:

{{< math >}} $$\epsilon \sim \mathcal{N}(0, I_d)$$ {{< /math >}}

where {{< math >}} $I_d$ {{< /math >}} is the identity matrix (meaning {{< math >}} $\epsilon$ {{< /math >}} is a vector of independent standard normal random variables). This step is random, but crucially, it does not depend on any model parameters.

### 2. Deterministic Transformation

Then, we deterministically transform this noise variable {{< math >}} $\epsilon$ {{< /math >}} using the parameters {{< math >}} $\mu$ {{< /math >}} and {{< math >}} $\Sigma$ {{< /math >}} (or often, the standard deviation {{< math >}} $\sigma = \sqrt{\Sigma}$ {{< /math >}}, or {{< math >}} $\exp(0.5 \log \Sigma)$ {{< /math >}} for element-wise operation) that are output by the encoder:

{{< math >}} $$z^* = \mu + \sigma \odot \epsilon$$ {{< /math >}}

(where {{< math >}} $\odot$ {{< /math >}} denotes element-wise multiplication if {{< math >}} $z$, $\mu$, $\sigma$, and $\epsilon$ {{< /math >}} are vectors).

## Why this allows gradient flow

Now, let's trace the flow for backpropagation:

**Loss dependence on {{< math >}} $z^*$ {{< /math >}}**: The downstream parts of the model (decoder, kinetic equations, likelihood calculation) compute a loss based on {{< math >}} $z^*$ {{< /math >}}. So, we can compute {{< math >}} $\frac{\partial L}{\partial z^*}$ {{< /math >}} using standard backpropagation rules.

**Chain Rule to {{< math >}} $\mu$ {{< /math >}} and {{< math >}} $\sigma$ {{< /math >}}**: Since {{< math >}} $z^*$ {{< /math >}} is now a deterministic function of {{< math >}} $\mu$, $\sigma$, and $\epsilon$ {{< /math >}}:

For {{< math >}} $\mu$ {{< /math >}}: We can apply the chain rule:
{{< math >}} $$\frac{\partial L}{\partial \mu} = \frac{\partial L}{\partial z^*} \times \frac{\partial z^*}{\partial \mu} = \frac{\partial L}{\partial z^*} \times 1$$ {{< /math >}}

(Since {{< math >}} $z^* = \mu + \sigma \odot \epsilon$ {{< /math >}}, the derivative of {{< math >}} $z^*$ {{< /math >}} with respect to {{< math >}} $\mu$ {{< /math >}} is simply 1).

For {{< math >}} $\sigma$ {{< /math >}}: Similarly:
{{< math >}} $$\frac{\partial L}{\partial \sigma} = \frac{\partial L}{\partial z^*} \times \frac{\partial z^*}{\partial \sigma} = \frac{\partial L}{\partial z^*} \times \epsilon$$ {{< /math >}}

(Since {{< math >}} $z^* = \mu + \sigma \odot \epsilon$ {{< /math >}}, the derivative of {{< math >}} $z^*$ {{< /math >}} with respect to {{< math >}} $\sigma$ {{< /math >}} is {{< math >}} $\epsilon$ {{< /math >}}). 

**Note**: If the encoder outputs {{< math >}} $\log \Sigma$ {{< /math >}} (log-variance) instead of {{< math >}} $\sigma$ {{< /math >}}, the chain rule would be applied one more step:
{{< math >}} $$\frac{\partial L}{\partial \log \Sigma} = \frac{\partial L}{\partial \sigma} \times \frac{\partial \sigma}{\partial \log \Sigma} = \frac{\partial L}{\partial \sigma} \times \frac{1}{2} \exp(0.5 \log \Sigma) = \frac{\partial L}{\partial \sigma} \times \frac{1}{2} \sigma$$ {{< /math >}}

**No Gradient to {{< math >}} $\epsilon$ {{< /math >}}**: Notice that {{< math >}} $\epsilon$ {{< /math >}} is a random variable, but we do not need to compute gradients with respect to {{< math >}} $\epsilon$ {{< /math >}} itself. {{< math >}} $\epsilon$ {{< /math >}} is fixed once sampled, and it contributes to the gradient flow in a deterministic way (as a multiplier for {{< math >}} $\frac{\partial L}{\partial \sigma}$ {{< /math >}}).

## Numerical Example

Let's use a very simplified numeric example to illustrate the gradient flow from {{< math >}} $z_1^*$ {{< /math >}} back to {{< math >}} $\mu_{z_1}$ {{< /math >}} and {{< math >}} $\Sigma_{z_1}$ {{< /math >}} using the reparameterization trick.

Imagine a highly simplified scenario for a single dimension of {{< math >}} $z_1$ {{< /math >}}:

We're focusing on just one dimension of {{< math >}} $z_1$ {{< /math >}}, say {{< math >}} $z_{1,\text{dim1}}$ {{< /math >}}.

The encoder network for Cell 1, given its {{< math >}} $(u_1, s_1)$ {{< /math >}} data, outputs:

- **Mean**: {{< math >}} $\mu_{z_{1,\text{dim1}}} = 0.5$ {{< /math >}}
- **Log-variance**: {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}} = -1.0$ {{< /math >}} (This implies {{< math >}} $\Sigma_{z_{1,\text{dim1}}} = e^{-1.0} \approx 0.3679$ {{< /math >}})
- **Standard deviation**: {{< math >}} $\sigma_{z_{1,\text{dim1}}} = \sqrt{e^{-1.0}} = e^{-0.5} \approx 0.6065$ {{< /math >}}

We sample a noise value from a standard normal distribution: {{< math >}} $\epsilon_{1,\text{dim1}} = 0.8$ {{< /math >}} (a single random draw).

### Forward Pass (Calculating {{< math >}} $z_1^*$ {{< /math >}} and its contribution to Loss)

**Calculate {{< math >}} $z_1^*$ {{< /math >}} using reparameterization**:
{{< math >}} $$z_{1,\text{dim1}}^* = \mu_{z_{1,\text{dim1}}} + \sigma_{z_{1,\text{dim1}}} \times \epsilon_{1,\text{dim1}}$$ {{< /math >}}
{{< math >}} $$z_{1,\text{dim1}}^* = 0.5 + 0.6065 \times 0.8$$ {{< /math >}}
{{< math >}} $$z_{1,\text{dim1}}^* = 0.5 + 0.4852 = 0.9852$$ {{< /math >}}

**Downstream Computation & Loss**: This {{< math >}} $z_{1,\text{dim1}}^*$ {{< /math >}} value (along with other dimensions of {{< math >}} $z_1$ {{< /math >}} and other factors) goes into the decoder, kinetic equations, and ultimately contributes to the overall loss {{< math >}} $L_{\text{velo}}$ {{< /math >}}.

**Assume we have the gradient of the loss with respect to {{< math >}} $z_1^*$ {{< /math >}}**: After the forward pass and computing the total loss {{< math >}} $L_{\text{velo}}$ {{< /math >}}, backpropagation would compute the gradient of the loss with respect to this sampled latent variable.

Let's assume, for simplicity, that:
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial z_{1,\text{dim1}}^*} = 2.0$$ {{< /math >}}

(This means if {{< math >}} $z_{1,\text{dim1}}^*$ {{< /math >}} were to increase by 1 unit, the loss would increase by 2.0 units, implying we want to decrease {{< math >}} $z_{1,\text{dim1}}^*$ {{< /math >}} to reduce loss).

### Backward Pass (Gradient Flow to {{< math >}} $\mu_{z_1}$ {{< /math >}} and {{< math >}} $\Sigma_{z_1}$ {{< /math >}} parameters)

Now, we use the chain rule to propagate this gradient back to {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}} and {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}}.

#### 1. Gradient for {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}}

We know {{< math >}} $z_{1,\text{dim1}}^* = \mu_{z_{1,\text{dim1}}} + \sigma_{z_{1,\text{dim1}}} \times \epsilon_{1,\text{dim1}}$ {{< /math >}}.

The partial derivative of {{< math >}} $z_{1,\text{dim1}}^*$ {{< /math >}} with respect to {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}} is:
{{< math >}} $$\frac{\partial z_{1,\text{dim1}}^*}{\partial \mu_{z_{1,\text{dim1}}}} = 1$$ {{< /math >}}

Now, applying the chain rule:
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \mu_{z_{1,\text{dim1}}}} = \frac{\partial L_{\text{velo}}}{\partial z_{1,\text{dim1}}^*} \times \frac{\partial z_{1,\text{dim1}}^*}{\partial \mu_{z_{1,\text{dim1}}}}$$ {{< /math >}}
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \mu_{z_{1,\text{dim1}}}} = 2.0 \times 1 = 2.0$$ {{< /math >}}

**Interpretation**: A positive gradient of 2.0 for {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}} means that increasing {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}} would increase the loss. Therefore, the optimizer will decrease {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}} in the update step.

#### 2. Gradient for {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}}

The partial derivative of {{< math >}} $z_{1,\text{dim1}}^*$ {{< /math >}} with respect to {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}} is:
{{< math >}} $$\frac{\partial z_{1,\text{dim1}}^*}{\partial \sigma_{z_{1,\text{dim1}}}} = \epsilon_{1,\text{dim1}}$$ {{< /math >}}
{{< math >}} $$\frac{\partial z_{1,\text{dim1}}^*}{\partial \sigma_{z_{1,\text{dim1}}}} = 0.8$$ {{< /math >}}

Applying the chain rule:
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \sigma_{z_{1,\text{dim1}}}} = \frac{\partial L_{\text{velo}}}{\partial z_{1,\text{dim1}}^*} \times \frac{\partial z_{1,\text{dim1}}^*}{\partial \sigma_{z_{1,\text{dim1}}}}$$ {{< /math >}}
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \sigma_{z_{1,\text{dim1}}}} = 2.0 \times 0.8 = 1.6$$ {{< /math >}}

**Interpretation**: A positive gradient of 1.6 for {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}} means increasing {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}} would increase the loss. So, the optimizer will decrease {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}}.

#### 3. Gradient for {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}} (if encoder outputs log-variance)

Often, neural networks output {{< math >}} $\log \Sigma$ {{< /math >}} for numerical stability and to ensure the variance is positive. We need one more chain rule step:

We know {{< math >}} $\sigma_{z_{1,\text{dim1}}} = e^{0.5 \times \log \Sigma_{z_{1,\text{dim1}}}}$ {{< /math >}}.

The partial derivative of {{< math >}} $\sigma_{z_{1,\text{dim1}}}$ {{< /math >}} with respect to {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}} is:
{{< math >}} $$\frac{\partial \sigma_{z_{1,\text{dim1}}}}{\partial \log \Sigma_{z_{1,\text{dim1}}}} = 0.5 \times e^{0.5 \times \log \Sigma_{z_{1,\text{dim1}}}} = 0.5 \times \sigma_{z_{1,\text{dim1}}}$$ {{< /math >}}
{{< math >}} $$\frac{\partial \sigma_{z_{1,\text{dim1}}}}{\partial \log \Sigma_{z_{1,\text{dim1}}}} = 0.5 \times 0.6065 = 0.30325$$ {{< /math >}}

Applying the chain rule:
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \log \Sigma_{z_{1,\text{dim1}}}} = \frac{\partial L_{\text{velo}}}{\partial \sigma_{z_{1,\text{dim1}}}} \times \frac{\partial \sigma_{z_{1,\text{dim1}}}}{\partial \log \Sigma_{z_{1,\text{dim1}}}}$$ {{< /math >}}
{{< math >}} $$\frac{\partial L_{\text{velo}}}{\partial \log \Sigma_{z_{1,\text{dim1}}}} = 1.6 \times 0.30325 = 0.4852$$ {{< /math >}}

**Interpretation**: A positive gradient of 0.4852 for {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}} means increasing {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}} (and thus {{< math >}} $\Sigma_{z_{1,\text{dim1}}}$ {{< /math >}}) would increase the loss. So, the optimizer will decrease {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}}.

### Parameter Update (Conceptual)

Let's assume a learning_rate of 0.01.

- **New {{< math >}} $\mu_{z_{1,\text{dim1}}}$ {{< /math >}}**: {{< math >}} $0.5 - 0.01 \times 2.0 = 0.5 - 0.02 = 0.48$ {{< /math >}}
- **New {{< math >}} $\log \Sigma_{z_{1,\text{dim1}}}$ {{< /math >}}**: {{< math >}} $-1.0 - 0.01 \times 0.4852 = -1.0 - 0.004852 = -1.004852$ {{< /math >}}

In the next iteration, the encoder network's weights and biases would be updated based on these gradients (and the gradients from all other parts of the network and all other samples in the mini-batch). When Cell 1's data is fed through the encoder again, it will produce a new mean of 0.48 and a new log-variance of -1.004852 for this dimension of {{< math >}} $z_1$ {{< /math >}}, which should ideally lead to a lower overall loss.