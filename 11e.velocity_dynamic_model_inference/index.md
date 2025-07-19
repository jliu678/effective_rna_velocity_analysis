---
title: 🧬 Dynamic RNA velocity model-- (2) parameter inference 
summary: Here derives the mathematics underpinning the parameter inference of dynamic RNA velocity model, which is the second installment of our blog series to effectively apply the dynamic model in revealing the RNA velocity of single-cell RNAseq.
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, parameter inference, EM algorithm
  - Dynamic model
  - scVelo
  - Differential equations
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
## Introduction
To effectively apply the dynamic model for revealing RNA velocity in single-cell RNA-seq data, this second installment of our blog series takes a deep dive into its parameter inference using a two-stage EM algorithm. In this approach, latent time is initially assigned using an explicit formula, and then refined through standard optimization during the "Expectation" step of the final EM iteration.

## Overview
- section A) and B) derives the equations (9-11) in the [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3), which together give the explicit formula for latent time {{< math >}} $\tau$ {{< /math >}} in terms of unspliced RNA {{< math >}} $u$ {{< /math >}}, spliced RNA {{< math >}} $s$ {{< /math >}}, and kinetic parameters {{< math >}} $\alpha^{(k)}$ {{< /math >}}, {{< math >}} $\beta$ {{< /math >}}, and {{< math >}} $\gamma$ {{< /math >}}
- section C) explains why the explicit formula is less accurate than exact optimization
- section D) compares the EM implementation for scVelo parameter inference with classical Gaussian Mixture Model (GMM) EM algorithm, and proposes a sketch of EM for RNA velocity with learnable time prior {{< math >}} \pi_j {{< /math >}}.

## Symbol definitions
Please refer to the [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3) for the definitions of the symbols used in this derivation.

## A) Derive equations (9-10) in [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3)

Let's derive the relationship $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$ from the provided information.

### Given Information

#### System of Differential Equations:

{{< math >}}
$$
\frac{du(t)}{dt} = \alpha^{(k)} - \beta u(t)
$$
{{< /math >}}

{{< math >}}
$$
\frac{ds(t)}{dt} = \beta u(t) - \gamma s(t)
$$
{{< /math >}}

#### Solutions for u(t) and s(t):

{{< math >}}
$$
u(t) = u_0 e^{-\beta\tau} + \frac{\alpha^{(k)}}{\beta} (1 - e^{-\beta\tau})
$$
{{< /math >}}

(This can be rewritten as $u(t) = \frac{\alpha^{(k)}}{\beta} + (u_0 - \frac{\alpha^{(k)}}{\beta})e^{-\beta\tau}$)

{{< math >}}
$$
s(t) = s_0 e^{-\gamma\tau} + \frac{\alpha^{(k)}}{\gamma} (1 - e^{-\gamma\tau}) + \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) (e^{-\gamma\tau} - e^{-\beta\tau})
$$
{{< /math >}}

(This can be expanded to $s(t) = \frac{\alpha^{(k)}}{\gamma} + (s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta})e^{-\gamma\tau} - (\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta})e^{-\beta\tau}$)

Where $\tau = t - t_0$.

#### Steady-State Values:

From setting the derivatives to zero:

{{< math >}}
$$
u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}
$$
{{< /math >}}

{{< math >}}
$$
s_\infty = \frac{\beta u_\infty^{(k)}}{\gamma} = \frac{\alpha^{(k)}}{\gamma}
$$
{{< /math >}}

#### New Variable Definitions:

{{< math >}}
$$
\tilde{\beta} := \frac{\beta}{\gamma - \beta}
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}(t) := s(t) - \tilde{\beta}u(t)
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}_\infty := s_\infty - \tilde{\beta}u_\infty^{(k)}
$$
{{< /math >}}

### Derivation Steps

Our goal is to demonstrate: $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$.

#### 1. Express u(t) in a convenient form:

Let's rewrite u(t) using its steady-state value $u_\infty^{(k)}$:

{{< math >}}
$$
u(t) = u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}
$$
{{< /math >}}

#### 2. Substitute s(t) and u(t) into the definition of $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \tilde{\beta}\left[u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}\right]
$$
{{< /math >}}

Now, substitute the definition of $\tilde{\beta} = \frac{\beta}{\gamma - \beta}$ and $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \frac{\beta}{\gamma - \beta}\left[\frac{\alpha^{(k)}}{\beta} + \left(u_0 - \frac{\alpha^{(k)}}{\beta}\right)e^{-\beta\tau}\right]
$$
{{< /math >}}

Let's simplify the terms multiplied by $\tilde{\beta}$:

{{< math >}}
$$
\tilde{\beta}u_\infty^{(k)} = \frac{\beta}{\gamma - \beta} \cdot \frac{\alpha^{(k)}}{\beta} = \frac{\alpha^{(k)}}{\gamma - \beta}
$$
{{< /math >}}

{{< math >}}
$$
\tilde{\beta}(u_0 - u_\infty^{(k)}) = \frac{\beta}{\gamma - \beta}\left(u_0 - \frac{\alpha^{(k)}}{\beta}\right) = \frac{\beta}{\gamma - \beta}\left(\frac{\beta u_0 - \alpha^{(k)}}{\beta}\right) = \frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}
$$
{{< /math >}}

Substitute these back into the expression for $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left[\frac{\alpha^{(k)}}{\gamma} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau} - \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\beta\tau}\right] - \left[\frac{\alpha^{(k)}}{\gamma - \beta} + \left(\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}\right)e^{-\beta\tau}\right]
$$
{{< /math >}}

#### 3. Cancel the $e^{-\beta\tau}$ terms:

Observe the coefficients of $e^{-\beta\tau}$:

{{< math >}}
$$
-\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) - \left(\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta}\right)
$$
{{< /math >}}

Since $\frac{\beta u_0 - \alpha^{(k)}}{\gamma - \beta} = -\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}$, this simplifies to:

{{< math >}}
$$
-\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) - \left(-\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) = -\left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) + \left(\frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right) = 0
$$
{{< /math >}}

The terms containing $e^{-\beta\tau}$ indeed cancel out, as intended by the choice of $\tilde{\beta}$.

#### 4. Simplify $\tilde{s}(t)$:

{{< math >}}
$$
\tilde{s}(t) = \left(\frac{\alpha^{(k)}}{\gamma} - \frac{\alpha^{(k)}}{\gamma - \beta}\right) + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

Let's simplify the constant term:

{{< math >}}
$$
\frac{\alpha^{(k)}}{\gamma} - \frac{\alpha^{(k)}}{\gamma - \beta} = \alpha^{(k)}\left(\frac{1}{\gamma} - \frac{1}{\gamma - \beta}\right) = \alpha^{(k)}\left(\frac{(\gamma - \beta) - \gamma}{\gamma(\gamma - \beta)}\right) = \alpha^{(k)}\left(\frac{-\beta}{\gamma(\gamma - \beta)}\right) = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

So, we have:

{{< math >}}
$$
\tilde{s}(t) = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

#### 5. Calculate $\tilde{s}_0$ and $\tilde{s}_\infty$:

Given {{< math >}} $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$ {{< /math >}}, {{< math >}} $ s_\infty = \frac{\beta u_\infty^{(k)}}{\gamma} = \frac{\alpha^{(k)}}{\gamma}
$ {{< /math >}} and the definition of $\tilde{s}_\infty$ and $\tilde{\beta}$, we derive below. You may also consider view $\tilde{s}_\infty$  as Steady-state of $\tilde{s}(t)$  ($\tau \to \infty$, $e^{-\gamma\tau} \to 0$ ).

{{< math >}}
$$
\tilde{s}_\infty = -\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} (= \lim_{\tau \to \infty} \tilde{s}(t) )
$$
{{< /math >}}

**$\tilde{s}_0$ (Initial value of $\tilde{s}(t)$ at $\tau = 0$):**

{{< math >}}
$$
\tilde{s}_0 = s(t_0) - \tilde{\beta}u(t_0) = s_0 - \tilde{\beta}u_0 = s_0 - \frac{\beta}{\gamma - \beta}u_0
$$
{{< /math >}}

#### 6. Verify the target relationship: $\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}$

**Left-Hand Side (LHS):**

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} + \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}\right) - \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right)
$$
{{< /math >}}

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = \left(s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta}\right)e^{-\gamma\tau}
$$
{{< /math >}}

**Right-Hand Side (RHS):**

First, calculate $(\tilde{s}_0 - \tilde{s}_\infty)$:

{{< math >}}
$$
(\tilde{s}_0 - \tilde{s}_\infty) = \left(s_0 - \frac{\beta u_0}{\gamma - \beta}\right) - \left(-\frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right) = s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

Now, multiply by $e^{-\gamma\tau}$:

{{< math >}}
$$
(\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau} = \left(s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}\right)e^{-\gamma\tau}
$$
{{< /math >}}

Finally, let's compare the coefficients of $e^{-\gamma\tau}$ from both the LHS and RHS. We need to confirm that:

{{< math >}}
$$
s_0 - \frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta} \text{ equals } s_0 - \frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

Let's simplify the terms without $s_0$:

**From LHS coefficient:**

{{< math >}}
$$
-\frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)} - \beta u_0}{\gamma - \beta} = -\frac{\alpha^{(k)}}{\gamma} + \frac{\alpha^{(k)}}{\gamma - \beta} - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{1}{\gamma - \beta} - \frac{1}{\gamma}\right) - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{\gamma - (\gamma - \beta)}{\gamma(\gamma - \beta)}\right) - \frac{\beta u_0}{\gamma - \beta}
$$
$$
= \alpha^{(k)}\left(\frac{\beta}{\gamma(\gamma - \beta)}\right) - \frac{\beta u_0}{\gamma - \beta} = \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)} - \frac{\beta u_0}{\gamma - \beta}
$$
{{< /math >}}

**From RHS coefficient:**

{{< math >}}
$$
-\frac{\beta u_0}{\gamma - \beta} + \frac{\alpha^{(k)} \beta}{\gamma(\gamma - \beta)}
$$
{{< /math >}}

As you can see, the simplified coefficients from both sides are indeed identical.

### Conclusion

By introducing the carefully chosen auxiliary variable $\tilde{\beta} = \frac{\beta}{\gamma - \beta}$, the transformation $\tilde{s}(t) = s(t) - \tilde{\beta}u(t)$ strategically eliminates the $e^{-\beta\tau}$ dependence from the solution of s(t). This simplification reveals a clean, exponential decay relationship for $\tilde{s}(t)$ towards its new steady-state value $\tilde{s}_\infty$, as expressed by:

{{< math >}}
$$
\tilde{s}(t) - \tilde{s}_\infty = (\tilde{s}_0 - \tilde{s}_\infty)e^{-\gamma\tau}
$$
{{< /math >}}

This shows how a more complex system can be transformed into a simpler, recognizable form by defining appropriate auxiliary variables.

## B) Derive equation (11) in [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3)

Above smoothly gives Inverse Relationship for $\tau$ in equation (10) in the [scVelo paper](https://www.nature.com/articles/s41587-020-0591-3). Although equation (11) looks not obvious, it actually requires even simpler reasoning.

{{< math >}}
$$
\tau = -\frac{1}{\beta}\ln\left(\frac{u - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)
$$
{{< /math >}}

This expression makes sense given the solution for {{< math >}} $u(t)$ {{< /math >}} can be rewritten as below by substituting {{< math >}} $u_\infty^{(k)} = \frac{\alpha^{(k)}}{\beta}$ {{< /math >}}:

{{< math >}}
$$
u(t) = u_\infty^{(k)} + (u_0 - u_\infty^{(k)})e^{-\beta\tau}
$$
{{< /math >}}

If we solve this for {{< math >}} $\tau$ {{< /math >}}, we get:

* Rearrange: {{< math >}} $u(t) - u_\infty^{(k)} = (u_0 - u_\infty^{(k)})e^{-\beta\tau}$ {{< /math >}}

* Divide both sides by {{< math >}} $(u_0 - u_\infty^{(k)})$ {{< /math >}}: {{< math >}} $\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}} = e^{-\beta\tau}$ {{< /math >}}

* Take natural log: {{< math >}} $\ln\left(\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right) = -\beta\tau$ {{< /math >}}

* Solve for {{< math >}} $\tau$ {{< /math >}}: {{< math >}} $\tau = -\frac{1}{\beta}\ln\left(\frac{u(t) - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)$ {{< /math >}}

## C) Explicit Formula is Less Accurate than Exact Optimization

The difference in accuracy between the explicit formula and exact optimization comes down to **model assumptions versus reality**.

### Why the Explicit Formula is Less Accurate

#### 1. Single-Gene Assumption

The explicit formula:

{{< math >}}
$$\tau_i = -\frac{1}{\beta}\ln\left(\frac{u_i - u_\infty^{(k)}}{u_0 - u_\infty^{(k)}}\right)$$
{{< /math >}}

- **Uses only unspliced RNA** {{< math >}} $(u)$ {{< /math >}} from one gene
- **Ignores spliced RNA** {{< math >}} $(s)$ {{< /math >}} information completely
- **Assumes perfect adherence** to the theoretical kinetic model for that single gene

#### 2. No Cross-Gene Integration

- Each gene would give a **different time estimate** for the same cell
- The formula doesn't reconcile these conflicting estimates
- **Gene-specific noise** and measurement errors aren't averaged out

#### 3. Model Violations

Real cells don't perfectly follow the kinetic equations because:

- **Transcriptional bursting** creates noise
- **Cell-to-cell variability** in kinetic parameters
- **Measurement noise** in RNA counts
- **Model approximations** (e.g., constant degradation rates)

### Why Exact Optimization is More Accurate

#### 1. Multi-Gene Integration

The optimization objective:

{{< math >}}
$$\tau_i^* = \arg\min_\tau \sum_j \left[ (u_{i,j} - u_j(\tau))^2 + (s_{i,j} - s_j(\tau))^2 \right]$$
{{< /math >}}

- **Uses both unspliced AND spliced** RNA information
- **Integrates evidence across ALL genes** simultaneously
- **Finds the time that best explains the entire transcriptome**

#### 2. Noise Averaging

- Gene-specific measurement errors **cancel out** when averaged
- **Outlier genes** have less impact on the final time estimate
- **Statistical power increases** with more observations

#### 3. Comprehensive Model Fitting

- Accounts for **both RNA species** {{< math >}} $u$ {{< /math >}} and {{< math >}} $s$ {{< /math >}} jointly
- **Balances competing evidence** from different genes
- **Robust to individual gene model violations**

### Computational Trade-off

**Explicit Formula**: {{< math >}} $O(1)$ {{< /math >}} per cell-- Fast but uses limited information

**Exact Optimization**: {{< math >}} $O(n_{iter} \times n_{genes})$ {{< /math >}} per cell-- Slower but uses all available information

### Mathematical Perspective

#### Information Content Comparison

**Explicit Formula Information**:
{{< math >}}
$$I_{explicit} = f(u_{i,j}) \text{ for single gene } j$$
{{< /math >}}

**Exact Optimization Information**:
{{< math >}}
$$I_{exact} = f(\{u_{i,j}, s_{i,j}\}_{j=1}^{N_{genes}})$$
{{< /math >}}

Where {{< math >}} $I_{exact} >> I_{explicit}$ {{< /math >}} in terms of information content.

#### Error Propagation

For the explicit formula, measurement error in a single gene directly affects the time estimate:

{{< math >}}
$$\sigma_{\tau}^2 \propto \frac{\sigma_u^2}{(u_i - u_\infty)^2}$$
{{< /math >}}

For exact optimization, errors are averaged across genes:

{{< math >}}
$$\sigma_{\tau}^2 \propto \frac{1}{N_{genes}} \sum_j \frac{\sigma_{u,j}^2 + \sigma_{s,j}^2}{(u_{i,j} - u_j(\tau))^2 + (s_{i,j} - s_j(\tau))^2}$$
{{< /math >}}

### Biological Analogy

Think of it like **estimating someone's age**:

**Explicit Formula** = Looking at just their hair color
- Quick, but hair color alone is unreliable

**Exact Optimization** = Looking at hair, skin, posture, clothing style, speech patterns, etc.
- Takes more time, but gives a much better estimate by integrating multiple sources of evidence

### Why the Two-Stage Strategy Works

The two-stage approach in scVelo is brilliant because:

1. **Stage 1 (Explicit)**: Gets "close enough" quickly when parameters are changing rapidly
2. **Stage 2 (Exact)**: Refines to high accuracy when parameters have stabilized

#### Convergence Analysis

Early iterations: {{< math >}} $|\theta^{(t)} - \theta^{(t-1)}|$ {{< /math >}} is large
- Approximate times sufficient since parameters will change significantly anyway

Later iterations: {{< math >}} $|\theta^{(t)} - \theta^{(t-1)}| \approx 0$ {{< /math >}}
- Exact times needed for final parameter refinement

#### Computational Efficiency

Total computational cost:

{{< math >}}
$$\text{Cost} = T_{approx} \times O(1) + T_{exact} \times O(n_{iter} \times n_{genes})$$
{{< /math >}}

Where {{< math >}} $T_{approx} >> T_{exact}$ {{< /math >}}, making the overall algorithm much faster than using exact optimization throughout.

### Key Insight

The explicit formula serves as an excellent **initialization and approximation** tool, while exact optimization provides the **precision needed for final convergence**. This hybrid approach captures the best of both worlds: computational efficiency and high accuracy.


## D) EM implementation

### scVelo EM
```python
initialize theta = (α, β, γ)  
# Initialize model parameters:
# α: transcription rate, β: splicing rate, γ: degradation rate

while not converged:  
    # Repeat the EM steps until parameters stabilize

    # E-step
    mu_traj = [solve_ODE(t, theta) for t in time_grid]  
    # Simulate the ODE using current parameters to get the predicted gene expression (mu_t)
    # at each time point t in the predefined time grid

    for x_i in data:  
        # Iterate over each observed cell x_i (expression vector)

        for j, mu_t in enumerate(mu_traj):  
            # For each time point t_j and corresponding model prediction mu_t

            logliks[j] = -||x_i - mu_t||^2 / σ  
            # Compute (unnormalized) log-likelihood of x_i under the model at time t_j
            # Assuming Gaussian noise: log P(x_i | t_j) ∝ -||x_i - mu_t||² / σ

        q_i = softmax(logliks)  
        # Convert log-likelihoods into a probability distribution over time
        # This is the posterior q_i[j] = P(t_j | x_i, θ) for each time point t_j

    # M-step
    def loss_fn(theta):  
        # Define a loss function to optimize the model parameters θ

        loss = 0  
        # Initialize total loss

        for x_i, q_i in zip(data, q_all):  
            # Iterate over each cell x_i and its posterior time distribution q_i

            for j, t_j in enumerate(time_grid):  
                # For each possible time point t_j

                mu_tj = solve_ODE(t_j, theta)  
                # Predict expression using current θ at time t_j

                loss += q_i[j] * ||x_i - mu_tj||^2  
                # Add the squared difference weighted by q_i[j] (soft assignment)
                # This is the expected reconstruction error under the posterior q_i

        return loss  
        # Return total expected loss to minimize

    theta = optimize(loss_fn)  
    # Update model parameters by minimizing the loss using an optimizer (e.g. gradient descent)
```
Compared to the below GMM EM, the scVelo EM has not incorporated the time prior {{< math >}} \pi_j {{< /math >}}.
### GMM EM
```python
initialize mu_k, sigma_k, pi_k for k in 1..K  
# Initialize the parameters for each of the K Gaussian components:
#   - mu_k: mean of component k
#   - sigma_k: standard deviation of component k
#   - pi_k: prior probability (mixing coefficient) of component k
# These are randomly initialized or via k-means, and will be updated iteratively.

while not converged:  
# Repeat EM steps until the parameters stop changing significantly

    # E-step
     for x_i in data:  
        # Loop over each observed data point x_i

        for k in 1..K:  
            # For each Gaussian component k

            log_p_x_given_z = -((x_i - mu_k)^2) / (2 * sigma_k^2)  
            # Compute log-likelihood:
            # log P(x_i | z=k) under Gaussian N(mu_k, sigma_k^2), ignoring constants
            # This measures how well component k explains data point x_i

            log_p_z = log(pi_k)  
            # Log of prior probability of choosing component k

            logliks[k] = log_p_x_given_z + log_p_z  
            # Compute log posterior (up to normalization): log P(z=k | x_i)

        q_i = softmax(logliks)  
        # Normalize to get posterior probabilities:
        # q_i[k] = P(z=k | x_i), the "responsibility" of component k for x_i


    # M-step
    for k in 1..K:  
        # For each component k

        N_k = sum_i q_i[k]  
        # Effective number of points assigned to component k (soft count)

        mu_k = sum_i q_i[k] * x_i / N_k  
        # Update the mean of component k as the weighted average of data points,
        # weighted by how responsible k is for each point (q_i[k])

        sigma_k = sqrt(sum_i q_i[k] * (x_i - mu_k)^2 / N_k)  
        # Update standard deviation as the weighted standard deviation
        # of the data points around the new mean

        pi_k = N_k / N  
        # Update the mixing coefficient (prior) for component k
        # as the proportion of total responsibility assigned to k
```
Incorporating the time prior {{< math >}} \pi_j {{< /math >}} into the scVelo EM algorithm allows it to leverage global temporal structure, improving the model's ability to infer latent time points for each cell. Especially when biological markers are available to inform the temporal structure, the prior acts as a regularization term, guiding the model towards more plausible time assignments based on the overall distribution of cells across time points. I will test the effectiveness of this approach some time later, but below is the pseudo-code for the scVelo EM with learnable time prior {{< math >}} \pi_j {{< /math >}}.

### scVelo EM with Learnable Time Prior π_j 
```python
initialize theta = (α, β, γ)  
# Initialize model parameters:
# α: transcription rate
# β: splicing rate
# γ: degradation rate
# These govern the gene expression dynamics through an ODE model.

initialize π = [1/T] * T  # uniform over time grid  
# Initialize time prior π as a uniform distribution over T discrete time points.
# This reflects no prior knowledge about the timing of each cell initially.

while not converged:  
    # Repeat E-step and M-step until model parameters (θ and π) converge.

    # --- E-step: Estimate latent time distribution for each cell ---
    mu_traj = [solve_ODE(t, theta) for t in time_grid]  
    # Simulate the dynamic model using current parameters θ
    # This gives the expected gene expression μ(t) at each time point t in the grid.

    for x_i in data:  
        # Loop over each observed cell x_i (expression vector)

        for j, mu_t in enumerate(mu_traj):  
            # For each time point t_j and its predicted expression μ_t:

            logliks[j] = -||x_i - mu_t||^2 / σ + log(π[j])  
            # Compute an unnormalized log-posterior score for time point t_j:
            #
            #   - The first term, -||x_i - mu_t||² / σ, is the log-likelihood.
            #     It measures how well the model at time t_j explains cell x_i,
            #     assuming Gaussian noise with variance σ.
            #
            #   - The second term, log(π[j]), is the log of the prior probability
            #     of time t_j (learned in the previous M-step).
            #
            # Together, this computes:
            #
            #   log P(t_j | x_i) ∝ log P(x_i | t_j) + log π[j]
            #
            # This is a direct application of Bayes’ rule and balances:
            #   - Goodness-of-fit (how well the model explains the data)
            #   - Global structure (how likely time t_j is overall)
            #
            # This prevents overfitting rare time points that happen to fit x_i well.

        q_i = softmax(logliks)  
        # Normalize the log scores across all j using the softmax function
        # to get a posterior distribution q_i[j] = P(t_j | x_i, θ, π).
        # This soft assignment tells us how likely cell x_i came from each time point.

    # --- M-step: Update π and θ to maximize expected data fit ---

    π[j] = sum_i q_i[j] / N  for all j  
    # Update the prior over time points:
    # For each time point t_j, π[j] becomes the average of q_i[j] over all cells.
    # This reflects how many cells the model currently thinks came from time t_j.
    # The prior π is now data-driven instead of uniform.

    def loss_fn(theta):  
        # Define the objective function for re-estimating model parameters θ.

        loss = 0  
        # Initialize the total loss to accumulate over all cells and times.

        for x_i, q_i in zip(data, q_all):  
            # For each observed cell x_i and its soft time assignment q_i:

            for j, t_j in enumerate(time_grid):  
                # For each time point t_j in the grid:

                mu_tj = solve_ODE(t_j, theta)  
                # Simulate the ODE to predict expression at time t_j using θ.

                loss += q_i[j] * ||x_i - mu_tj||^2  
                # Add weighted squared error between observed x_i and predicted mu_tj.
                # The weight q_i[j] reflects the model’s belief that x_i came from t_j.
                # This gives a soft, expectation-based loss function.

        return loss  
        # Total expected reconstruction error across all cells and times.

    theta = optimize(loss_fn)  
    # Optimize model parameters θ to minimize the expected error (loss).
    # This is typically done using gradient descent or similar techniques.
```
