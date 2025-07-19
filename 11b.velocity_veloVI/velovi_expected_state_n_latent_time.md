---
title: "Velovi Expected State N Latent Time"
date: 2025-07-10
draft: True
---

# Expected Value of Kinetic State Probability Vector

This document explains how to calculate the expected value of the kinetic state probability vector {{< math >}}$\pi_{ng}${{< /math >}}, where this expectation is taken twice: first conditional on a latent variable {{< math >}}$z_n${{< /math >}}, and then averaged over the possible values of {{< math >}}$z_n${{< /math >}} for a given cell.

The main expression we're evaluating is:

{{< math >}}
$$
\mathbb{E}_{q_\phi(z_n|u_n,s_n)} \left[ \mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [\pi_{ng}] \right] \tag{1}
$$
{{< /math >}}

This represents VeloVI's average belief about the kinetic state probabilities ({{< math >}}$\pi_{ng}${{< /math >}}) for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$n${{< /math >}}, considering all possible latent representations ({{< math >}}$z_n${{< /math >}}) that could explain the observed data ({{< math >}}$u_n, s_n${{< /math >}}).

## Understanding the Components

- **{{< math >}}$\pi_{ng}${{< /math >}}**: The vector of probabilities for the 4 kinetic states: {{< math >}}$(\pi_{ng,1}, \pi_{ng,2}, \pi_{ng,3}, \pi_{ng,4})${{< /math >}}
- **{{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}**: The approximate posterior distribution for {{< math >}}$\pi_{ng}${{< /math >}} given a specific latent variable {{< math >}}$z_n${{< /math >}}. This is a Dirichlet distribution with concentration parameters {{< math >}}$\alpha_{q,ng}(z_n) = (\alpha_{q,ng,1}(z_n), \ldots, \alpha_{q,ng,K}(z_n))${{< /math >}}
- **{{< math >}}$q_\phi(z_n|u_n,s_n)${{< /math >}}**: The approximate posterior distribution for the latent variable {{< math >}}$z_n${{< /math >}} given the observed counts. This is a Gaussian distribution {{< math >}}$\mathcal{N}(\mu_{z_n}, \Sigma_{z_n})${{< /math >}}

## Step 1: Calculate the Inner Expectation

The inner expectation {{< math >}}$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [\pi_{ng}]${{< /math >}} is the expected value of a Dirichlet-distributed random variable.

If {{< math >}}$\pi_{ng} \sim \text{Dirichlet}(\alpha_{q,ng,1}(z_n), \ldots, \alpha_{q,ng,K}(z_n))${{< /math >}}, then the expected value of its {{< math >}}$i${{< /math >}}-th component is:

{{< math >}}
$$
\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [\pi_{ng,i}] = \frac{\alpha_{q,ng,i}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)} \tag{2}
$$
{{< /math >}}

The result of the inner expectation is a vector of expected probabilities:

{{< math >}}
$$
\bar{\pi}_{ng}(z_n) = \left( \frac{\alpha_{q,ng,1}(z_n)}{\sum_j \alpha_{q,ng,j}(z_n)}, \ldots, \frac{\alpha_{q,ng,K}(z_n)}{\sum_j \alpha_{q,ng,j}(z_n)} \right) \tag{3}
$$
{{< /math >}}

### Numeric Example for Cell 1, Gene A

Let's assume we sample {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}}.

The encoder outputs Dirichlet concentration parameters for Gene A:
{{< math >}}$\alpha_{q,1A}(z_1^*) = (5.0, 0.1, 0.2, 0.1)${{< /math >}}

Calculate the expected {{< math >}}$\pi_{1A}${{< /math >}} for this {{< math >}}$z_1^*${{< /math >}}:

- Sum of concentrations: {{< math >}}$5.0 + 0.1 + 0.2 + 0.1 = 5.4${{< /math >}}
- {{< math >}}$\mathbb{E}[\pi_{1A,1}|z_1^*] = 5.0/5.4 \approx 0.9259${{< /math >}}
- {{< math >}}$\mathbb{E}[\pi_{1A,2}|z_1^*] = 0.1/5.4 \approx 0.0185${{< /math >}}
- {{< math >}}$\mathbb{E}[\pi_{1A,3}|z_1^*] = 0.2/5.4 \approx 0.0370${{< /math >}}
- {{< math >}}$\mathbb{E}[\pi_{1A,4}|z_1^*] = 0.1/5.4 \approx 0.0185${{< /math >}}

So {{< math >}}$\bar{\pi}_{1A}(z_1^*) = (0.9259, 0.0185, 0.0370, 0.0185)${{< /math >}}.

## Step 2: Calculate the Outer Expectation

Now we need to average {{< math >}}$\bar{\pi}_{ng}(z_n)${{< /math >}} over the distribution {{< math >}}$q_\phi(z_n|u_n,s_n) = \mathcal{N}(\mu_{z_n}, \Sigma_{z_n})${{< /math >}}. Since {{< math >}}$\bar{\pi}_{ng}(z_n)${{< /math >}} involves a neural network mapping, this expectation is intractable and we use Monte Carlo approximation.

### Monte Carlo Approximation

For each sample {{< math >}}$m \in \{1, \ldots, M\}${{< /math >}}:

1. Draw random noise: {{< math >}}$\epsilon_z^{(m)} \sim \mathcal{N}(0, I_d)${{< /math >}}
2. Calculate sample: {{< math >}}$z_n^{*(m)} = \mu_{z_n} + \sqrt{\Sigma_{z_n}} \odot \epsilon_z^{(m)}${{< /math >}} (reparameterization trick)

### Numeric Example Continued

**Known parameters for Cell 1:**
- {{< math >}}$\mu_{z_1} = (0.5, -0.3)${{< /math >}}
- {{< math >}}$\log\Sigma_{z_1} = (-0.5, -0.5) \rightarrow \Sigma_{z_1} = \text{diag}(0.6065, 0.6065)${{< /math >}}

**Sample 1:**
- {{< math >}}$\epsilon_z^{(1)} = (0.2, -0.1)${{< /math >}}
- {{< math >}}$z_1^{*(1)} = (0.5, -0.3) + (0.7788, 0.7788) \odot (0.2, -0.1) = (0.6558, -0.3779)${{< /math >}}

**Sample 2:**
- {{< math >}}$\epsilon_z^{(2)} = (-0.5, 0.3)${{< /math >}}
- {{< math >}}$z_1^{*(2)} = (0.5, -0.3) + (0.7788, 0.7788) \odot (-0.5, 0.3) = (0.1106, -0.0663)${{< /math >}}

For {{< math >}}$z_1^{*(2)}${{< /math >}}, assume encoder outputs {{< math >}}$\alpha_{q,1A}(z_1^{*(2)}) = (4.0, 0.3, 0.1, 0.1)${{< /math >}}:

- Sum of concentrations: {{< math >}}$4.0 + 0.3 + 0.1 + 0.1 = 4.5${{< /math >}}
- {{< math >}}$\bar{\pi}_{1A}^{(2)} \approx (0.8889, 0.0667, 0.0222, 0.0222)${{< /math >}}

**Final average (with M=2):**

{{< math >}}
$$
\mathbb{E}_{q_\phi(z_n|u_n,s_n)} [\bar{\pi}_{ng}(z_n)] \approx \frac{1}{M} \sum_{m=1}^M \bar{\pi}_{ng}^{(m)} \tag{4}
$$
{{< /math >}}

- Component 1: {{< math >}}$(0.9259 + 0.8889)/2 = 0.9074${{< /math >}}
- Component 2: {{< math >}}$(0.0185 + 0.0667)/2 = 0.0426${{< /math >}}
- Component 3: {{< math >}}$(0.0370 + 0.0222)/2 = 0.0296${{< /math >}}
- Component 4: {{< math >}}$(0.0185 + 0.0222)/2 = 0.0204${{< /math >}}

**Final result:** {{< math >}}$\approx (0.9074, 0.0426, 0.0296, 0.0204)${{< /math >}}

---

# Expected Latent Time Calculation

Now we calculate the average latent time {{< math >}}$t_{ng}${{< /math >}} for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$n${{< /math >}}, accounting for uncertainty in both the cell's latent representation and the kinetic state.

The expression is:

{{< math >}}
$$
\mathbb{E}_{q_\phi(z_n|u_n,s_n)} \left[ \mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [t^{(k_{ng})}_{ng}] \right] \tag{5}
$$
{{< /math >}}

where {{< math >}}$t^{(k_{ng})}_{ng}${{< /math >}} represents the latent time for gene {{< math >}}$g${{< /math >}} in cell {{< math >}}$n${{< /math >}} if it's in kinetic state {{< math >}}$k_{ng}${{< /math >}}.

## Derive the Inner Expectation for Time

The inner expectation over kinetic states is:

{{< math >}}
$$
\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [t^{(k_{ng})}_{ng}] = \sum_{k=1}^K P(k_{ng} = k | z_n) \cdot t_{ng}^{(k)}(z_n) \tag{6}
$$
{{< /math >}}

Where:
- {{< math >}}$P(k_{ng} = k | z_n) = \frac{\alpha_{q,ng,k}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)}${{< /math >}}
- {{< math >}}$t_{ng}^{(k)}(z_n) = \rho_{ng}^{(k)}(z_n) \times t_g^s${{< /math >}} (decoder output times switching time)


In our expression {{< math >}}$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [t^{(k_{ng})}_{ng}]${{< /math >}}:

- **The random variable:** The random variable is {{< math >}}$t^{(k_{ng})}_{ng}${{< /math >}}. This variable's value *depends* on the discrete kinetic state {{< math >}}$k_{ng}${{< /math >}}.

- **The Possible Values of the Random Variable**
  - The kinetic state {{< math >}}$k_{ng}${{< /math >}} can take on {{< math >}}$K${{< /math >}} discrete values (1, 2, 3, 4 for VeloVI's four states).
  - For each specific state {{< math >}}$k${{< /math >}}, the corresponding value of our random variable is {{< math >}}$t_{ng}^{(k)}(z_n)${{< /math >}}.
  - This notation {{< math >}}$t_{ng}^{(k)}(z_n)${{< /math >}} explicitly shows that the latent time is specific to state {{< math >}}$k${{< /math >}} and is determined by {{< math >}}$z_n${{< /math >}}.

- **The Probability Distribution**
  - The expectation is taken with respect to {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}.
  - This distribution tells us the probabilities of being in each kinetic state {{< math >}}$k${{< /math >}}.
  - Specifically, {{< math >}}$\pi_{ng}${{< /math >}} is a vector of probabilities for the states, and {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}} is a Dirichlet distribution that *governs* these probabilities.

We know that {{< math >}}$k_{ng} \sim \text{Categorical}(\pi_{ng})${{< /math >}}. This means, if we knew the exact value of {{< math >}}$\pi_{ng}${{< /math >}}, then {{< math >}}$P(k_{ng}=k|\pi_{ng}) = \pi_{ng,k}${{< /math >}}.

However, {{< math >}}$\pi_{ng}${{< /math >}} itself is a random variable distributed according to {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}.

So, to get the marginal probability {{< math >}}$P(k_{ng}=k|z_n)${{< /math >}}, we need to integrate over {{< math >}}$\pi_{ng}${{< /math >}}:

{{< math >}}
$$
P(k_{ng}=k|z_n) = \int P(k_{ng}=k|\pi_{ng}) q_\phi(\pi_{ng}|z_n) d\pi_{ng} \tag{6.1}
$$
{{< /math >}}

{{< math >}}
$$
P(k_{ng}=k|z_n) = \int \pi_{ng,k} q_\phi(\pi_{ng}|z_n) d\pi_{ng} \tag{6.2}
$$
{{< /math >}}

This integral is simply the definition of the expected value of the {{< math >}}$k${{< /math >}}-th component of {{< math >}}$\pi_{ng}${{< /math >}} under the distribution {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}:

{{< math >}}
$$
P(k_{ng}=k|z_n) = \mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [\pi_{ng,k}] \tag{6.3}
$$
{{< /math >}}

For a Dirichlet distribution {{< math >}}$\text{Dir}(\alpha_{q,ng,1}(z_n), \ldots, \alpha_{q,ng,K}(z_n))${{< /math >}}, this expectation has a closed form:

{{< math >}}
$$
\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [\pi_{ng,k}] = \frac{\alpha_{q,ng,k}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)} \tag{6.4}
$$
{{< /math >}}

So, {{< math >}}$P(k_{ng}=k|z_n)${{< /math >}} is precisely this fraction.

Therefore:

{{< math >}}
$$
\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [t^{(k_{ng})}_{ng}] = \sum_{k=1}^K \left( t_{ng}^{(k)}(z_n) \cdot P(k_{ng}=k|z_n) \right) \tag{6}
$$
{{< /math >}}

Substituting equation (4) into equation (6):

{{< math >}}
$$
\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} [t^{(k_{ng})}_{ng}] = \sum_{k=1}^K \left( t_{ng}^{(k)}(z_n) \cdot \frac{\alpha_{q,ng,k}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)} \right) \tag{6.5}
$$
{{< /math >}}

This is exactly the equation we use in practice, explicitly showing that the expectation is a **weighted sum** where:
- Each possible latent time {{< math >}}$t_{ng}^{(k)}(z_n)${{< /math >}} is weighted by 
- The probability of its corresponding state {{< math >}}$k${{< /math >}} occurring, as determined by the approximate posterior {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}.

## Numeric Example for Expected Time

This expression asks for the average latent time {{< math >}}$t_{ng}${{< /math >}} for gene g in cell n, taking into account both:

- The uncertainty in the cell's underlying latent representation {{< math >}}$z_n${{< /math >}}
- The uncertainty in which kinetic state {{< math >}}$k_{ng}${{< /math >}} gene g is in

Similar to the previous expectation, this involves an intractable integral, so we will use Monte Carlo sampling for approximation.

The expression is:

{{< math >}}
$$\mathbb{E}_{q_\phi(z_n|u_n,s_n)} \left[\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} \left[t(k_{ng})_{ng}\right]\right] \tag{1}$$
{{< /math >}}

Where {{< math >}}$t(k_{ng})_{ng}${{< /math >}} represents the latent time for gene g in cell n if it's in kinetic state {{< math >}}$k_{ng}${{< /math >}}. Note that this latent time {{< math >}}$t_{ng}${{< /math >}} is itself a function of {{< math >}}$z_n${{< /math >}} and the specific kinetic state k.

Let's break down the calculation in steps with a numeric example for Cell 1 and Gene A.

### What We Assume We Have (after VeloVI training)
---
**Observed Data for Cell 1, Gene A:** {{< math >}}$u_{1A} = 0.8, s_{1A} = 0.2${{< /math >}}

**Learned Biophysical Parameters (θ) for Gene A:**
- Switching Time {{< math >}}$(t_A^s) = 8.0${{< /math >}}
- (Other θ like α,β,γ are not directly used in this specific calculation of expected time, but they would be used to calculate expected RNA levels from this time)

**Learned Neural Network Parameters (φ) for Encoder and Decoder:**

*Encoder outputs for* {{< math >}}$z_1${{< /math >}} *(given* {{< math >}}$u_1, s_1${{< /math >}}*):* Defines {{< math >}}$q_\phi(z_1|u_1,s_1) = \mathcal{N}(\mu_{z_1}, \Sigma_{z_1})${{< /math >}}
- {{< math >}}$\mu_{z_1} = (0.5, -0.3)${{< /math >}}
- {{< math >}}$\Sigma_{z_1} = \text{diag}(e^{-0.5}, e^{-0.5}) \approx \text{diag}(0.6065, 0.6065)${{< /math >}} (for {{< math >}}$z_n${{< /math >}} as 2D)

*Encoder outputs for* {{< math >}}$\pi_{1A}${{< /math >}} *(given* {{< math >}}$z_1${{< /math >}}*):* Defines {{< math >}}$q_\phi(\pi_{1A}|z_1) = \text{Dirichlet}(\alpha_{q,1A}(z_1))${{< /math >}}

*Decoder outputs for* {{< math >}}$\rho_{1A}^{(k)}${{< /math >}} *(given* {{< math >}}$z_1${{< /math >}}*):* For each state k, the decoder outputs a scaled latent time {{< math >}}$\rho_{1A}^{(k)}${{< /math >}}, where {{< math >}}$t_{1A}^{(k)} = \rho_{1A}^{(k)} \times t_A^s${{< /math >}}.

### Step 1: Calculate the Inner Expectation
---
#### **Recall the Inner Expectation**
{{< math >}}
$$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} \left[t(k_{ng})_{ng}\right] \tag{2}$$
{{< /math >}}

This inner expectation is over the kinetic states {{< math >}}$k_{ng}${{< /math >}}, weighted by their probabilities. The latent time {{< math >}}$t(k_{ng})_{ng}${{< /math >}} is specific to the state {{< math >}}$k_{ng}${{< /math >}}. We can write this expectation as a sum:

{{< math >}}
$$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} \left[t(k_{ng})_{ng}\right] = \sum_{k=1}^K P(k_{ng} = k|z_n) \cdot t_{ng}^{(k)}(z_n) \tag{3}$$
{{< /math >}}

Where:

- K is the number of kinetic states (4 in VeloVI)
- {{< math >}}$P(k_{ng} = k|z_n)${{< /math >}} is the probability that gene g in cell n is in state k, given {{< math >}}$z_n${{< /math >}}. This is the k-th component of the expected {{< math >}}$\pi_{ng}${{< /math >}} vector from the Dirichlet distribution: 

{{< math >}}
$$P(k_{ng} = k|z_n) = \frac{\alpha_{q,ng,k}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)} \tag{4}$$
{{< /math >}}

- {{< math >}}$t_{ng}^{(k)}(z_n)${{< /math >}} is the latent time for gene g in cell n if it is in state k. This is a deterministic output of the decoder, given {{< math >}}$z_n${{< /math >}}: {{< math >}}$t_{ng}^{(k)}(z_n) = \text{Decoder}_k(z_n) \times t_g^s = \rho_{ng}^{(k)}(z_n) \times t_g^s${{< /math >}}.

The result of this inner expectation will be a single scalar value, let's call it {{< math >}}$\bar{t}_{ng}(z_n)${{< /math >}}, which is still a function of {{< math >}}$z_n${{< /math >}}.

#### **Numeric Example for Cell 1, Gene A**
Here we take Cell 1, Gene A for a Numeric Example: assume we have a sample {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}} (we'll show how to get this in Step 2).

1. **From** {{< math >}}$z_1^*${{< /math >}}**, get Dirichlet concentration parameters:**
   - Assume encoder outputs: {{< math >}}$\alpha_{q,1A}(z_1^*) = (5.0, 0.1, 0.2, 0.1)${{< /math >}}
   - Sum of concentrations = 5.4

2. **Calculate** {{< math >}}$P(k_{1A} = k|z_1^*)${{< /math >}} **for each state:**
   - {{< math >}}$P(k_{1A} = 1|z_1^*) = 5.0/5.4 \approx 0.9259${{< /math >}}
   - {{< math >}}$P(k_{1A} = 2|z_1^*) = 0.1/5.4 \approx 0.0185${{< /math >}}
   - {{< math >}}$P(k_{1A} = 3|z_1^*) = 0.2/5.4 \approx 0.0370${{< /math >}}
   - {{< math >}}$P(k_{1A} = 4|z_1^*) = 0.1/5.4 \approx 0.0185${{< /math >}}

3. **From** {{< math >}}$z_1^*${{< /math >}}**, get latent times** {{< math >}}$t_{1A}^{(k)}(z_1^*)${{< /math >}} **for each state:**
   (These are specific outputs from the decoder for this {{< math >}}$z_1^*${{< /math >}})
   
   - Assume decoder outputs {{< math >}}$\rho_{1A}^{(1)}(z_1^*) = 0.5${{< /math >}} (for Induction)
     {{< math >}}$t_{1A}^{(1)}(z_1^*) = 0.5 \times t_A^s = 0.5 \times 8.0 = 4.0${{< /math >}}
   - Assume decoder outputs {{< math >}}$\rho_{1A}^{(2)}(z_1^*) = 0.3${{< /math >}} (for Repression)
     {{< math >}}$t_{1A}^{(2)}(z_1^*) = 0.3 \times t_A^s = 0.3 \times 8.0 = 2.4${{< /math >}}
   - Assume decoder outputs {{< math >}}$\rho_{1A}^{(3)}(z_1^*) = 0.8${{< /math >}} (for Steady-State Induction)
     {{< math >}}$t_{1A}^{(3)}(z_1^*) = 0.8 \times t_A^s = 0.8 \times 8.0 = 6.4${{< /math >}}
   - Assume decoder outputs {{< math >}}$\rho_{1A}^{(4)}(z_1^*) = 0.1${{< /math >}} (for Steady-State Repression)
     {{< math >}}$t_{1A}^{(4)}(z_1^*) = 0.1 \times t_A^s = 0.1 \times 8.0 = 0.8${{< /math >}}

4. **Calculate the weighted sum** {{< math >}}$\bar{t}_{1A}(z_1^*)${{< /math >}}**:**

{{< math >}}
$$\bar{t}_{1A}(z_1^*) = (0.9259 \times 4.0) + (0.0185 \times 2.4) + (0.0370 \times 6.4) + (0.0185 \times 0.8) \tag{5}$$
{{< /math >}}

{{< math >}}
$$\bar{t}_{1A}(z_1^*) = 3.7036 + 0.0444 + 0.2368 + 0.0148 \tag{6}$$
{{< /math >}}

{{< math >}}
$$\bar{t}_{1A}(z_1^*) \approx 3.9996 \text{ (This is a scalar value for this specific } z_1^*) \tag{7}$$
{{< /math >}}

### Step 2: Monte Carlo Approximation Calculates the Outer Expectation
---
{{< math >}}
$$\mathbb{E}_{q_\phi(z_n|u_n,s_n)} \left[\bar{t}_{ng}(z_n)\right] \tag{8}$$
{{< /math >}}

Now we need to average this {{< math >}}$\bar{t}_{ng}(z_n)${{< /math >}} value over the distribution of {{< math >}}$z_n${{< /math >}}. Since {{< math >}}$\bar{t}_{ng}(z_n)${{< /math >}} is a complex, non-linear function of {{< math >}}$z_n${{< /math >}} (involving neural network outputs), we use Monte Carlo approximation.

1. **Draw M samples of** {{< math >}}$z_n^*${{< /math >}} **from** {{< math >}}$q_\phi(z_n|u_n,s_n)${{< /math >}}**:**

   Recall {{< math >}}$q_\phi(z_1|u_1,s_1) = \mathcal{N}(\mu_{z_1}, \Sigma_{z_1})${{< /math >}} with {{< math >}}$\mu_{z_1} = (0.5, -0.3)${{< /math >}} and {{< math >}}$\Sigma_{z_1} = \text{diag}(0.6065, 0.6065)${{< /math >}}.

2. **For each sample** {{< math >}}$m \in \{1, ..., M\}${{< /math >}}**:**
   - Draw a random noise vector {{< math >}}$\varepsilon_z^{(m)} \sim \mathcal{N}(0, I_d)${{< /math >}}
   - Calculate {{< math >}}$z_n^{*(m)} = \mu_{z_n} + \sqrt{\Sigma_{z_n}} \odot \varepsilon_z^{(m)}${{< /math >}} (This is the reparameterization trick)

Just two samples for brevity, usually M=1000+:

*Sample 1:*
- {{< math >}}$\varepsilon_z^{(1)} = (0.2, -0.1)${{< /math >}}
- {{< math >}}$z_1^{*(1)} = (0.5, -0.3) + (\sqrt{0.6065}, \sqrt{0.6065}) \odot (0.2, -0.1) = (0.6558, -0.3779)${{< /math >}}

*Sample 2:*
- {{< math >}}$\varepsilon_z^{(2)} = (-0.5, 0.3)${{< /math >}}
- {{< math >}}$z_1^{*(2)} = (0.5, -0.3) + (\sqrt{0.6065}, \sqrt{0.6065}) \odot (-0.5, 0.3) = (0.1106, -0.0663)${{< /math >}}

3. **For each sampled** {{< math >}}$z_n^{*(m)}${{< /math >}}**, calculate** {{< math >}}$\bar{t}_{ng}(z_n^{*(m)})${{< /math >}}**:**

   Follow the steps for the inner expectation (from above) for each {{< math >}}$z_n^{*(m)}${{< /math >}}.


*For* {{< math >}}$z_1^{*(1)} = (0.6558, -0.3779)${{< /math >}}*:*
We already calculated: {{< math >}}$\bar{t}_{1A}^{(1)} \approx 3.9996${{< /math >}}

*For* {{< math >}}$z_1^{*(2)} = (0.1106, -0.0663)${{< /math >}}*:*
- Assume encoder outputs {{< math >}}$\alpha_{q,1A}(z_1^{*(2)}) = (4.0, 0.3, 0.1, 0.1)${{< /math >}} (differs from first sample). Sum = 4.5
- {{< math >}}$P(k_{1A} = 1|z_1^{*(2)}) = 4.0/4.5 \approx 0.8889${{< /math >}}
- {{< math >}}$P(k_{1A} = 2|z_1^{*(2)}) = 0.3/4.5 \approx 0.0667${{< /math >}}
- {{< math >}}$P(k_{1A} = 3|z_1^{*(2)}) = 0.1/4.5 \approx 0.0222${{< /math >}}
- {{< math >}}$P(k_{1A} = 4|z_1^{*(2)}) = 0.1/4.5 \approx 0.0222${{< /math >}}

Assume decoder outputs (different for this {{< math >}}$z_1^{*(2)}${{< /math >}}):
- {{< math >}}$\rho_{1A}^{(1)}(z_1^{*(2)}) = 0.4 \rightarrow t_{1A}^{(1)}(z_1^{*(2)}) = 0.4 \times 8.0 = 3.2${{< /math >}}
- {{< math >}}$\rho_{1A}^{(2)}(z_1^{*(2)}) = 0.2 \rightarrow t_{1A}^{(2)}(z_1^{*(2)}) = 0.2 \times 8.0 = 1.6${{< /math >}}
- {{< math >}}$\rho_{1A}^{(3)}(z_1^{*(2)}) = 0.7 \rightarrow t_{1A}^{(3)}(z_1^{*(2)}) = 0.7 \times 8.0 = 5.6${{< /math >}}
- {{< math >}}$\rho_{1A}^{(4)}(z_1^{*(2)}) = 0.05 \rightarrow t_{1A}^{(4)}(z_1^{*(2)}) = 0.05 \times 8.0 = 0.4${{< /math >}}

Calculate weighted sum {{< math >}}$\bar{t}_{1A}(z_1^{*(2)})${{< /math >}} for this sample:

{{< math >}}
$$\bar{t}_{1A}(z_1^{*(2)}) = (0.8889 \times 3.2) + (0.0667 \times 1.6) + (0.0222 \times 5.6) + (0.0222 \times 0.4) \tag{9}$$
{{< /math >}}

{{< math >}}
$$\bar{t}_{1A}(z_1^{*(2)}) = 2.8445 + 0.1067 + 0.1243 + 0.0089 \tag{10}$$
{{< /math >}}

{{< math >}}
$$\bar{t}_{1A}(z_1^{*(2)}) \approx 3.0844 \tag{11}$$
{{< /math >}}

4. **Average the results:**

The overall expectation is approximated by the average of these M scalar values:

{{< math >}}
$$\mathbb{E}_{q_\phi(z_n|u_n,s_n)} [\bar{t}_{ng}(z_n)] \approx \frac{1}{M} \sum_{m=1}^M \bar{t}_{ng}^{(m)} \tag{12}$$
{{< /math >}}

For our example with M=2:

Average latent time for Gene A in Cell 1:

{{< math >}}
$$\frac{3.9996 + 3.0844}{2} = 3.542 \tag{13}$$
{{< /math >}}

This final scalar value (e.g., 3.542) represents VeloVI's average estimated latent time for Gene A in Cell 1, taking into account all the uncertainties modeled by the variational distributions.

## Summary Table

| Item | Status | Role | Example Value |
|------|--------|------|---------------|
| {{< math >}}$u_n, s_n${{< /math >}} | Known (observed) | Input data | Cell 1 observations |
| {{< math >}}$\phi${{< /math >}} | Learned & Fixed | Neural network parameters | (Internal to networks) |
| {{< math >}}$\mu_{z_n}, \Sigma_{z_n}${{< /math >}} | Calculated by encoder | Parameters of {{< math >}}$q_\phi(z_n\|u_n,s_n)${{< /math >}} | {{< math >}}$(0.5, -0.3)${{< /math >}}, {{< math >}}$\text{diag}(0.6065, 0.6065)${{< /math >}} |
| {{< math >}}$z_n^{*(m)}${{< /math >}} | Sampled | Sample from {{< math >}}$q_\phi(z_n\|u_n,s_n)${{< /math >}} | {{< math >}}$(0.6558, -0.3779)${{< /math >}} |
| {{< math >}}$\alpha_{q,ng}(z_n^{*(m)})${{< /math >}} | Calculated by encoder | Dirichlet parameters | {{< math >}}$(5.0, 0.1, 0.2, 0.1)${{< /math >}} |
| {{< math >}}$\bar{\pi}_{ng}(z_n^{*(m)})${{< /math >}} | Calculated | Inner expectation result | {{< math >}}$(0.9259, 0.0185, 0.0370, 0.0185)${{< /math >}} |
| Final average | Calculated | Overall approximation | {{< math >}}$(0.9074, 0.0426, 0.0296, 0.0204)${{< /math >}} |


# Calculate Expected Velocity for a sampled {{< math >}}$z_n${{< /math >}}
> **Note: Posterior Velocity Distribution, i.e. the full posterior predictive velocity distribution, is just the distribution of all the expected velocity values for the sample {{< math >}}$z_n^*${{< /math >}} (drawn from {{< math >}}$q_\phi(z_n|u_n, s_n)${{< /math >}}).**

---

{{< math >}}
$$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} \left[v^{(g)}(t(k_{ng})_{ng}, k_{ng})\right] \tag{1}$$
{{< /math >}}

This term represents the expected velocity for gene g in cell n, given a specific latent variable {{< math >}}$z_n${{< /math >}}, averaging over the uncertainty in its kinetic state {{< math >}}$k_{ng}${{< /math >}}.

We will use the same logic as the previous expectation for latent time, as the structure is identical: we are calculating the expectation of a quantity ({{< math >}}$v^{(g)}${{< /math >}}) that depends on a discrete random variable ({{< math >}}$k_{ng}${{< /math >}}), whose probabilities are determined by the Dirichlet distribution {{< math >}}$q_\phi(\pi_{ng}|z_n)${{< /math >}}.

## 1. Recall the General Expectation for a Discrete Variable

If X is a discrete random variable that can take on values {{< math >}}$x_1, x_2, \ldots, x_K${{< /math >}} with corresponding probabilities {{< math >}}$P(X=x_1), P(X=x_2), \ldots, P(X=x_K)${{< /math >}}, then:

{{< math >}}
$$\mathbb{E}[X] = \sum_{i=1}^K x_i \cdot P(X=x_i) \tag{2}$$
{{< /math >}}

## 2. Applying to Our Velocity Case

**The Random Variable:** Our random variable is {{< math >}}$v^{(g)}(t(k_{ng})_{ng}, k_{ng})${{< /math >}}. Its value depends on the specific kinetic state {{< math >}}$k_{ng}${{< /math >}}.

**The Possible Values:** For each state {{< math >}}$k \in \{1,2,3,4\}${{< /math >}}, the value of the velocity is {{< math >}}$v^{(g)}(t_{ng}^{(k)}(z_n), k)${{< /math >}}. This velocity depends on the latent time {{< math >}}$t_{ng}^{(k)}(z_n)${{< /math >}} (which is a function of {{< math >}}$z_n${{< /math >}} and k) and the state k itself.

**The Probabilities:** The probability of being in state k is {{< math >}}$P(k_{ng} = k|z_n)${{< /math >}}. As derived previously, this is given by the expected value of the k-th component of {{< math >}}$\pi_{ng}${{< /math >}} under the Dirichlet distribution:

{{< math >}}
$$P(k_{ng} = k|z_n) = \mathbb{E}_{q_\phi(\pi_{ng}|z_n)}[\pi_{ng,k}] = \frac{\alpha_{q,ng,k}(z_n)}{\sum_{j=1}^K \alpha_{q,ng,j}(z_n)} \tag{3}$$
{{< /math >}}

Putting it together, the inner expectation for velocity becomes:

{{< math >}}
$$\mathbb{E}_{q_\phi(\pi_{ng}|z_n)} \left[v^{(g)}(t(k_{ng})_{ng}, k_{ng})\right] = \sum_{k=1}^K P(k_{ng} = k|z_n) \cdot v^{(g)}(t_{ng}^{(k)}(z_n), k) \tag{4}$$
{{< /math >}}

This will result in a single scalar value for each gene and cell, given a specific {{< math >}}$z_n${{< /math >}}.

## Numeric Example for Cell 1, Gene A

Let's use the same example inputs from our previous discussions:

**Assumed Sampled** {{< math >}}$z_1^*${{< /math >}}: {{< math >}}$(0.6558, -0.3779)${{< /math >}} (This comes from the first step in the "posterior predictive velocity distribution" process).

**Learned Biophysical Parameters (θ) for Gene A:**
- {{< math >}}$\alpha_A = 1.2${{< /math >}} (Transcription rate in Induction)
- {{< math >}}$\beta_A = 0.2${{< /math >}} (Splicing rate)
- {{< math >}}$\gamma_A = 0.1${{< /math >}} (Degradation rate)
- {{< math >}}$t_A^s = 8.0${{< /math >}} (Switching time)

**Encoder Output for** {{< math >}}$\pi_{1A}${{< /math >}} **(given** {{< math >}}$z_1^*${{< /math >}}**):**
- {{< math >}}$\alpha_{q,1A}(z_1^*) = (5.0, 0.1, 0.2, 0.1)${{< /math >}}
- Sum = 5.4

This gives us the state probabilities {{< math >}}$P(k_{1A} = k|z_1^*)${{< /math >}}:
- {{< math >}}$P(k_{1A} = 1|z_1^*) = 5.0/5.4 \approx 0.9259${{< /math >}}
- {{< math >}}$P(k_{1A} = 2|z_1^*) = 0.1/5.4 \approx 0.0185${{< /math >}}
- {{< math >}}$P(k_{1A} = 3|z_1^*) = 0.2/5.4 \approx 0.0370${{< /math >}}
- {{< math >}}$P(k_{1A} = 4|z_1^*) = 0.1/5.4 \approx 0.0185${{< /math >}}

**Decoder Output for** {{< math >}}$\rho_{1A}^{(k)}${{< /math >}} **(given** {{< math >}}$z_1^*${{< /math >}}**) and corresponding latent times** {{< math >}}$t_{1A}^{(k)}${{< /math >}}**:**
- {{< math >}}$\rho_{1A}^{(1)} = 0.5 \Rightarrow t_{1A}^{(1)} = 0.5 \times 8.0 = 4.0${{< /math >}}
- {{< math >}}$\rho_{1A}^{(2)} = 0.3 \Rightarrow t_{1A}^{(2)} = 0.3 \times 8.0 = 2.4${{< /math >}}
- {{< math >}}$\rho_{1A}^{(3)} = 0.8 \Rightarrow t_{1A}^{(3)} = 0.8 \times 8.0 = 6.4${{< /math >}}
- {{< math >}}$\rho_{1A}^{(4)} = 0.1 \Rightarrow t_{1A}^{(4)} = 0.1 \times 8.0 = 0.8${{< /math >}}

**The Velocity Formula:**

{{< math >}}
$$v^{(g)}(t^{(k)}, k) = \beta_g \bar{u}^{(g)}(t^{(k)}, k) - \gamma_g \bar{s}^{(g)}(t^{(k)}, k) \tag{5}$$
{{< /math >}}

We need the expected unspliced ({{< math >}}$\bar{u}${{< /math >}}) and spliced ({{< math >}}$\bar{s}${{< /math >}}) counts for each state at its corresponding latent time. These are given by the kinetic equations in VeloVI.

### Calculations for each state k and its corresponding {{< math >}}$t_{1A}^{(k)}${{< /math >}}:

**State 1 (Induction, k=1): t=4.0**

{{< math >}}
$$\bar{u}^{(A)}(4.0, 1) = \frac{\alpha_A}{\beta_A}(1 - e^{-\beta_A t}) = \frac{1.2}{0.2}(1 - e^{-0.2 \times 4.0}) = 6(1 - e^{-0.8}) \approx 6(1 - 0.4493) = 3.3042 \tag{6}$$
{{< /math >}}

{{< math >}}
$$\bar{s}^{(A)}(4.0, 1) = \frac{\alpha_A}{\gamma_A}\left(1 - \frac{\gamma_A}{\gamma_A - \beta_A}e^{-\beta_A t} + \frac{\beta_A}{\gamma_A - \beta_A}e^{-\gamma_A t}\right) \tag{7}$$
{{< /math >}}

{{< math >}}
$$= \frac{1.2}{0.1}\left(1 - \frac{0.1}{0.1 - 0.2}e^{-0.8} + \frac{0.2}{0.1 - 0.2}e^{-0.4}\right) = 12(1 + e^{-0.8} - 2e^{-0.4}) \tag{8}$$
{{< /math >}}

{{< math >}}
$$\approx 12(1 + 0.4493 - 2 \times 0.6703) = 12(1.4493 - 1.3406) \approx 1.3044 \tag{9}$$
{{< /math >}}

{{< math >}}
$$v^{(A)}(4.0, 1) = \beta_A \bar{u}^{(A)}(4.0, 1) - \gamma_A \bar{s}^{(A)}(4.0, 1) = (0.2 \times 3.3042) - (0.1 \times 1.3044) = 0.66084 - 0.13044 = 0.5304 \tag{10}$$
{{< /math >}}

**State 2 (Repression, k=2): t=2.4**

(For repression, we need the values at {{< math >}}$t_A^s = 8.0${{< /math >}} (end of induction) as initial conditions for decay.)

{{< math >}}
$$\bar{u}^{(A)}(t_A^s, 1) = 6(1 - e^{-0.2 \times 8.0}) \approx 6(1 - 0.2019) = 4.7886 \tag{11}$$
{{< /math >}}

{{< math >}}
$$\bar{s}^{(A)}(t_A^s, 1) = 12(1 + e^{-1.6} - 2e^{-0.8}) \approx 12(1 + 0.2019 - 2 \times 0.4493) = 3.6396 \tag{12}$$
{{< /math >}}

{{< math >}}
$$\bar{u}^{(A)}(2.4, 2) = \bar{u}^{(A)}(t_A^s, 1)e^{-\beta_A t_{1A}^{(2)}} = 4.7886 \times e^{-0.2 \times 2.4} \approx 4.7886 \times 0.6188 = 2.9634 \tag{13}$$
{{< /math >}}

{{< math >}}
$$\bar{s}^{(A)}(2.4, 2) = \bar{s}^{(A)}(t_A^s, 1)e^{-\gamma_A t_{1A}^{(2)}} + \frac{\beta_A \bar{u}^{(A)}(t_A^s, 1)}{\gamma_A - \beta_A}(e^{-\beta_A t_{1A}^{(2)}} - e^{-\gamma_A t_{1A}^{(2)}}) \tag{14}$$
{{< /math >}}

{{< math >}}
$$= 3.6396e^{-0.1 \times 2.4} + \frac{0.2 \times 4.7886}{0.1 - 0.2}(e^{-0.2 \times 2.4} - e^{-0.1 \times 2.4}) \tag{15}$$
{{< /math >}}

{{< math >}}
$$\approx 3.6396 \times 0.7866 + (-9.5772)(0.6188 - 0.7866) \approx 2.8617 + 1.6072 = 4.4689 \tag{16}$$
{{< /math >}}

{{< math >}}
$$v^{(A)}(2.4, 2) = (0.2 \times 2.9634) - (0.1 \times 4.4689) = 0.59268 - 0.44689 = 0.1458 \tag{17}$$
{{< /math >}}

**State 3 (Steady-State Induction, k=3): t=6.4**

{{< math >}}
$$\bar{u}^{(A)}(6.4, 3) = \alpha_A/\beta_A = 1.2/0.2 = 6.0 \tag{18}$$
{{< /math >}}

{{< math >}}
$$\bar{s}^{(A)}(6.4, 3) = \alpha_A/\gamma_A = 1.2/0.1 = 12.0 \tag{19}$$
{{< /math >}}

{{< math >}}
$$v^{(A)}(6.4, 3) = \beta_A \bar{u}^{(A)}(6.4, 3) - \gamma_A \bar{s}^{(A)}(6.4, 3) = (0.2 \times 6.0) - (0.1 \times 12.0) = 1.2 - 1.2 = 0.0 \tag{20}$$
{{< /math >}}

**State 4 (Steady-State Repression, k=4): t=0.8**

{{< math >}}
$$\bar{u}^{(A)}(0.8, 4) = 0.0 \tag{21}$$
{{< /math >}}

{{< math >}}
$$\bar{s}^{(A)}(0.8, 4) = 0.0 \tag{22}$$
{{< /math >}}

{{< math >}}
$$v^{(A)}(0.8, 4) = \beta_A \bar{u}^{(A)}(0.8, 4) - \gamma_A \bar{s}^{(A)}(0.8, 4) = (0.2 \times 0.0) - (0.1 \times 0.0) = 0.0 \tag{23}$$
{{< /math >}}

### Calculate the weighted sum:

Now, we multiply each calculated velocity by its corresponding state probability {{< math >}}$P(k_{1A} = k|z_1^*)${{< /math >}} and sum them up:

{{< math >}}
$$\mathbb{E}_{q_\phi(\pi_{1A}|z_1^*)} \left[v^{(A)}(t(k_{1A})_{1A}, k_{1A})\right] = \tag{24}$$
{{< /math >}}

{{< math >}}
$$(0.9259 \times 0.5304) + (0.0185 \times 0.1458) + (0.0370 \times 0.0) + (0.0185 \times 0.0) \tag{25}$$
{{< /math >}}

{{< math >}}
$$= 0.4900 + 0.0027 + 0.0 + 0.0 \approx 0.4927 \tag{26}$$
{{< /math >}}

## Get Posterior Predictive Velocity Distribution

The value **0.4927** is the expected velocity for Gene A in Cell 1, given the specific sampled latent variable {{< math >}}$z_1^* = (0.6558, -0.3779)${{< /math >}}.

This value is one sample of the posterior predictive velocity. To get the full posterior predictive velocity distribution, you would repeat this entire calculation for many different samples of {{< math >}}$z_n^*${{< /math >}} (drawn from {{< math >}}$q_\phi(z_n|u_n, s_n)${{< /math >}}) and then analyze the distribution of these resulting expected velocity values.
