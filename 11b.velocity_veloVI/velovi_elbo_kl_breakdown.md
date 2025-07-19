---
title: "Velovi Elbo Kl Breakdown"
date: 2025-07-10
draft: True
---

# Breaking KL Terms into Manageable, Conditional Pieces

**break the KL terms into manageable, conditional pieces** in the context of your ELBO--

{{< math >}} $$\mathcal{L}_{\text{elbo}}(\theta, \phi; u, s) = \sum_n - \mathbb{E}_{q_\phi(z_n, \pi_n | u_n, s_n)} [\log p_\theta(u_n, s_n | z_n, \pi_n)]$$ {{< /math >}}
{{< math >}} $$+ \text{KL} (q_\phi(z_n | u_n, s_n) \| p(z))$$ {{< /math >}}
{{< math >}} $$+ \mathbb{E}_{q_\phi(z_n | u_n, s_n)} \left[\sum_g \text{KL} (q_\phi(\pi_{ng} | z_n) \| p(\pi_{ng}))\right], \tag{26}$$ {{< /math >}}

## 🔁 Starting Point: What's the KL term in a basic ELBO?

In standard variational inference, we approximate the true posterior {{< math >}} $p(z|x)$ {{< /math >}} with a simpler distribution {{< math >}} $q(z|x)$ {{< /math >}}. The **Evidence Lower Bound (ELBO)** becomes:

{{< math >}} $$ \text{ELBO} = \mathbb{E}_{q(z)}[\log p(x | z)] - \text{KL}(q(z) \| p(z)) $$ {{< /math >}}

This KL divergence:
{{< math >}} $$ \text{KL}(q(z) \| p(z)) $$ {{< /math >}}

measures how far the approximate posterior {{< math >}} $q(z)$ {{< /math >}} is from the true prior {{< math >}} $p(z)$ {{< /math >}}. That's the **standard, single-latent variable case**.

## 🌿 Now Enter: Hierarchical Latents (like z, π)

In your ELBO, the model has *multiple* latent variables:
* {{< math >}} $z_n$ {{< /math >}} (global, per-cell latent)
* {{< math >}} $\pi_{ng}$ {{< /math >}} (local, per-gene latent for each cell)

We now have an **approximate joint posterior**:
{{< math >}} $$ q(z_n, \pi_n | u_n, s_n) $$ {{< /math >}}

This could be intractable or too complex. So we factor it to **make it more manageable**:
{{< math >}} $$ q(z_n, \pi_n | u_n, s_n) = q(z_n | u_n, s_n) \prod_g q(\pi_{ng} | z_n) $$ {{< /math >}}

This means:
* First infer {{< math >}} $z_n$ {{< /math >}} from data.
* Then, infer each {{< math >}} $\pi_{ng}$ {{< /math >}} from {{< math >}} $z_n$ {{< /math >}}, **not directly from the full data**.

This makes the math simpler and lets us **break the big KL divergence**:
{{< math >}} $$ \text{KL}(q(z_n, \pi_n | u_n, s_n) \| p(z_n, \pi_n)) $$ {{< /math >}}

into **two smaller, manageable pieces**:

1. **KL between {{< math >}} $q(z_n)$ {{< /math >}} and prior {{< math >}} $p(z_n)$ {{< /math >}}:**
   {{< math >}} $$ \text{KL}(q(z_n | u_n, s_n) \| p(z_n)) $$ {{< /math >}}

2. **KLs between each gene-level posterior and its prior, conditioned on {{< math >}} $z_n$ {{< /math >}}:**
   {{< math >}} $$ \mathbb{E}_{q(z_n)} \left[ \sum_g \text{KL}(q(\pi_{ng} | z_n) \| p(\pi_{ng})) \right] $$ {{< /math >}}

These KL terms are easier to compute, optimize, and interpret.