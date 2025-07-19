---
title: "Unraveled Eq15 16"
date: 2025-07-10
draft: True
---

# Stochastic Gene Expression Analysis

## 🔍 The Expression

The time-averaged probability generating function is given by:

{{< math >}} 
$$E_f[G(z,t)] = \sum_{x=0}^{\infty} \int P(x,t)z^x \, df = \frac{1}{T}\int_0^T \sum_{x=0}^{\infty} P(x,t)z^x \, dt = \frac{1}{T}\int_0^T G(z,t) \, dt \quad (12)$$
{{< /math >}}

Let's break it down.

## 🔣 What Does Each Term Mean?

| Term | Meaning |
|------|---------|
| {{< math >}} $G(z,t) = \sum_{x=0}^{\infty} P(x,t)z^x$ {{< /math >}} | The probability generating function at time {{< math >}} $t$ {{< /math >}} |
| {{< math >}} $df = \frac{1}{T}dt$ {{< /math >}} | Assumes cells are sampled uniformly across time interval {{< math >}} $[0,T]$ {{< /math >}} — desynchronization assumption |
| {{< math >}} $E_f[G(z,t)]$ {{< /math >}} | Expectation of the generating function under the time-uniform measure {{< math >}} $f$ {{< /math >}} |

## 🧠 Interpretation

### 1. Desynchronized Observations

In scRNA-seq, each cell is assumed to be a random snapshot of a stochastic process — but not all at the same time. So instead of analyzing {{< math >}} $G(z,t)$ {{< /math >}} at a fixed {{< math >}} $t$ {{< /math >}}, we average it over the full observable process window {{< math >}} $[0,T]$ {{< /math >}}.

This leads to integrating over time using a uniform distribution: {{< math >}} $df = \frac{1}{T}dt$ {{< /math >}}, meaning all time points in {{< math >}} $[0,T]$ {{< /math >}} are equally likely.

### 2. Swapping Sum and Integral

The first expression:

{{< math >}} 
$$\sum_{x=0}^{\infty} \int P(x,t)z^x \, df$$
{{< /math >}}

is re-written as:

{{< math >}} 
$$\frac{1}{T}\int_0^T \sum_{x=0}^{\infty} P(x,t)z^x \, dt$$
{{< /math >}}

Here, you're swapping the sum and integral, which is allowed under standard convergence conditions (Fubini's theorem), because all terms are non-negative and well-behaved.

### 3. Recovering the Time-Averaged PGF

The inner sum is just the definition of {{< math >}} $G(z,t)$ {{< /math >}}, so the final expression becomes:

{{< math >}} 
$$E_f[G(z,t)] = \frac{1}{T}\int_0^T G(z,t) \, dt \quad (13)$$
{{< /math >}}

This gives you the expected generating function across time, assuming you observe cells uniformly sampled over the interval {{< math >}} $[0,T]$ {{< /math >}}.

## 📌 Summary

This derivation is important in modeling time-averaged behavior of stochastic systems like gene expression, especially in snapshot data such as scRNA-seq, where you don't observe a full trajectory, but a desynchronized ensemble.

So this final quantity:

{{< math >}} 
$$E_f[G(z,t)] = \frac{1}{T}\int_0^T G(z,t) \, dt$$
{{< /math >}}

gives a time-averaged probability generating function, which can then be used to compute moments (mean, variance) or to infer the steady-state distribution by taking the limit {{< math >}} $T \to \infty$ {{< /math >}}.


# Equation 16 in ["RNA velocity unraveled"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010492)
## 🔹 1. The Mean (Expected Value) Over Time

{{< math >}} 
$$y(t) = \left(y_0 - \frac{\alpha}{\beta}\right)e^{-\beta t} + \frac{\alpha}{\beta} \equiv ce^{-\beta t} + \frac{\alpha}{\beta} \quad (14)$$
{{< /math >}}

**Interpretation:**
- {{< math >}} $y(t) = E[X_t]$ {{< /math >}}: the expected number of transcripts at time {{< math >}} $t$ {{< /math >}}
- {{< math >}} $\alpha$ {{< /math >}}: production (birth) rate
- {{< math >}} $\beta$ {{< /math >}}: degradation (death) rate
- {{< math >}} $y_0$ {{< /math >}}: initial expected number of molecules at {{< math >}} $t = 0$ {{< /math >}}
- {{< math >}} $c = y_0 - \frac{\alpha}{\beta}$ {{< /math >}}: shorthand for compact expression

This is the solution to the ODE:

{{< math >}} 
$$\frac{dy}{dt} = \alpha - \beta y \quad (15)$$
{{< /math >}}

which models exponential approach to the steady-state mean {{< math >}} $\alpha/\beta$ {{< /math >}}. It's a standard result from linear first-order ODEs.

## 🔹 2. The Probability Generating Function (PGF)

{{< math >}} 
$$G(z,t) = E[z^{X_t}] = e^{(z-1)y(t)} \equiv e^{uy(t)} \quad (16)$$
{{< /math >}}

**Interpretation:**

The PGF is a compact way to encode the full distribution {{< math >}} $P(X_t = x)$ {{< /math >}} in a generating function.

{{< math >}} $u = z - 1$ {{< /math >}}: just a substitution to simplify notation.

This specific form,

{{< math >}} 
$$G(z,t) = e^{(z-1)y(t)},$$
{{< /math >}}

tells us the distribution of {{< math >}} $X_t$ {{< /math >}} is Poisson-distributed with mean {{< math >}} $y(t)$ {{< /math >}}, since:

The PGF of a Poisson random variable with mean {{< math >}} $\lambda$ {{< /math >}} is:

{{< math >}} 
$$G(z) = e^{(z-1)\lambda} \quad (17)$$
{{< /math >}}

So:

{{< math >}} 
$$X_t \sim \text{Poisson}(y(t))$$
{{< /math >}}

This is a key classical result for the linear birth-death process with constant rates and Poisson initial condition — the distribution remains Poisson at all {{< math >}} $t$ {{< /math >}}, with mean evolving over time via {{< math >}} $y(t)$ {{< /math >}}.

### Derivation

Suppose {{< math >}} $P(x,t)$ {{< /math >}} is Poisson-distributed with time-dependent mean {{< math >}} $y(t)$ {{< /math >}}:

{{< math >}} 
$$P(x,t) = \frac{y(t)^x e^{-y(t)}}{x!} \quad (18)$$
{{< /math >}}

Then the PGF is:

{{< math >}} 
$$\begin{align}
G(z,t) &= \sum_{x=0}^{\infty} \frac{y(t)^x e^{-y(t)}}{x!} z^x \\
&= e^{-y(t)} \sum_{x=0}^{\infty} \frac{(y(t)z)^x}{x!} \\
&= e^{-y(t)} e^{y(t)z} \\
&= e^{y(t)(z-1)} \quad (19)
\end{align}$$
{{< /math >}}

## 📌 Interim Summary

- The expression for {{< math >}} $y(t)$ {{< /math >}} comes from solving a first-order linear ODE describing mean dynamics of the birth-death process.

- The generating function {{< math >}} $G(z,t) = e^{(z-1)y(t)}$ {{< /math >}} reveals the full time-dependent distribution {{< math >}} $X_t \sim \text{Poisson}(y(t))$ {{< /math >}}.

- At {{< math >}} $t \to \infty$ {{< /math >}}, this converges to {{< math >}} $\text{Poisson}(\alpha/\beta)$ {{< /math >}}, the steady-state.

## Time-averaged Generating Function

You then average the PGF over time to account for the fact that in single-cell RNA-seq, cells are observed asynchronously across the time interval {{< math >}} $[0,T]$ {{< /math >}}. That is:

{{< math >}}
$$H(z) = E_f[G(z,t)] = \frac{1}{T}\int_0^T G(z,t)\,dt = \frac{1}{T}\int_0^T e^{uy(t)}dt \quad (12)$$
{{< /math >}}

Now substitute the form of {{< math >}} $y(t) = ce^{-\beta t} + \frac{\alpha}{\beta}$ {{< /math >}}:

{{< math >}}
$$= \frac{1}{T}\int_0^T \exp\left(u\left(ce^{-\beta t} + \frac{\alpha}{\beta}\right)\right)dt = \frac{1}{T}e^{u\alpha/\beta}\int_0^T \exp(uce^{-\beta t})dt \quad (13)$$
{{< /math >}}

## Change of Variables & Exponential Integral

Make substitution:
Let {{< math >}} $s = e^{-\beta t} \Rightarrow t = -\frac{1}{\beta}\ln s \Rightarrow dt = -\frac{1}{\beta s}ds$ {{< /math >}}

As {{< math >}} $t$ {{< /math >}} goes from {{< math >}} $0 \to T$ {{< /math >}}, {{< math >}} $s$ {{< /math >}} goes from {{< math >}} $1 \to e^{-\beta T}$ {{< /math >}}

So:

{{< math >}}
$$\int_0^T \exp(uce^{-\beta t})dt = \frac{1}{\beta}\int_{e^{-\beta T}}^1 \frac{1}{s}\exp(ucs)ds \quad (14)$$
{{< /math >}}

This is a standard exponential integral:

{{< math >}}
$$\int \frac{e^{as}}{s}ds = \text{Ei}(as) \quad (15)$$
{{< /math >}}

denoted:

{{< math >}}
$$\text{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t}\,dt \quad \text{(for real } x < 0\text{)} \quad (16)$$
{{< /math >}}

or equivalently (for {{< math >}} $x > 0$ {{< /math >}}):

{{< math >}}
$$\text{Ei}(x) = -\int_{-x}^\infty \frac{e^{-t}}{t}\,dt \quad (17)$$
{{< /math >}}

Thus:

{{< math >}}
$$H(z) = \frac{e^{u\alpha/\beta}}{T\beta}[\text{Ei}(uc) - \text{Ei}(uce^{-\beta T})] \quad (18)$$
{{< /math >}}

Or compactly:

{{< math >}}
$$H(z) = \frac{e^{u\alpha/\beta}}{\beta T}(\text{Ei}(uc) - \text{Ei}(uce^{-\beta T})) \quad (19)$$
{{< /math >}}

## ✅ Interpretation

{{< math >}} $H(z)$ {{< /math >}} is the time-averaged generating function of the count distribution observed in single-cell RNA-seq experiments.

It accounts for **desynchronization of cells** — i.e., cells are randomly captured along the time axis, not in lockstep.

The resulting distribution is **not strictly Poisson**, but a time-average of Poisson distributions — which can exhibit **overdispersion**, consistent with what's observed in real data.

The appearance of the **exponential integral Ei** is due to integrating the exponential of an exponential — common in birth-death analyses.