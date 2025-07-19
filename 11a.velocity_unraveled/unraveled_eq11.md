---
title: "Unraveled Eq11"
date: 2025-07-10
draft: True
---

# A) Proof: A Linear Birth-Death Process Preserves Poisson Distribution

We aim to prove that if a system characterized by constant production (α) and first-order degradation (β) starts with a Poisson-distributed number of molecules, it will remain Poisson-distributed for all subsequent times, and its mean will evolve according to a simple deterministic ODE.

## 1. Framework and Definitions

We consider a single molecular species X with its copy number denoted by x. The reactions are:

- **Production**: {{< math >}} $0 \xrightarrow{\alpha} X$ {{< /math >}} (molecules are created at a constant rate α, independent of current x)
- **Degradation**: {{< math >}} $X \xrightarrow{\beta} 0$ {{< /math >}} (molecules are destroyed at a rate proportional to their current number, βx)

The probability P(x;t) of having x molecules at time t is governed by the Chemical Master Equation (CME):

{{< math >}} 
$$ \frac{dP(x;t)}{dt} = \alpha[P(x-1;t) - P(x;t)] + \beta[(x+1)P(x+1;t) - xP(x;t)] \tag{1} $$
{{< /math >}}

**Interpretation of (1):**

The first term {{< math >}} $\alpha[P(x-1;t) - P(x;t)]$ {{< /math >}} accounts for changes due to production:
- {{< math >}} $\alpha P(x-1;t)$ {{< /math >}}: Probability flux into state x from state x−1 (a molecule is produced).
- {{< math >}} $-\alpha P(x;t)$ {{< /math >}}: Probability flux out of state x to state x+1 (a molecule is produced from state x).

The second term {{< math >}} $\beta[(x+1)P(x+1;t) - xP(x;t)]$ {{< /math >}} accounts for changes due to degradation:
- {{< math >}} $\beta(x+1)P(x+1;t)$ {{< /math >}}: Probability flux into state x from state x+1 (a molecule degrades from state x+1). The (x+1) factor represents the rate at which degradation occurs when there are x+1 molecules.
- {{< math >}} $-\beta xP(x;t)$ {{< /math >}}: Probability flux out of state x to state x−1 (a molecule degrades from state x). The x factor represents the rate at which degradation occurs when there are x molecules.

**Boundary Conditions:** For x=0, {{< math >}} $P(x-1;t)$ {{< /math >}} and {{< math >}} $(x+1)P(x+1;t)$ {{< /math >}} terms need careful handling (e.g., {{< math >}} $P(-1,t)=0$ {{< /math >}}).

## 2. Poisson Ansatz

We propose a solution of the Poisson form for the probability distribution P(x;t):

{{< math >}} 
$$ P(x;t) = \frac{y(t)^x e^{-y(t)}}{x!} \tag{2} $$
{{< /math >}}

where y(t) is the time-dependent mean of the distribution ({{< math >}} $y(t) = \langle x \rangle_t$ {{< /math >}}).

## 3. Moment Analysis (Key Innovation)

We analyze how the moments of the distribution evolve.

### 3.1 First Moment Equation (Mean Evolution)

We multiply the CME (1) by x and sum over all possible values of x (from 0 to ∞) to derive the differential equation for the mean {{< math >}} $\langle x \rangle_t = \sum_{x=0}^{\infty} x P(x;t)$ {{< /math >}}.

{{< math >}} 
$$ \frac{d\langle x \rangle}{dt} = \sum_{x=0}^{\infty} x \frac{dP(x;t)}{dt} $$
{{< /math >}}

{{< math >}} 
$$ = \sum_{x=0}^{\infty} x[\alpha P(x-1;t) - \alpha P(x;t) + \beta(x+1)P(x+1;t) - \beta x P(x;t)] $$
{{< /math >}}

We break this into four summations:

{{< math >}} 
$$ \frac{d\langle x \rangle}{dt} = \alpha \sum_{x=0}^{\infty} x P(x-1;t) - \alpha \sum_{x=0}^{\infty} x P(x;t) + \beta \sum_{x=0}^{\infty} x(x+1)P(x+1;t) - \beta \sum_{x=0}^{\infty} x^2 P(x;t) $$
{{< /math >}}

Let's simplify each sum:

**First term**: {{< math >}} $\alpha \sum_{x=0}^{\infty} x P(x-1;t)$ {{< /math >}}. Let {{< math >}} $k = x-1$ {{< /math >}}, so {{< math >}} $x = k+1$ {{< /math >}}. When {{< math >}} $x = 0$ {{< /math >}}, {{< math >}} $k = -1$ {{< /math >}}, but {{< math >}} $P(-1;t) = 0$ {{< /math >}}. So the sum effectively starts from {{< math >}} $k = 0$ {{< /math >}}.

{{< math >}} 
$$ \alpha \sum_{k=0}^{\infty} (k+1)P(k;t) = \alpha\left(\sum_{k=0}^{\infty} k P(k;t) + \sum_{k=0}^{\infty} P(k;t)\right) = \alpha(\langle x \rangle + 1) $$
{{< /math >}}

**Second term**: {{< math >}} $-\alpha \sum_{x=0}^{\infty} x P(x;t) = -\alpha \langle x \rangle$ {{< /math >}}.

**Third term**: {{< math >}} $\beta \sum_{x=0}^{\infty} x(x+1)P(x+1;t)$ {{< /math >}}. Let {{< math >}} $k = x+1$ {{< /math >}}, so {{< math >}} $x = k-1$ {{< /math >}}. When {{< math >}} $x = 0$ {{< /math >}}, {{< math >}} $k = 1$ {{< /math >}}.

{{< math >}} 
$$ \beta \sum_{k=1}^{\infty} (k-1)k P(k;t) = \beta\left(\sum_{k=1}^{\infty} k^2 P(k;t) - \sum_{k=1}^{\infty} k P(k;t)\right) = \beta(\langle x^2 \rangle - \langle x \rangle) $$
{{< /math >}}

**Fourth term**: {{< math >}} $-\beta \sum_{x=0}^{\infty} x^2 P(x;t) = -\beta \langle x^2 \rangle$ {{< /math >}}.

Combining these:

{{< math >}} 
$$ \frac{d\langle x \rangle}{dt} = \alpha(\langle x \rangle + 1) - \alpha \langle x \rangle + \beta(\langle x^2 \rangle - \langle x \rangle) - \beta \langle x^2 \rangle $$
{{< /math >}}

{{< math >}} 
$$ \frac{d\langle x \rangle}{dt} = \alpha + \alpha \langle x \rangle - \alpha \langle x \rangle + \beta \langle x^2 \rangle - \beta \langle x \rangle - \beta \langle x^2 \rangle $$
{{< /math >}}

{{< math >}} 
$$ \frac{d\langle x \rangle}{dt} = \alpha - \beta \langle x \rangle $$
{{< /math >}}

Since we defined {{< math >}} $y(t) = \langle x \rangle$ {{< /math >}}, we get the mean evolution equation:

{{< math >}} 
$$ \frac{dy}{dt} = \alpha - \beta y(t) \tag{3} $$
{{< /math >}}

This matches the deterministic rate equation for the average copy number.

### 3.2 Second Moment Verification

The derivation for the second moment from the CME is:

{{< math >}} 
$$ \frac{d\langle x^2 \rangle}{dt} = \alpha(1 + 2\langle x \rangle) - \beta(2\langle x^2 \rangle - \langle x \rangle) $$
{{< /math >}}

For a Poisson distribution, the variance is equal to the mean ({{< math >}} $\sigma^2 = y$ {{< /math >}}). The second uncentered moment is given by {{< math >}} $\langle x^2 \rangle = \text{Var}[x] + \langle x \rangle^2 = y(t) + y(t)^2$ {{< /math >}}.

If we substitute {{< math >}} $\langle x \rangle = y$ {{< /math >}} and {{< math >}} $\langle x^2 \rangle = y^2 + y$ {{< /math >}} into this equation, and {{< math >}} $\frac{d\langle x^2 \rangle}{dt} = \frac{d}{dt}(y^2 + y) = 2y \frac{dy}{dt} + \frac{dy}{dt} = (2y+1) \frac{dy}{dt}$ {{< /math >}}:

{{< math >}} 
$$ (2y+1) \frac{dy}{dt} = \alpha(1+2y) - \beta(2(y^2 + y) - y) $$
{{< /math >}}

{{< math >}} 
$$ (2y+1)(\alpha - \beta y) = \alpha(1+2y) - \beta(2y^2 + 2y - y) $$
{{< /math >}}

{{< math >}} 
$$ (2y+1)(\alpha - \beta y) = \alpha(1+2y) - \beta(2y^2 + y) $$
{{< /math >}}

{{< math >}} 
$$ \alpha + 2\alpha y - \beta y - 2\beta y^2 = \alpha + 2\alpha y - 2\beta y^2 - \beta y $$
{{< /math >}}

The two sides are equal, confirming consistency. The key is that the Poisson ansatz self-consistently satisfies the higher moment equations when the mean follows (3).

## 4. Consistency Check

This step directly verifies that if P(x,t) is Poisson with mean y(t) that satisfies {{< math >}} $\frac{dy}{dt} = \alpha - \beta y$ {{< /math >}}, then P(x,t) is indeed a solution to the CME.

**Time derivative of Poisson PMF (LHS of CME):**

{{< math >}} 
$$ \frac{dP(x;t)}{dt} = \frac{d}{dt}\left(\frac{y(t)^x e^{-y(t)}}{x!}\right) $$
{{< /math >}}

Using product rule and chain rule (recall {{< math >}} $\frac{d}{dt} y(t)^x = x y(t)^{x-1} \frac{dy}{dt}$ {{< /math >}} and {{< math >}} $\frac{d}{dt} e^{-y(t)} = -e^{-y(t)} \frac{dy}{dt}$ {{< /math >}}):

{{< math >}} 
$$ \frac{dP(x;t)}{dt} = \frac{1}{x!} \left(x y^{x-1} e^{-y} \frac{dy}{dt} - y^x e^{-y} \frac{dy}{dt}\right) $$
{{< /math >}}

{{< math >}} 
$$ \frac{dP(x;t)}{dt} = \frac{1}{x!} y^x e^{-y} \left(\frac{x}{y} \frac{dy}{dt} - \frac{dy}{dt}\right) $$
{{< /math >}}

{{< math >}} 
$$ \frac{dP(x;t)}{dt} = P(x;t)\left(\frac{x}{y} - 1\right) \frac{dy}{dt} \tag{5} $$
{{< /math >}}

**Substitute Poisson PMF into the CME (RHS of CME):**

{{< math >}} 
$$ \alpha\left[\frac{y^{x-1} e^{-y}}{(x-1)!} - \frac{y^x e^{-y}}{x!}\right] + \beta\left[(x+1) \frac{y^{x+1} e^{-y}}{(x+1)!} - x \frac{y^x e^{-y}}{x!}\right] $$
{{< /math >}}

{{< math >}} 
$$ = \alpha \frac{y^x e^{-y}}{x!} \left[\frac{x}{y} - 1\right] + \beta \frac{y^x e^{-y}}{x!} \left[y - x\right] $$
{{< /math >}}

{{< math >}} 
$$ = P(x;t)\left[\alpha\left(\frac{x}{y} - 1\right) + \beta(y - x)\right] \tag{4} $$
{{< /math >}}

**Equate (4) and (5):**

{{< math >}} 
$$ P(x;t)\left(\frac{x}{y} - 1\right) \frac{dy}{dt} = P(x;t)\left[\alpha\left(\frac{x}{y} - 1\right) + \beta(y - x)\right] $$
{{< /math >}}

Divide by {{< math >}} $P(x;t)$ {{< /math >}} (assuming {{< math >}} $P(x;t) \neq 0$ {{< /math >}}):

{{< math >}} 
$$ \left(\frac{x}{y} - 1\right) \frac{dy}{dt} = \alpha\left(\frac{x}{y} - 1\right) + \beta(y - x) $$
{{< /math >}}

{{< math >}} 
$$ \left(\frac{x-y}{y}\right) \frac{dy}{dt} = \alpha\left(\frac{x-y}{y}\right) - \beta(x-y) $$
{{< /math >}}

Multiply by y:

{{< math >}} 
$$ (x-y) \frac{dy}{dt} = \alpha(x-y) - \beta y(x-y) $$
{{< /math >}}

{{< math >}} 
$$ (x-y) \frac{dy}{dt} = (x-y)(\alpha - \beta y) $$
{{< /math >}}

For this to hold for all x (not just for {{< math >}} $x = y$ {{< /math >}}):

{{< math >}} 
$$ \frac{dy}{dt} = \alpha - \beta y $$
{{< /math >}}

This is consistent.

## 5. Solution and Verification

The ODE {{< math >}} $\frac{dy}{dt} = \alpha - \beta y$ {{< /math >}} is a first-order linear differential equation. Its solution with initial condition {{< math >}} $y(0) = y_0$ {{< /math >}} is:

{{< math >}} 
$$ y(t) = \frac{\alpha}{\beta} + \left(y_0 - \frac{\alpha}{\beta}\right)e^{-\beta t} \tag{6} $$
{{< /math >}}

**Steady-state:** As {{< math >}} $t \to \infty$ {{< /math >}}, {{< math >}} $e^{-\beta t} \to 0$ {{< /math >}}, so {{< math >}} $y(t) \to \frac{\alpha}{\beta}$ {{< /math >}}. This is the steady-state mean.

**Verification:**
- **Mean:** Confirmed, as it's the solution to (3).
- **Variance:** For a Poisson distribution, {{< math >}} $\sigma^2(t) = y(t)$ {{< /math >}}, which is consistent.
- **Higher moments:** All cumulants {{< math >}} $\kappa_n(t)$ {{< /math >}} of a Poisson distribution are equal to its mean y(t). The fact that these equations are consistently satisfied by y(t) reinforces that the Poisson distribution is maintained.

## 6. Conclusion

The proof establishes that the system (constant birth, linear death) preserves its Poisson character over time. This is because:

1. The mean evolution equation derived from the CME matches the deterministic rate equation (3).
2. The Poisson ansatz precisely satisfies the CME (Equation 1) when its mean follows Equation (3).
3. All higher moments automatically remain consistent with the Poisson statistics as dictated by y(t).

**Final Solution:**

The probability distribution for x at time t is Poisson with mean y(t):

{{< math >}} 
$$ P(x;t) = \frac{y(t)^x e^{-y(t)}}{x!} $$
{{< /math >}}

where y(t) is given by the solution to the ODE:

{{< math >}} 
$$ y(t) = \frac{\alpha}{\beta}(1 - e^{-\beta t}) + y_0 e^{-\beta t} $$
{{< /math >}}

**Physical Interpretation:**

- **Production (α) and degradation (β):** These two processes compete to determine the average number of molecules y(t) over time.

- **Poisson fluctuations preserved:** This crucial property arises because the production is a simple zeroth-order process (like adding a constant trickle of water to a bucket), and degradation is a first-order process (each molecule has an equal, independent chance of disappearing, like removing water from the bucket at a rate proportional to how much is there). These linear kinetics do not introduce additional complexity (like "bursting" where production itself is stochastic) that would lead to super-Poissonian noise.

# B) Explain Second Uncentered Moment of Poisson Distribution

This is a standard property of the Poisson distribution.

## Definition of Expected Value for a Discrete Random Variable

For a discrete random variable {{< math >}}$X${{< /math >}} with probability mass function {{< math >}}$P(x)${{< /math >}}, the expected value of {{< math >}}$g(X)${{< /math >}} is:

{{< math >}}
$$E[g(X)]=\sum_{x} g(x)P(x) \quad (1)$$
{{< /math >}}

## Definition of Poisson Distribution

A random variable {{< math >}}$X${{< /math >}} follows a Poisson distribution with parameter {{< math >}}$y${{< /math >}} (which is its mean), if its probability mass function is:

{{< math >}}
$$P(x)= \frac{y^x e^{-y}}{x!} \text{ for } x=0,1,2,\ldots \quad (2)$$
{{< /math >}}

## Goal

Prove {{< math >}}$\langle x^2 \rangle = E[X^2] = y^2 + y${{< /math >}}

## Proof Strategy

We will use the definition of {{< math >}}$E[X^2]${{< /math >}} and substitute the Poisson PMF. We can also use the relationship between variance and second moment: {{< math >}}$\text{Var}[X] = E[X^2] - (E[X])^2${{< /math >}}. If we know the variance and mean of a Poisson distribution, we can derive {{< math >}}$E[X^2]${{< /math >}}.

## Method 1: Using the Relationship between Variance and Mean

### Recall the mean of a Poisson distribution:

{{< math >}}
$$E[X] = \sum_{x=0}^{\infty} x \frac{y^x e^{-y}}{x!} \quad (3)$$
{{< /math >}}

{{< math >}}
$$E[X] = e^{-y} \sum_{x=1}^{\infty} \frac{xy^x}{x!} = e^{-y} \sum_{x=1}^{\infty} \frac{y^x}{(x-1)!} \text{ (since } \frac{x}{x!} = \frac{1}{(x-1)!}\text{)} \quad (4)$$
{{< /math >}}

Let {{< math >}}$k = x - 1${{< /math >}}. When {{< math >}}$x = 1${{< /math >}}, {{< math >}}$k = 0${{< /math >}}.

{{< math >}}
$$E[X] = e^{-y} \sum_{k=0}^{\infty} \frac{y^{k+1}}{k!} = e^{-y} y \sum_{k=0}^{\infty} \frac{y^k}{k!} \quad (5)$$
{{< /math >}}

Since {{< math >}}$\sum_{k=0}^{\infty} \frac{y^k}{k!} = e^y${{< /math >}} (Taylor series for {{< math >}}$e^y${{< /math >}}),

{{< math >}}
$$E[X] = e^{-y} ye^y = y \quad (6)$$
{{< /math >}}

So, the mean {{< math >}}$\langle x \rangle = y${{< /math >}}.

### Recall the variance of a Poisson distribution:

A fundamental property of the Poisson distribution is that its variance is equal to its mean.

{{< math >}}
$$\text{Var}[X] = y \quad (7)$$
{{< /math >}}

### Relate variance to the second moment:

The definition of variance is {{< math >}}$\text{Var}[X] = E[X^2] - (E[X])^2${{< /math >}}.

Rearranging to solve for {{< math >}}$E[X^2]${{< /math >}}:

{{< math >}}
$$E[X^2] = \text{Var}[X] + (E[X])^2 \quad (8)$$
{{< /math >}}

### Substitute the mean and variance of the Poisson distribution:

{{< math >}}
$$E[X^2] = y + (y)^2 \quad (9)$$
{{< /math >}}

{{< math >}}
$$E[X^2] = y^2 + y \quad (10)$$
{{< /math >}}

Therefore, {{< math >}}$\langle x^2 \rangle = y^2 + y${{< /math >}}.

## Method 2: Direct Calculation of {{< math >}}$E[X^2]${{< /math >}}

This method is more involved but directly calculates the sum.

{{< math >}}
$$E[X^2] = \sum_{x=0}^{\infty} x^2 P(x) = \sum_{x=0}^{\infty} x^2 \frac{y^x e^{-y}}{x!} \quad (11)$$
{{< /math >}}

We can rewrite {{< math >}}$x^2${{< /math >}} as {{< math >}}$x(x-1) + x${{< /math >}}. This trick is common for moments of discrete distributions, as it helps simplify factorials.

{{< math >}}
$$E[X^2] = \sum_{x=0}^{\infty} (x(x-1) + x) \frac{y^x e^{-y}}{x!} \quad (12)$$
{{< /math >}}

{{< math >}}
$$E[X^2] = \sum_{x=0}^{\infty} x(x-1) \frac{y^x e^{-y}}{x!} + \sum_{x=0}^{\infty} x \frac{y^x e^{-y}}{x!} \quad (13)$$
{{< /math >}}

Let's evaluate each sum:

### First Sum: {{< math >}}$\sum_{x=0}^{\infty} x(x-1) \frac{y^x e^{-y}}{x!}${{< /math >}}

For {{< math >}}$x = 0${{< /math >}} and {{< math >}}$x = 1${{< /math >}}, the term {{< math >}}$x(x-1)${{< /math >}} is zero, so we can start the sum from {{< math >}}$x = 2${{< /math >}}.

For {{< math >}}$x \geq 2${{< /math >}}, {{< math >}}$\frac{x(x-1)}{x!} = \frac{x(x-1)}{x(x-1)(x-2)!} = \frac{1}{(x-2)!}${{< /math >}}.

Thus, the first sum becomes:

{{< math >}}
$$e^{-y} \sum_{x=2}^{\infty} \frac{y^x}{(x-2)!} \quad (14)$$
{{< /math >}}

Let {{< math >}}$k = x - 2${{< /math >}}. When {{< math >}}$x = 2${{< /math >}}, {{< math >}}$k = 0${{< /math >}}.

{{< math >}}
$$e^{-y} \sum_{k=0}^{\infty} \frac{y^{k+2}}{k!} = e^{-y} y^2 \sum_{k=0}^{\infty} \frac{y^k}{k!} \quad (15)$$
{{< /math >}}

Since {{< math >}}$\sum_{k=0}^{\infty} \frac{y^k}{k!} = e^y${{< /math >}}:

{{< math >}}
$$e^{-y} y^2 e^y = y^2 \quad (16)$$
{{< /math >}}

### Second Sum: {{< math >}}$\sum_{x=0}^{\infty} x \frac{y^x e^{-y}}{x!}${{< /math >}}

This is simply the definition of {{< math >}}$E[X]${{< /math >}}. As shown in Method 1, this sum evaluates to {{< math >}}$y${{< /math >}}.

### Combining the sums:

{{< math >}}
$$E[X^2] = y^2 + y \quad (17)$$
{{< /math >}}

Therefore, for a Poisson distribution with mean {{< math >}}$y${{< /math >}}, {{< math >}}$\langle x^2 \rangle = y^2 + y${{< /math >}}.

This derivation is fundamental to understanding the baseline noise in gene expression and why more advanced models (like scVelo's, which handles bursting) are needed to explain the overdispersion observed in real single-cell RNA-seq data.

# C) Time Dependence of $e^{sx}$ in Probability Generating Functions

## 1. Is {{< math >}} $e^{sx}$ {{< /math >}} a Function of {{< math >}} $t$ {{< /math >}}?

The term {{< math >}} $e^{sx}$ {{< /math >}} depends on {{< math >}} $x$ {{< /math >}}, but not directly on time {{< math >}} $t$ {{< /math >}}.

The probability distribution {{< math >}} $P(x,t)$ {{< /math >}} changes over time, but the values {{< math >}} $x$ {{< /math >}} themselves remain the same discrete states.

Thus, while {{< math >}} $x$ {{< /math >}} indirectly depends on {{< math >}} $t$ {{< /math >}} through {{< math >}} $P(x,t)$ {{< /math >}}, {{< math >}} $e^{sx}$ {{< /math >}} itself is not explicitly time-dependent.

## 2. Why Is {{< math >}} $e^{sx}$ {{< /math >}} Treated as Constant?

When differentiating {{< math >}} $M(s,t)$ {{< /math >}}:

{{< math >}} 
$$ M(s,t) = \sum_{x=0}^{\infty} e^{sx} P(x,t) $$
{{< /math >}}

with respect to {{< math >}} $t$ {{< /math >}}, we apply the Leibniz rule for differentiation under summation:

{{< math >}} 
$$ \frac{d}{dt} M(s,t) = \sum_{x=0}^{\infty} e^{sx} \frac{d}{dt} P(x,t) $$
{{< /math >}}

Since {{< math >}} $e^{sx}$ {{< /math >}} does not explicitly depend on {{< math >}} $t$ {{< /math >}}, it remains constant during differentiation, allowing it to be factored out.

## 3. Intuition from a Physical Perspective

Think of {{< math >}} $e^{sx}$ {{< /math >}} as a weighting function applied to each discrete state {{< math >}} $x$ {{< /math >}}.

The system evolves by shifting probability mass between states {{< math >}} $P(x,t)$ {{< /math >}}.

Because {{< math >}} $x$ {{< /math >}} represents fixed molecular counts, and {{< math >}} $e^{sx}$ {{< /math >}} is defined over these counts, its value does not change dynamically.

# D) Derivation of the Second Uncentered Moment from the Chemical Master Equation

## Process Definition

We consider a simple birth-death process defined by:

- **Production (birth)**: {{< math >}} $0 \xrightarrow{\alpha} X$ {{< /math >}} (constant rate {{< math >}} $\alpha$ {{< /math >}})
- **Degradation (death)**: {{< math >}} $X \xrightarrow{\beta} 0$ {{< /math >}} (first-order rate {{< math >}} $\beta x$ {{< /math >}})

## Chemical Master Equation

The Chemical Master Equation (CME) for the probability {{< math >}} $P(x,t)$ {{< /math >}} of having {{< math >}} $x$ {{< /math >}} molecules at time {{< math >}} $t$ {{< /math >}} is:

{{< math >}} 
$$ \frac{dP(x,t)}{dt} = \alpha P(x-1,t) - \alpha P(x,t) + \beta(x+1)P(x+1,t) - \beta x P(x,t) \tag{1} $$
{{< /math >}}

*Note: For {{< math >}} $x=0$ {{< /math >}}, we assume {{< math >}} $P(-1,t) = 0$ {{< /math >}}. For {{< math >}} $x=0$ {{< /math >}}, the first term {{< math >}} $\alpha P(x-1,t)$ {{< /math >}} is considered {{< math >}} $\alpha P(-1,t) = 0$ {{< /math >}}. The last term {{< math >}} $\beta x P(x,t)$ {{< /math >}} is also 0 when {{< math >}} $x=0$ {{< /math >}}.*

## Derivation of {{< math >}} $\frac{d\langle x^2 \rangle}{dt}$ {{< /math >}}

The time derivative of the {{< math >}} $n$ {{< /math >}}-th moment {{< math >}} $\langle x^n \rangle$ {{< /math >}} is defined as:

{{< math >}} 
$$ \frac{d\langle x^n \rangle}{dt} = \sum_{x=0}^{\infty} x^n \frac{dP(x,t)}{dt} \tag{2} $$
{{< /math >}}

For the second moment, {{< math >}} $n=2$ {{< /math >}}:

{{< math >}} 
$$ \frac{d\langle x^2 \rangle}{dt} = \sum_{x=0}^{\infty} x^2 [\alpha P(x-1,t) - \alpha P(x,t) + \beta(x+1)P(x+1,t) - \beta x P(x,t)] \tag{3} $$
{{< /math >}}

We'll break this sum into terms related to {{< math >}} $\alpha$ {{< /math >}} (production) and terms related to {{< math >}} $\beta$ {{< /math >}} (degradation).

### 1. Production Terms (Terms with {{< math >}} $\alpha$ {{< /math >}})

The {{< math >}} $\alpha$ {{< /math >}} terms are:

{{< math >}} 
$$ T_\alpha = \sum_{x=0}^{\infty} x^2 [\alpha P(x-1,t) - \alpha P(x,t)] \tag{4} $$
{{< /math >}}

{{< math >}} 
$$ T_\alpha = \alpha \left[ \sum_{x=0}^{\infty} x^2 P(x-1,t) - \sum_{x=0}^{\infty} x^2 P(x,t) \right] \tag{5} $$
{{< /math >}}

The second sum is straightforward: {{< math >}} $\sum_{x=0}^{\infty} x^2 P(x,t) = \langle x^2 \rangle$ {{< /math >}}.

For the first sum, {{< math >}} $\sum_{x=0}^{\infty} x^2 P(x-1,t)$ {{< /math >}}:

Since {{< math >}} $P(-1,t) = 0$ {{< /math >}}, the term for {{< math >}} $x=0$ {{< /math >}} is 0. We can start the sum from {{< math >}} $x=1$ {{< /math >}}.

Let {{< math >}} $j = x-1$ {{< /math >}}. Then {{< math >}} $x = j+1$ {{< /math >}}. As {{< math >}} $x$ {{< /math >}} goes from 1 to {{< math >}} $\infty$ {{< /math >}}, {{< math >}} $j$ {{< /math >}} goes from 0 to {{< math >}} $\infty$ {{< /math >}}.

{{< math >}} 
$$ \sum_{x=1}^{\infty} x^2 P(x-1,t) = \sum_{j=0}^{\infty} (j+1)^2 P(j,t) \tag{6} $$
{{< /math >}}

{{< math >}} 
$$ = \sum_{j=0}^{\infty} (j^2 + 2j + 1)P(j,t) \tag{7} $$
{{< /math >}}

{{< math >}} 
$$ = \sum_{j=0}^{\infty} j^2 P(j,t) + 2\sum_{j=0}^{\infty} j P(j,t) + \sum_{j=0}^{\infty} P(j,t) \tag{8} $$
{{< /math >}}

Recognizing the definitions of moments:

{{< math >}} 
$$ = \langle x^2 \rangle + 2\langle x \rangle + 1 \tag{9} $$
{{< /math >}}

(Since {{< math >}} $\sum P(j,t) = 1$ {{< /math >}} because it's a probability distribution).

Combining the {{< math >}} $\alpha$ {{< /math >}} terms:

{{< math >}} 
$$ T_\alpha = \alpha[(\langle x^2 \rangle + 2\langle x \rangle + 1) - \langle x^2 \rangle] \tag{10} $$
{{< /math >}}

{{< math >}} 
$$ T_\alpha = \alpha(1 + 2\langle x \rangle) \tag{11} $$
{{< /math >}}

### 2. Degradation Terms (Terms with {{< math >}} $\beta$ {{< /math >}})

The {{< math >}} $\beta$ {{< /math >}} terms are:

{{< math >}} 
$$ T_\beta = \sum_{x=0}^{\infty} x^2 [\beta(x+1)P(x+1,t) - \beta x P(x,t)] \tag{12} $$
{{< /math >}}

{{< math >}} 
$$ T_\beta = \beta \left[ \sum_{x=0}^{\infty} x^2 (x+1)P(x+1,t) - \sum_{x=0}^{\infty} x^3 P(x,t) \right] \tag{13} $$
{{< /math >}}

The second sum is: {{< math >}} $\sum_{x=0}^{\infty} x^3 P(x,t) = \langle x^3 \rangle$ {{< /math >}}.

For the first sum, {{< math >}} $\sum_{x=0}^{\infty} x^2 (x+1)P(x+1,t)$ {{< /math >}}:

Let {{< math >}} $k = x+1$ {{< /math >}}. Then {{< math >}} $x = k-1$ {{< /math >}}. As {{< math >}} $x$ {{< /math >}} goes from 0 to {{< math >}} $\infty$ {{< /math >}}, {{< math >}} $k$ {{< /math >}} goes from 1 to {{< math >}} $\infty$ {{< /math >}}.

{{< math >}} 
$$ \sum_{k=1}^{\infty} (k-1)^2 k P(k,t) \tag{14} $$
{{< /math >}}

{{< math >}} 
$$ = \sum_{k=1}^{\infty} (k^2 - 2k + 1)k P(k,t) \tag{15} $$
{{< /math >}}

{{< math >}} 
$$ = \sum_{k=1}^{\infty} (k^3 - 2k^2 + k)P(k,t) \tag{16} $$
{{< /math >}}

Recognizing the definitions of moments (summing from {{< math >}} $k=0$ {{< /math >}} or {{< math >}} $k=1$ {{< /math >}} is equivalent for moments since terms at {{< math >}} $k=0$ {{< /math >}} with {{< math >}} $k, k^2, k^3$ {{< /math >}} are zero):

{{< math >}} 
$$ = \langle x^3 \rangle - 2\langle x^2 \rangle + \langle x \rangle \tag{17} $$
{{< /math >}}

Combining the {{< math >}} $\beta$ {{< /math >}} terms:

{{< math >}} 
$$ T_\beta = \beta[(\langle x^3 \rangle - 2\langle x^2 \rangle + \langle x \rangle) - \langle x^3 \rangle] \tag{18} $$
{{< /math >}}

{{< math >}} 
$$ T_\beta = \beta(-2\langle x^2 \rangle + \langle x \rangle) \tag{19} $$
{{< /math >}}

{{< math >}} 
$$ T_\beta = -\beta(2\langle x^2 \rangle - \langle x \rangle) \tag{20} $$
{{< /math >}}

### 3. Combining All Terms

Finally, add the simplified {{< math >}} $\alpha$ {{< /math >}} and {{< math >}} $\beta$ {{< /math >}} terms:

{{< math >}} 
$$ \frac{d\langle x^2 \rangle}{dt} = T_\alpha + T_\beta \tag{21} $$
{{< /math >}}

{{< math >}} 
$$ \frac{d\langle x^2 \rangle}{dt} = \alpha(1 + 2\langle x \rangle) - \beta(2\langle x^2 \rangle - \langle x \rangle) \tag{22} $$
{{< /math >}}