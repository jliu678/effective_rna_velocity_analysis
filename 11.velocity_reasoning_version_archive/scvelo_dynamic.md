---
title: "Scvelo Dynamic"
date: 2025-07-10
draft: True
---

# Dynamic Solutions for RNA Velocity Model

To derive the solution for s(t) from the given system of differential equations, we will follow a two-step process:

1. Solve the first differential equation for u(t) using the initial condition u(t₀)=u₀.
2. Substitute the expression for u(t) into the second differential equation and then solve for s(t) using the initial condition s(t₀)=s₀.

## Step 1: Solve for u(t)

The first differential equation is:

{{< math >}}  
$$
\frac{du(t)}{dt}=\alpha-\beta u(t)
$$  
{{< /math >}}

Rearrange it into the standard form for a first-order linear ODE:

{{< math >}}  
$$
\frac{du(t)}{dt}+\beta u(t)=\alpha
$$  
{{< /math >}}

This is a linear first-order differential equation of the form $\frac{dy}{dx}+P(x)y=Q(x)$, where P(t)=β and Q(t)=α.
The integrating factor (IF) is $e^{\int P(t)dt}=e^{\int \beta dt}=e^{\beta t}$.

Multiply the entire equation by the integrating factor:

{{< math >}}  
$$
e^{\beta t}\frac{du(t)}{dt}+\beta u(t)e^{\beta t}=\alpha e^{\beta t}
$$  
{{< /math >}}

The left side is the derivative of the product $(u(t)e^{\beta t})$:

{{< math >}}  
$$
\frac{d}{dt}(u(t)e^{\beta t})=\alpha e^{\beta t}
$$  
{{< /math >}}

Integrate both sides with respect to t:

{{< math >}}  
$$
\int \frac{d}{dt}(u(t)e^{\beta t})dt=\int \alpha e^{\beta t}dt
$$
$$
u(t)e^{\beta t}=\frac{\alpha}{\beta}e^{\beta t}+C_1
$$  
{{< /math >}}

Now, solve for u(t):

{{< math >}}  
$$
u(t)=\frac{\alpha}{\beta}+C_1e^{-\beta t}
$$  
{{< /math >}}

Use the initial condition u(t₀)=u₀ to find C₁:

{{< math >}}  
$$
u_0=\frac{\alpha}{\beta}+C_1e^{-\beta t_0}
$$
$$
C_1e^{-\beta t_0}=u_0-\frac{\alpha}{\beta}
$$
$$
C_1=(u_0-\frac{\alpha}{\beta})e^{\beta t_0}
$$  
{{< /math >}}


Substitute C₁ back into the equation for u(t):

{{< math >}}  
$$
u(t)=\frac{\alpha}{\beta}+(u_0-\frac{\alpha}{\beta})e^{\beta t_0}e^{-\beta t}
$$
$$
u(t)=\frac{\alpha}{\beta}+(u_0-\frac{\alpha}{\beta})e^{-\beta(t-t_0)}
$$  
{{< /math >}}

## Step 2: Solve for s(t)

Now that we have u(t), we can solve the second differential equation:

{{< math >}}  
$$
\frac{ds(t)}{dt}=\beta u(t)-\gamma s(t)
$$  
{{< /math >}}

Substitute the expression for u(t):

{{< math >}}  
$$
\frac{ds(t)}{dt}=\beta[\frac{\alpha}{\beta}+(u_0-\frac{\alpha}{\beta})e^{-\beta(t-t_0)}]-\gamma s(t)
$$
$$
\frac{ds(t)}{dt}+\gamma s(t)=\alpha+(u_0\beta-\alpha)e^{-\beta(t-t_0)}
$$  
{{< /math >}}

This is again a linear first-order differential equation. The integrating factor is $e^{\gamma t}$.

Multiply both sides by the integrating factor:

{{< math >}}  
$$
\frac{d}{dt}(s(t)e^{\gamma t})=\alpha e^{\gamma t}+(u_0\beta-\alpha)e^{\gamma t}e^{-\beta(t-t_0)}
$$  
{{< /math >}}

Integrate both sides:

{{< math >}}  
$$
s(t)e^{\gamma t}=\frac{\alpha}{\gamma}e^{\gamma t}+(u_0\beta-\alpha)\int e^{\gamma t}e^{-\beta(t-t_0)}dt+C_2
$$  
{{< /math >}}

Simplify the integral:

{{< math >}}  
$$
s(t)e^{\gamma t}=\frac{\alpha}{\gamma}e^{\gamma t}+(u_0\beta-\alpha)e^{\beta t_0}\int e^{(\gamma-\beta)t}dt+C_2
$$
$$
s(t)e^{\gamma t}=\frac{\alpha}{\gamma}e^{\gamma t}+(u_0\beta-\alpha)e^{\beta t_0}\frac{e^{(\gamma-\beta)t}}{\gamma-\beta}+C_2
$$  
{{< /math >}}

Solve for s(t):

{{< math >}}  
$$
s(t)=\frac{\alpha}{\gamma}+(u_0\beta-\alpha)e^{\beta t_0}\frac{e^{-\beta t}}{\gamma-\beta}+C_2e^{-\gamma t}
$$  
{{< /math >}}

Use the initial condition s(t₀)=s₀ to find C₂:

{{< math >}}  
$$
s_0=\frac{\alpha}{\gamma}+(u_0\beta-\alpha)\frac{e^{-\beta t_0+\beta t_0}}{\gamma-\beta}+C_2e^{-\gamma t_0}
$$
$$
C_2e^{-\gamma t_0}=s_0-\frac{\alpha}{\gamma}-(u_0\beta-\alpha)\frac{1}{\gamma-\beta}
$$
$$
C_2=[s_0-\frac{\alpha}{\gamma}-(u_0\beta-\alpha)\frac{1}{\gamma-\beta}]e^{\gamma t_0}
$$  
{{< /math >}}

The complete solution for s(t) is:

{{< math >}}  
$$
s(t)=\frac{\alpha}{\gamma}+\frac{u_0\beta-\alpha}{\gamma-\beta}e^{-\beta(t-t_0)}+[s_0-\frac{\alpha}{\gamma}-\frac{u_0\beta-\alpha}{\gamma-\beta}]e^{-\gamma(t-t_0)}
$$  
{{< /math >}}

### Special case: If γ=β

If γ=β, the integral becomes:

{{< math >}}  
$$
\int(βu_0-α)e^{(γ-β)t+βt_0}dt=(βu_0-α)e^{βt_0}\int dt=(βu_0-α)e^{βt_0}t
$$  
{{< /math >}}

Then,

{{< math >}}  
$$
s(t)e^{γt}=\frac{α}{γ}e^{γt}+(βu_0-α)e^{βt_0}t+C_2
$$
$$
s(t)=\frac{α}{γ}+(βu_0-α)te^{-γt+βt_0}+C_2e^{-γt}
$$  
{{< /math >}}

Since γ=β:

{{< math >}}  
$$
s(t)=\frac{α}{γ}+(γu_0-α)te^{-γ(t-t_0)}+C_2e^{-γt}
$$
$$
s(t)=\frac{α}{γ}+(γu_0-α)(t_0+τ)e^{-γτ}+C_2e^{-γ(t_0+τ)}
$$  
{{< /math >}}

This scenario is not covered by the target equation, which explicitly has a denominator (γ−β), implying γ≠β.



