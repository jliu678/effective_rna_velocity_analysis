---
title: "Gillespie Ssa"
date: 2025-07-10
draft: True
---

# Gillespie's SSA for Cyclic Reactions

## 1. Theory: Gillespie's SSA for Cyclic Reactions

### System Setup

Consider a cyclic set of chemical reactions:

{{< math >}} 
$$A \xrightarrow{k_1} B \xrightarrow{k_2} C \xrightarrow{k_3} A \tag{1}$$
{{< /math >}}

where:
- {{< math >}} $k_1, k_2, k_3$ {{< /math >}} are rate constants
- Molecule counts at time {{< math >}} $t$ {{< /math >}} are {{< math >}} $X_A(t), X_B(t), X_C(t)$ {{< /math >}}

### Propensity Functions

Each reaction {{< math >}} $R_j$ {{< /math >}} has a propensity function {{< math >}} $a_j(X)$ {{< /math >}}, which gives the instantaneous probability per unit time of that reaction firing, given state {{< math >}} $X$ {{< /math >}}:

{{< math >}} 
$$\begin{cases}
a_1(X) = k_1 X_A \\
a_2(X) = k_2 X_B \\
a_3(X) = k_3 X_C
\end{cases} \tag{2}$$
{{< /math >}}

### Gillespie's SSA: Key Ideas

The system evolves as a continuous-time Markov jump process.

Time to next reaction is exponentially distributed with rate:

{{< math >}} 
$$a_0 = a_1 + a_2 + a_3 \tag{3}$$
{{< /math >}}

Probability that next reaction is reaction {{< math >}} $j$ {{< /math >}} is:

{{< math >}} 
$$\frac{a_j}{a_0} \tag{4}$$
{{< /math >}}

### Sampling Time to Next Reaction

The waiting time {{< math >}} $\tau$ {{< /math >}} until the next reaction follows:

{{< math >}} 
$$P(\tau > t) = e^{-a_0 t} \tag{5}$$
{{< /math >}}

Sample {{< math >}} $\tau$ {{< /math >}} by inverse transform sampling:

{{< math >}} 
$$\tau = \frac{1}{a_0} \ln\left(\frac{1}{r_1}\right), \quad r_1 \sim U(0,1) \tag{6}$$
{{< /math >}}

### Sampling Which Reaction Fires

Choose reaction {{< math >}} $R_j$ {{< /math >}} so that:

{{< math >}} 
$$\sum_{i=1}^{j-1} a_i < r_2 a_0 \leq \sum_{i=1}^{j} a_i, \quad r_2 \sim U(0,1) \tag{7}$$
{{< /math >}}

### State Update

When reaction {{< math >}} $R_j$ {{< /math >}} fires, update molecule counts:

{{< math >}} 
$$\begin{align}
X_A &\rightarrow X_A + \nu_{A,j} \tag{8a}\\
X_B &\rightarrow X_B + \nu_{B,j} \tag{8b}\\
X_C &\rightarrow X_C + \nu_{C,j} \tag{8c}
\end{align}$$
{{< /math >}}

where {{< math >}} $\nu_{\cdot,j}$ {{< /math >}} is the stoichiometric change vector for reaction {{< math >}} $j$ {{< /math >}}:

{{< math >}} 
$$\begin{cases}
R_1: (X_A, X_B, X_C) \rightarrow (X_A - 1, X_B + 1, X_C) \\
R_2: (X_A, X_B, X_C) \rightarrow (X_A, X_B - 1, X_C + 1) \\
R_3: (X_A, X_B, X_C) \rightarrow (X_A + 1, X_B, X_C - 1)
\end{cases} \tag{9}$$
{{< /math >}}

### Interim Summary: Gillespie Algorithm for the Cycle

1. Initialize {{< math >}} $t = 0$ {{< /math >}}, state {{< math >}} $X$ {{< /math >}}
2. Compute propensities {{< math >}} $a_1, a_2, a_3$ {{< /math >}} and total {{< math >}} $a_0$ {{< /math >}}
3. Sample {{< math >}} $\tau = \frac{1}{a_0} \ln(1/r_1)$ {{< /math >}}, update time {{< math >}} $t \leftarrow t + \tau$ {{< /math >}}
4. Sample {{< math >}} $r_2$ {{< /math >}}, select reaction {{< math >}} $j$ {{< /math >}}
5. Update state by stoichiometry of {{< math >}} $R_j$ {{< /math >}}
6. Repeat until {{< math >}} $t$ {{< /math >}} reaches max time or stopping condition

## 2. Proof Why SSA Works for Cyclic Reactions

### Markov Jump Process Foundation

The system's state {{< math >}} $X$ {{< /math >}} evolves randomly with jump rates given by {{< math >}} $a_j(X)$ {{< /math >}}.

#### What is "Jump"?

In a stochastic chemical reaction system, the system's state is defined by the numbers of molecules of each species.

A jump means the system changes state by a reaction firing.

For example, if reaction {{< math >}} $R_1$ {{< /math >}} turns one molecule of A into one molecule of B, when {{< math >}} $R_1$ {{< /math >}} fires, the state "jumps" from {{< math >}} $(A,B,C)$ {{< /math >}} to {{< math >}} $(A-1,B+1,C)$ {{< /math >}}.

#### What is a Jump Rate?

The jump rate is the instantaneous probability per unit time that a particular reaction will occur, causing a jump from the current state.

Formally, for reaction {{< math >}} $R_j$ {{< /math >}}, the jump rate at state {{< math >}} $X$ {{< /math >}} is the propensity function {{< math >}} $a_j(X)$ {{< /math >}}.

{{< math >}} 
$$\text{Jump rate for reaction } R_j = a_j(X) \tag{1}$$
{{< /math >}}

This is sometimes called the reaction rate in stochastic modeling (not to be confused with deterministic rates).

#### Why "Jump"?

The system jumps from one discrete state (counts of molecules) to another when a reaction occurs.

Between jumps, the state is constant.

The jump rate governs the timing of these jumps probabilistically.

#### Example

For reaction {{< math >}} $R_1: A \rightarrow B$ {{< /math >}} with rate constant {{< math >}} $k_1$ {{< /math >}}, if you have {{< math >}} $X_A$ {{< /math >}} molecules of A:

{{< math >}} 
$$\text{Jump rate} = a_1(X) = k_1 \times X_A \tag{2}$$
{{< /math >}}

This means:

The probability that {{< math >}} $R_1$ {{< /math >}} fires in the next tiny time {{< math >}} $dt$ {{< /math >}} is approximately {{< math >}} $a_1(X)dt$ {{< /math >}}.

{{< math >}} 
$$P(\text{reaction } R_1 \text{ fires in time } dt) \approx a_1(X) \cdot dt \tag{3}$$
{{< /math >}}

Waiting time until next jump is exponentially distributed because the minimum of independent exponential waiting times (each reaction) is exponential with rate {{< math >}} $a_0$ {{< /math >}}.

The probability the minimum corresponds to reaction {{< math >}} $j$ {{< /math >}} is proportional to {{< math >}} $a_j$ {{< /math >}}.

### Intuition

Each reaction channel {{< math >}} $j$ {{< /math >}} can be imagined as an independent Poisson process with rate {{< math >}} $a_j$ {{< /math >}}.

The system "waits" for the first of these Poisson events.

Time until next reaction = minimum of all {{< math >}} $\tau_j \sim \text{Exponential}(a_j)$ {{< /math >}}.

The minimum of exponentials is exponential with rate {{< math >}} $a_0$ {{< /math >}}.

The probability reaction {{< math >}} $j$ {{< /math >}} occurs is the chance {{< math >}} $\tau_j$ {{< /math >}} is the minimum → proportional to {{< math >}} $a_j$ {{< /math >}}.

This general logic holds regardless of network topology (including cycles).

## 3. Numeric Example: Step-by-Step Simulation

### Parameters and Initial Conditions

| Reaction | Rate {{< math >}} $k$ {{< /math >}} | Initial Molecules |
|----------|-------------------------------------|-------------------|
| {{< math >}} $R_1: A \rightarrow B$ {{< /math >}} | 0.01 | {{< math >}} $A = 100$ {{< /math >}} |
| {{< math >}} $R_2: B \rightarrow C$ {{< /math >}} | 0.02 | {{< math >}} $B = 0$ {{< /math >}} |
| {{< math >}} $R_3: C \rightarrow A$ {{< /math >}} | 0.015 | {{< math >}} $C = 0$ {{< /math >}} |

### Step 1: Initial State

{{< math >}} $t = 0$ {{< /math >}}, state {{< math >}} $(A,B,C) = (100,0,0)$ {{< /math >}}

Calculate reaction propensities:

{{< math >}} 
$$a_1 = 0.01 \times 100 = 1.0 \tag{1}$$
{{< /math >}}

{{< math >}} 
$$a_2 = 0.02 \times 0 = 0 \tag{2}$$
{{< /math >}}

{{< math >}} 
$$a_3 = 0.015 \times 0 = 0 \tag{3}$$
{{< /math >}}

{{< math >}} 
$$a_0 = 1.0 \tag{4}$$
{{< /math >}}

Sample random numbers:

{{< math >}} 
$$r_1 = 0.5 \Rightarrow \tau = \frac{1}{1.0} \ln\left(\frac{1}{0.5}\right) = 0.693 \tag{5}$$
{{< /math >}}

{{< math >}} 
$$r_2 = 0.3 \Rightarrow r_2 a_0 = 0.3 \tag{6}$$
{{< /math >}}

Since {{< math >}} $0.3 < a_1 = 1.0$ {{< /math >}}, fire {{< math >}} $R_1$ {{< /math >}}:

New state: {{< math >}} $(A,B,C) = (99,1,0)$ {{< /math >}}

Time {{< math >}} $t = 0.693$ {{< /math >}}

### Step 2: Second Iteration

{{< math >}} $t = 0.693$ {{< /math >}}, state {{< math >}} $(99,1,0)$ {{< /math >}}

Calculate reaction propensities:

{{< math >}} 
$$a_1 = 0.01 \times 99 = 0.99 \tag{7}$$
{{< /math >}}

{{< math >}} 
$$a_2 = 0.02 \times 1 = 0.02 \tag{8}$$
{{< /math >}}

{{< math >}} 
$$a_3 = 0.015 \times 0 = 0 \tag{9}$$
{{< /math >}}

{{< math >}} 
$$a_0 = 1.01 \tag{10}$$
{{< /math >}}

Sample random numbers:

{{< math >}} 
$$r_1 = 0.1 \Rightarrow \tau = \frac{1}{1.01} \ln\left(\frac{1}{0.1}\right) \approx 2.28 \tag{11}$$
{{< /math >}}

{{< math >}} 
$$r_2 = 0.995 \Rightarrow r_2 a_0 = 1.00 \tag{12}$$
{{< /math >}}

Cumulative propensities:

{{< math >}} 
$$a_1 = 0.99 \tag{13}$$
{{< /math >}}

{{< math >}} 
$$a_1 + a_2 = 1.01 \tag{14}$$
{{< /math >}}

Since {{< math >}} $1.00$ {{< /math >}} is between {{< math >}} $0.99$ {{< /math >}} and {{< math >}} $1.01$ {{< /math >}}, fire {{< math >}} $R_2$ {{< /math >}}:

New state: {{< math >}} $(A,B,C) = (99,0,1)$ {{< /math >}}

Time {{< math >}} $t = 0.693 + 2.28 = 2.973$ {{< /math >}}

### Step 3: Third Iteration

{{< math >}} $t = 2.973$ {{< /math >}}, state {{< math >}} $(99,0,1)$ {{< /math >}}

Calculate reaction propensities:

{{< math >}} 
$$a_1 = 0.01 \times 99 = 0.99 \tag{15}$$
{{< /math >}}

{{< math >}} 
$$a_2 = 0.02 \times 0 = 0 \tag{16}$$
{{< /math >}}

{{< math >}} 
$$a_3 = 0.015 \times 1 = 0.015 \tag{17}$$
{{< /math >}}

{{< math >}} 
$$a_0 = 1.005 \tag{18}$$
{{< /math >}}

Sample random numbers:

{{< math >}} 
$$r_1 = 0.7 \Rightarrow \tau = \frac{1}{1.005} \ln\left(\frac{1}{0.7}\right) \approx 0.357 \tag{19}$$
{{< /math >}}

{{< math >}} 
$$r_2 = 0.98 \Rightarrow r_2 a_0 = 0.985 \tag{20}$$
{{< /math >}}

Cumulative sums:

{{< math >}} 
$$a_1 = 0.99 \tag{21}$$
{{< /math >}}

{{< math >}} 
$$a_1 + a_3 = 1.005 \tag{22}$$
{{< /math >}}

Since {{< math >}} $0.985 < 0.99$ {{< /math >}}, fire {{< math >}} $R_1$ {{< /math >}}:

New state: {{< math >}} $(A,B,C) = (98,1,1)$ {{< /math >}}

Time {{< math >}} $t = 2.973 + 0.357 = 3.33$ {{< /math >}}

### Step 4: Fourth Iteration

{{< math >}} $t = 3.33$ {{< /math >}}, state {{< math >}} $(98,1,1)$ {{< /math >}}

Calculate reaction propensities:

{{< math >}} 
$$a_1 = 0.01 \times 98 = 0.98 \tag{23}$$
{{< /math >}}

{{< math >}} 
$$a_2 = 0.02 \times 1 = 0.02 \tag{24}$$
{{< /math >}}

{{< math >}} 
$$a_3 = 0.015 \times 1 = 0.015 \tag{25}$$
{{< /math >}}

{{< math >}} 
$$a_0 = 1.015 \tag{26}$$
{{< /math >}}

Sample random numbers:

{{< math >}} 
$$r_1 = 0.2 \Rightarrow \tau = \frac{1}{1.015} \ln\left(\frac{1}{0.2}\right) \approx 1.59 \tag{27}$$
{{< /math >}}

{{< math >}} 
$$r_2 = 0.99 \Rightarrow r_2 a_0 = 1.004 \tag{28}$$
{{< /math >}}

Cumulative sums:

{{< math >}} 
$$a_1 = 0.98 \tag{29}$$
{{< /math >}}

{{< math >}} 
$$a_1 + a_2 = 1.00 \tag{30}$$
{{< /math >}}

{{< math >}} 
$$a_1 + a_2 + a_3 = 1.015 \tag{31}$$
{{< /math >}}

Since {{< math >}} $1.004$ {{< /math >}} is between {{< math >}} $1.00$ {{< /math >}} and {{< math >}} $1.015$ {{< /math >}}, fire {{< math >}} $R_3$ {{< /math >}}:

New state: {{< math >}} $(A,B,C) = (99,1,0)$ {{< /math >}}

Time {{< math >}} $t = 3.33 + 1.59 = 4.92$ {{< /math >}}

### Summary Table

| Step | Time {{< math >}} $t$ {{< /math >}} | {{< math >}} $A$ {{< /math >}} | {{< math >}} $B$ {{< /math >}} | {{< math >}} $C$ {{< /math >}} | Reaction Fired |
|------|-------------------------------------|--------------------------------|--------------------------------|--------------------------------|----------------|
| 1    | 0.693                               | 99                             | 1                              | 0                              | {{< math >}} $R_1$ {{< /math >}} |
| 2    | 2.973                               | 99                             | 0                              | 1                              | {{< math >}} $R_2$ {{< /math >}} |
| 3    | 3.33                                | 98                             | 1                              | 1                              | {{< math >}} $R_1$ {{< /math >}} |
| 4    | 4.92                                | 99                             | 1                              | 0                              | {{< math >}} $R_3$ {{< /math >}} |

### Optional: Python Snippet to Simulate

```python
import numpy as np

k = [0.01, 0.02, 0.015]
X = np.array([100, 0, 0])
t = 0.0
tmax = 10.0

while t < tmax:
    a = np.array([k[0]*X[0], k[1]*X[1], k[2]*X[2]])
    a0 = a.sum()
    if a0 == 0:
        break
    r1, r2 = np.random.rand(2)
    tau = (1/a0)*np.log(1/r1)
    t += tau
    cumsum = np.cumsum(a)
    reaction = np.searchsorted(cumsum, r2*a0)
    # Update state
    if reaction == 0:
        X += np.array([-1, +1, 0])
    elif reaction == 1:
        X += np.array([0, -1, +1])
    else:
        X += np.array([+1, 0, -1])
    print(f"t={t:.3f}, state={X}, reaction=R{reaction+1}")
```


# Why No Events Before Time t in a Poisson Process Has Probability e^(-λt)

## The Question

Why is the probability that no event occurs before time {{< math >}} $t$ {{< /math >}} in a Poisson process:

{{< math >}} 
$$P(\text{no event before } t) = e^{-\lambda t} \tag{1}$$
{{< /math >}}

Let's explain why this is true, both intuitively and mathematically.

## 🔁 Intuition: What is a Poisson Process?

A Poisson process models events that:

- Happen randomly in time
- Are independent of each other  
- Happen with a constant average rate {{< math >}} $\lambda$ {{< /math >}} per unit time

We're asking: What is the probability that no events occur in time interval {{< math >}} $[0,t]$ {{< /math >}}?

## ✅ Mathematical Derivation

**Fact:** In a Poisson process with rate {{< math >}} $\lambda$ {{< /math >}}, the number of events in time {{< math >}} $t$ {{< /math >}} follows a Poisson distribution:

{{< math >}} 
$$P(N(t) = k) = \frac{(\lambda t)^k}{k!} e^{-\lambda t} \tag{2}$$
{{< /math >}}

So the probability of zero events by time {{< math >}} $t$ {{< /math >}} is:

{{< math >}} 
$$P(N(t) = 0) = \frac{(\lambda t)^0}{0!} e^{-\lambda t} = e^{-\lambda t} \tag{3}$$
{{< /math >}}

That's it! That's where the formula comes from.

## 🧠 Interpretation

{{< math >}} $e^{-\lambda t}$ {{< /math >}} is the probability that nothing happens in time {{< math >}} $t$ {{< /math >}}

This is also the probability that the first event happens after time {{< math >}} $t$ {{< /math >}}

## 🔄 Relationship to Waiting Time

Let {{< math >}} $\tau$ {{< /math >}} be the waiting time to the first event. Then:

{{< math >}} 
$$P(\tau > t) = P(\text{no event in } [0,t]) = e^{-\lambda t} \tag{4}$$
{{< /math >}}

So the waiting time to the first event is exponentially distributed:

{{< math >}} 
$$f_\tau(t) = \lambda e^{-\lambda t} \tag{5}$$
{{< /math >}}

## ⏳ Real-World Analogy

Imagine a radioactive atom decaying. The chance that it survives without decaying for time {{< math >}} $t$ {{< /math >}} is:

{{< math >}} 
$$\text{Survival probability} = e^{-\lambda t} \tag{6}$$
{{< /math >}}

This is why the exponential distribution is also used to model lifetimes and failure times.

# How to sample τ from an Exponential Distribution

We want to sample the time {{< math >}} $\tau$ {{< /math >}} until the next event, where {{< math >}} $\tau \sim \text{Exponential}(a_0)$ {{< /math >}}, and {{< math >}} $a_0$ {{< /math >}} is the total propensity (i.e. the total rate at which any reaction occurs).

The cumulative distribution function (CDF) of the exponential distribution is:

{{< math >}} 
$$P(\tau \leq t) = 1 - e^{-a_0 t} \tag{7}$$
{{< /math >}}

So the probability that τ is greater than t is:

{{< math >}} 
$$P(\tau > t) = e^{-a_0 t} \tag{8}$$
{{< /math >}}

This tells us: The longer we wait, the less likely it is that no reaction has occurred.

We use inverse transform sampling, which works like this:

**Step 1:** Generate a uniform random number:

{{< math >}} 
$$r_1 \sim \text{Uniform}(0,1) \tag{9}$$
{{< /math >}}

**Step 2:** Set it equal to the survival function:

{{< math >}} 
$$r_1 = P(\tau > t) = e^{-a_0 t} \tag{10}$$
{{< /math >}}

**Step 3:** Solve for {{< math >}} $t$ {{< /math >}}:

{{< math >}} 
$$\ln r_1 = -a_0 t \Rightarrow \tau = \frac{1}{a_0} \ln\left(\frac{1}{r_1}\right) \tag{11}$$
{{< /math >}}

This gives a sample {{< math >}} $\tau$ {{< /math >}} from the correct exponential distribution.

## ✅ Why This Is Reasonable

- The exponential distribution exactly describes the waiting time to the first event in a Poisson process (which is what Gillespie's SSA models).
- {{< math >}} $r_1 \sim \text{Uniform}(0,1)$ {{< /math >}} is used to generate a random waiting time from that distribution.
- The math guarantees that the generated {{< math >}} $\tau$ {{< /math >}} has the right probability density:

{{< math >}} 
$$f(\tau) = a_0 e^{-a_0 \tau} \tag{12}$$
{{< /math >}}

# Proof: Sampling Formula Gives Exponential Distribution

## 🎯 Goal

Prove that sampling:

{{< math >}} 
$$\tau = \frac{1}{a_0} \ln\left(\frac{1}{r_1}\right) \text{ where } r_1 \sim \text{Uniform}(0,1) \tag{13}$$
{{< /math >}}

gives a random variable that follows the exponential distribution with rate {{< math >}} $a_0$ {{< /math >}}:

{{< math >}} 
$$f_\tau(t) = a_0 e^{-a_0 t} \text{ for } t \geq 0 \tag{14}$$
{{< /math >}}

## 🧮 Step 1: Start with the Transformation

Let:

{{< math >}} 
$$\tau = -\frac{1}{a_0} \ln(r_1) \text{ where } r_1 \sim \text{Uniform}(0,1) \tag{15}$$
{{< /math >}}

We'll derive the probability density function (PDF) of {{< math >}} $\tau$ {{< /math >}} using change of variables.

## 🔁 Step 2: Compute the CDF of τ

Let's compute the cumulative distribution function {{< math >}} $F_\tau(t)$ {{< /math >}}:

{{< math >}} 
$$F_\tau(t) = P(\tau \leq t) = P\left(-\frac{1}{a_0} \ln(r_1) \leq t\right) \tag{16}$$
{{< /math >}}

{{< math >}} 
$$= P(\ln(r_1) \geq -a_0 t) = P(r_1 \geq e^{-a_0 t}) \tag{17}$$
{{< /math >}}

Since {{< math >}} $r_1 \sim \text{Uniform}(0,1)$ {{< /math >}}, we have:

{{< math >}} 
$$P(r_1 \geq e^{-a_0 t}) = 1 - e^{-a_0 t} \tag{18}$$
{{< /math >}}

So the CDF of {{< math >}} $\tau$ {{< /math >}} is:

{{< math >}} 
$$F_\tau(t) = 1 - e^{-a_0 t}, \quad t \geq 0 \tag{19}$$
{{< /math >}}

## 🧠 Step 3: Differentiate to Get the PDF

Now differentiate the CDF to get the PDF:

{{< math >}} 
$$f_\tau(t) = \frac{d}{dt} F_\tau(t) = \frac{d}{dt}(1 - e^{-a_0 t}) = a_0 e^{-a_0 t} \tag{20}$$
{{< /math >}}

Which is exactly the exponential distribution with rate {{< math >}} $a_0$ {{< /math >}}!

## ✅ Conclusion

We've shown:

1. The transformation {{< math >}} $\tau = -\frac{1}{a_0} \ln(r_1)$ {{< /math >}} gives a variable with exponential distribution.

2. This matches the required PDF:

{{< math >}} 
$$f(\tau) = a_0 e^{-a_0 \tau} \tag{21}$$
{{< /math >}}

Therefore, Gillespie's formula for sampling waiting time is mathematically justified.