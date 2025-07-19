---
title: "Dyngen Propensity N Transcriptor"
date: 2025-07-10
draft: True
---

# A) Generalized Transcription Propensity Function for N Transcription Factors

This generalizes the transcription propensity function to the case where an arbitrary number {{< math >}}$N${{< /math >}} of transcription factors (TFs) regulate a gene by binding to its promoter. Here's a step-by-step breakdown:

## 🔹 Background

When N transcription factors can bind a promoter, there are {{< math >}}$2^N${{< /math >}} possible promoter states, depending on which combination of TFs are bound. Each state {{< math >}}$S_j${{< /math >}} has:

- A probability {{< math >}}$P(S_j)${{< /math >}}: how likely the promoter is in that state.
- A relative activation {{< math >}}$\alpha_j \in [0,1]${{< /math >}}: how active the promoter is in that state.

The total transcriptional propensity is then:

{{< math >}}
$$f(y_1, y_2, \ldots, y_N) = x_{\text{pr}} \cdot \sum_{j=0}^{2^N-1} \alpha_j \cdot P(S_j) \tag{1}$$
{{< /math >}}

Where {{< math >}}$x_{\text{pr}}${{< /math >}} is the maximal transcription rate.

## 🔹 Defining {{< math >}}$P(S_j)${{< /math >}}

Let:
- {{< math >}}$w_j${{< /math >}}: the unnormalized weight of state {{< math >}}$S_j${{< /math >}}, reflecting how likely the TFs bind in that configuration.

Then:

{{< math >}}
$$P(S_j) = \frac{w_j}{\sum_{k=0}^{2^N-1} w_k} \tag{2}$$
{{< /math >}}

## 🧠 How to compute {{< math >}}$w_j${{< /math >}}

Each TF {{< math >}}$H_i${{< /math >}} has a binding weight (or strength):

{{< math >}}
$$\nu_i = \left(\frac{y_i}{k_i}\right)^{n_i} \tag{3}$$
{{< /math >}}

where:
- {{< math >}}$y_i${{< /math >}}: concentration of TF {{< math >}}$H_i${{< /math >}}
- {{< math >}}$k_i${{< /math >}}: half-occupation constant
- {{< math >}}$n_i${{< /math >}}: Hill coefficient

Now for any state {{< math >}}$S_j${{< /math >}}, where {{< math >}}$j \in [0, 2^N-1]${{< /math >}}, you check which TFs are bound:

- Use binary digits of {{< math >}}$j${{< /math >}}: if the {{< math >}}$i${{< /math >}}-th bit is 1, TF {{< math >}}$H_i${{< /math >}} is bound

Define:

{{< math >}}
$$w_j = \prod_{i=1}^{N} \begin{cases}
\nu_i & \text{if TF } H_i \text{ is bound in state } S_j \\
1 & \text{otherwise}
\end{cases} \tag{4}$$
{{< /math >}}

Or more compactly using the notation in the question:

{{< math >}}
$$w_j = \prod_{i=1}^{N} \left[\text{if bit}_i(j) = 1 \Rightarrow \nu_i \text{ else } 1\right] \tag{5}$$
{{< /math >}}

So e.g.:
- {{< math >}}$j = 0${{< /math >}}: no TFs bound → {{< math >}}$w_0 = 1${{< /math >}}
- {{< math >}}$j = 1${{< /math >}} (binary 001): only {{< math >}}$H_1${{< /math >}} bound → {{< math >}}$w_1 = \nu_1${{< /math >}}
- {{< math >}}$j = 3${{< /math >}} (binary 011): {{< math >}}$H_1${{< /math >}} and {{< math >}}$H_2${{< /math >}} bound → {{< math >}}$w_3 = \nu_1 \cdot \nu_2${{< /math >}}

Then normalize:

{{< math >}}
$$P(S_j) = \frac{w_j}{\sum_{k=0}^{2^N-1} w_k} \tag{6}$$
{{< /math >}}

In the fully independent case, the denominator simplifies to:

{{< math >}}
$$\sum_j w_j = \prod_{i=1}^{N} (1 + \nu_i) \tag{7}$$
{{< /math >}}

So:

{{< math >}}
$$P(S_j) = \frac{w_j}{\prod_{i=1}^{N} (1 + \nu_i)} \tag{8}$$
{{< /math >}}

## 🔹 Final Form of Propensity

Now plug this into Equation (1):

{{< math >}}
$$f(y_1, \ldots, y_N) = x_{\text{pr}} \cdot \sum_{j=0}^{2^N-1} \alpha_j \cdot \frac{w_j}{\prod_{i=1}^{N} (1 + \nu_i)} \tag{9}$$
{{< /math >}}

Or:

{{< math >}}
$$f(y_1, \ldots, y_N) = x_{\text{pr}} \cdot \frac{1}{\prod_{i=1}^{N} (1 + \nu_i)} \cdot \sum_{j=0}^{2^N-1} \alpha_j \cdot w_j \tag{10}$$
{{< /math >}}

## 🔹 Biological Interpretation

- {{< math >}}$\alpha_j${{< /math >}}: models how active the promoter is in each binding state — could be additive, synergistic, etc.
- {{< math >}}$w_j${{< /math >}}: captures the binding probabilities of TFs to the promoter.

This setup flexibly models how combinations of TFs (with different concentrations and binding strengths) influence gene expression.


# B) Inductive Proof that {{< math >}} $$Z_N = \sum_{j=0}^{2^N-1} w_j = \prod_{i=1}^N (1+\nu_i)$$ {{< /math >}}

## 🧠 Definition Recap

Let:
- {{< math >}} $N$ {{< /math >}} be the number of transcription factors.

For any subset {{< math >}} $B_j \subseteq \{1,2,\ldots,N\}$ {{< /math >}}, define:

{{< math >}} 
$$w_j = \prod_{i \in B_j} \nu_i \tag{1}$$
{{< /math >}}

The total sum over all promoter states is:

{{< math >}} 
$$Z_N = \sum_{j=0}^{2^N-1} w_j = \sum_{B_j \subseteq \{1,\ldots,N\}} \prod_{i \in B_j} \nu_i \tag{2}$$
{{< /math >}}

## ✅ Base Case: {{< math >}} $N = 1$ {{< /math >}}

The two subsets of {{< math >}} $\{1\}$ {{< /math >}} are:

- {{< math >}} $\emptyset$ {{< /math >}} → {{< math >}} $w_0 = 1$ {{< /math >}}
- {{< math >}} $\{1\}$ {{< /math >}} → {{< math >}} $w_1 = \nu_1$ {{< /math >}}

Then:

{{< math >}} 
$$Z_1 = 1 + \nu_1 = \prod_{i=1}^1 (1+\nu_i) \tag{3}$$
{{< /math >}}

✔️ Base case holds.

## 🔁 Inductive Step

Assume the formula holds for {{< math >}} $N = k$ {{< /math >}}, i.e.:

{{< math >}} 
$$Z_k = \sum_{B_j \subseteq \{1,\ldots,k\}} \prod_{i \in B_j} \nu_i = \prod_{i=1}^k (1+\nu_i) \tag{4}$$
{{< /math >}}

### Add {{< math >}} $\nu_{k+1}$ {{< /math >}} (i.e., one more transcription factor)

Every state in {{< math >}} $Z_{k+1}$ {{< /math >}} corresponds to:

- either a state from {{< math >}} $Z_k$ {{< /math >}} (no TF {{< math >}} $k+1$ {{< /math >}}),
- or a state from {{< math >}} $Z_k$ {{< /math >}} with TF {{< math >}} $k+1$ {{< /math >}} also bound.

So:

{{< math >}} 
$$Z_{k+1} = \sum_{B \subseteq \{1,\ldots,k\}} \left(\prod_{i \in B} \nu_i\right) + \sum_{B \subseteq \{1,\ldots,k\}} \left(\prod_{i \in B} \nu_i\right) \cdot \nu_{k+1} \tag{5}$$
{{< /math >}}

Factor this:

{{< math >}} 
$$Z_{k+1} = \left(\sum_{B \subseteq \{1,\ldots,k\}} \prod_{i \in B} \nu_i\right) \cdot (1 + \nu_{k+1}) = Z_k \cdot (1 + \nu_{k+1}) \tag{6}$$
{{< /math >}}

Now plug in the inductive hypothesis:

{{< math >}} 
$$Z_{k+1} = \left(\prod_{i=1}^k (1+\nu_i)\right)(1+\nu_{k+1}) = \prod_{i=1}^{k+1} (1+\nu_i) \tag{7}$$
{{< /math >}}

✔️ Inductive step holds.

## ✅ Conclusion

By mathematical induction, for all {{< math >}} $N$ {{< /math >}},

{{< math >}} 
$$Z_N = \sum_{j=0}^{2^N-1} w_j = \prod_{i=1}^N (1+\nu_i) \tag{8}$$
{{< /math >}}

Let me know if you'd like to visualize this using binary state enumeration or connect it with logic gates or regulatory graphs.