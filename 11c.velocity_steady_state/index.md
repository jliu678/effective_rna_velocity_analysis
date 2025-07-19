---
title: 🧬 Math derivation for steady-state RNA velocity model
summary: The steady‑state model was the first to enable a mathematical estimation of RNA velocity, and most subsequent methods are modified versions of it or its generalization (the dynamic model in `scVelo`; see our other blogs). It has its limitations and a solid understanding of its underlying mathematics is needed to apply the model effectively. Here, we derive the steady-state model in `scVelo` and `velocyto`.
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, Math
  - steady-state model
  - scVelo
  - velocyto
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---
The steady‑state model was the first to enable a mathematical estimation of RNA velocity, and most subsequent methods are modified versions of it or its generalization (the dynamic model in `scVelo`; see our other blogs). It has its limitations and a solid understanding of its underlying mathematics is needed to apply the model effectively. Here, we derive the steady-state model in `scVelo` and `velocyto`.

### A. Estimating RNA velocity Using Steady-State Ratio in `scVelo` and `velocyto`

For steady-state models in `scVelo` and `velocyto`, RNA velocities are computed as deviations from this steady-state ratio:

{{< math >}} $$ \nu_i = u_i - \gamma_0 s_i  $$ {{< /math >}}

where 
{{< math >}}  
$$
\gamma_0 = \frac{\beta}{\gamma}
$$  
{{< /math >}}  
(i.e. the expected ratio of unspliced to spliced mRNA) is estimated analytically using least squares regression. Lower and upper quantiles in phase space, that is, where mRNA levels
reach minimum and maximum expression, respectively were assumed to be steady states and used for the estimation. Hence, the ratio can be approximated by a linear regression on these extreme quantiles.

Specifically, {{< math >}}$ \gamma_0 ${{< /math >}} is estimated as :

{{< math >}}  
$$
\gamma_0 = \frac{\mathbf{u}^\top \mathbf{s}}{\lVert \mathbf{s} \rVert^2}
$$  
{{< /math >}}

This is a least-squares fit of the form:

{{< math >}}  
$$
\mathbf{u} \approx \gamma_0 \cdot \mathbf{s}
$$  
{{< /math >}}

This means they model unspliced mRNA ({{< math >}}$\mathbf{u}${{< /math >}}) as being linearly dependent on spliced mRNA ({{< math >}}$\mathbf{s}${{< /math >}}), and solve for the slope {{< math >}}$\gamma_0${{< /math >}}.

---

#### Symbol Definition

- {{< math >}}$\mathbf{u} = (u_1, \ldots, u_n)${{< /math >}}: a vector of size-normalized unspliced mRNA counts for a single gene across a subset of cells.  
- {{< math >}}$\mathbf{s} = (s_1, \ldots, s_n)${{< /math >}}: a vector of corresponding spliced mRNA counts for that gene in the same subset of cells.  
- {{< math >}}$\mathbf{u}^\top \mathbf{s}${{< /math >}}: the dot product between unspliced and spliced vectors.  
- {{< math >}}$\lVert \mathbf{s} \rVert^2 = \mathbf{s}^\top \mathbf{s}${{< /math >}}: the squared norm of the spliced vector.  
- {{< math >}}$\gamma_0${{< /math >}}: the slope of the best-fit line through the origin relating {{< math >}}$\mathbf{u}${{< /math >}} to {{< math >}}$\mathbf{s}${{< /math >}}, i.e., the steady-state ratio.  

---

### B. Derive Least Squares Solution for Estimating the Steady-State Ratio $\gamma_0$

Let’s walk through the full derivation of the least squares solution for estimating the steady-state ratio {{< math >}}$\gamma_0${{< /math >}}, which minimizes the squared difference between the unspliced and spliced counts scaled by {{< math >}}$\gamma_0${{< /math >}}.

---

#### 🔧 Problem Setup

Given:

- {{< math >}}$u \in \mathbb{R}^n${{< /math >}}: vector of unspliced counts across {{< math >}}$n${{< /math >}} cells  
- {{< math >}}$s \in \mathbb{R}^n${{< /math >}}: vector of spliced counts across the same {{< math >}}$n${{< /math >}} cells

We model the unspliced counts as linearly proportional to spliced counts:

{{< math >}}
$$
u \approx \gamma_0 s
$$
{{< /math >}}

We want to find {{< math >}}$\gamma_0${{< /math >}} that minimizes the squared error:

{{< math >}}
$$
\min_{\gamma_0} \| u - \gamma_0 s \|^2
$$
{{< /math >}}

---

#### 🧮 Expand the Norm

{{< math >}}
$$
\| u - \gamma_0 s \|^2 = (u - \gamma_0 s)^\top (u - \gamma_0 s)
$$
{{< /math >}}

Expanding this expression:

{{< math >}}
$$
= u^\top u - 2\gamma_0 u^\top s + \gamma_0^2 s^\top s
$$
{{< /math >}}

---

#### 📉 Minimize with Respect to $\gamma_0$

Take derivative with respect to {{< math >}}$\gamma_0${{< /math >}} and set it to zero:

{{< math >}}
$$
\frac{d}{d\gamma_0} \left( \| u - \gamma_0 s \|^2 \right) = -2 u^\top s + 2\gamma_0 s^\top s = 0
$$
{{< /math >}}

Solve for {{< math >}}$\gamma_0${{< /math >}}:

{{< math >}}
$$
\gamma_0 = \frac{u^\top s}{\|s\|^2}
$$
{{< /math >}}

Where:

- {{< math >}}$u^\top s${{< /math >}} is the dot product  
- {{< math >}}$\|s\|^2 = s^\top s${{< /math >}} is the squared norm

---

#### ✅ Interpretation

This is the slope of the best-fit line (through the origin) predicting unspliced from spliced counts. It's equivalent to projecting {{< math >}}$u${{< /math >}} onto {{< math >}}$s${{< /math >}} in vector space.

In RNA velocity, this slope {{< math >}}$\gamma_0${{< /math >}} serves as a reference ratio under the steady-state assumption:

**Expected:**

{{< math >}}
$$
u = \gamma_0 s
$$
{{< /math >}}

Any deviation from this in other cells suggests that the gene is being upregulated or downregulated.

### C. Extending the Steady-State Model with Offset

#### 🧩 Original Model (No Offset)

Previously, we assumed:

{{< math >}}
$$u_i \approx \gamma_0 \cdot s_i$$
{{< /math >}}

But this forces the regression line to go through the origin, which isn't always biologically realistic.

#### ✅ Generalized Model (With Offset)

Now we model:

{{< math >}}
$$u_i \approx \gamma_0 \cdot s_i + o$$
{{< /math >}}

Where:

- {{< math >}}$\gamma_0${{< /math >}} is still the slope (steady-state ratio)
- {{< math >}}$o${{< /math >}} is a constant offset (intercept), modeling basal transcription

#### 🔍 Least Squares Solution with Offset

To solve for both {{< math >}}$\gamma_0${{< /math >}} and {{< math >}}$o${{< /math >}}, we use ordinary least squares (OLS) for linear regression.

Given data {{< math >}}$u=(u_1,\ldots,u_n)${{< /math >}}, {{< math >}}$s=(s_1,\ldots,s_n)${{< /math >}}, we estimate:

{{< math >}}
$$\hat{u} = \gamma_0 s + o$$
{{< /math >}}

OLS gives:


##### Slope (Steady-State Ratio):

{{< math >}}
$$\gamma_0 = \frac{\text{Cov}(u,s)}{\text{Var}(s)}$$
{{< /math >}}

Where:

{{< math >}}
$$\text{Cov}(u,s) = \frac{1}{n}\sum_{i=1}^n (u_i-\bar{u})(s_i-\bar{s})$$
{{< /math >}}

{{< math >}}
$$\text{Var}(s) = \frac{1}{n}\sum_{i=1}^n (s_i-\bar{s})^2$$
{{< /math >}}

##### Offset:

{{< math >}}
$$o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}

This centers the regression line at the mean of the data points.

### D. Derive Least Squares Solution with Intercept

Let's derive the least squares solution with an intercept (offset) step by step. Our goal is to estimate both:

- {{< math >}}$\gamma_0${{< /math >}}: the slope (steady-state ratio), and
- {{< math >}}$o${{< /math >}}: the offset (basal transcription level)

given that we model:

{{< math >}}
$$u_i = \gamma_0 s_i + o + \varepsilon_i$$
{{< /math >}}

for each cell {{< math >}}$i${{< /math >}}, where {{< math >}}$\varepsilon_i${{< /math >}} is the residual.

#### 🧩 1. Problem Setup

Let {{< math >}}$u = [u_1, u_2, ..., u_n]^T${{< /math >}}, and {{< math >}}$s = [s_1, s_2, ..., s_n]^T${{< /math >}}

We want to solve for {{< math >}}$\gamma_0${{< /math >}} and {{< math >}}$o${{< /math >}} that minimize the squared residuals:

{{< math >}}
$$\min_{\gamma_0,o} \sum_{i=1}^n (u_i - \gamma_0 s_i - o)^2$$
{{< /math >}}

This is a linear least squares regression with intercept.

#### 🧮 2. Matrix Form

Rewriting the model:

{{< math >}}
$$u = \gamma_0 s + o \cdot 1 + \varepsilon$$
{{< /math >}}

Let's construct a design matrix {{< math >}}$X${{< /math >}}:

{{< math >}}
$$X = \begin{bmatrix} s_1 & 1 \\ s_2 & 1 \\ \vdots & \vdots \\ s_n & 1 \end{bmatrix} \in \mathbb{R}^{n \times 2}$$
{{< /math >}}

Then the model becomes:

{{< math >}}
$$u = X \cdot \begin{bmatrix} \gamma_0 \\ o \end{bmatrix} + \varepsilon$$
{{< /math >}}

The least squares solution is given by:

{{< math >}}
$$\begin{bmatrix} \gamma_0 \\ o \end{bmatrix} = (X^T X)^{-1} X^T u$$
{{< /math >}}

We can compute this explicitly to get formulas for {{< math >}}$\gamma_0${{< /math >}} and {{< math >}}$o${{< /math >}} as below.

##### Slope:

{{< math >}}
$$\gamma_0 = \frac{\text{Cov}(u,s)}{\text{Var}(s)} = \frac{\sum_{i=1}^n (u_i-\bar{u})(s_i-\bar{s})}{\sum_{i=1}^n (s_i-\bar{s})^2}$$
{{< /math >}}

##### Offset:

{{< math >}}
$$o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}

These are the well-known results from simple linear regression that center the regression line at the mean of the data points.

#### 🎯 3. Matrix Derivation of the Least Squares Solution

Let's explicitly derive the solution:

{{< math >}}
$$\begin{bmatrix} \gamma_0 \\ o \end{bmatrix} = (X^\top X)^{-1} X^\top u$$
{{< /math >}}

for the model:

{{< math >}}
$$u = \gamma_0 s + o \cdot 1 + \varepsilon$$
{{< /math >}}

##### 3.1 Define the Design Matrix and Vectors

Let's say we have {{< math >}}$n${{< /math >}} samples (cells), with:

- {{< math >}}$s \in \mathbb{R}^n${{< /math >}}: spliced counts
- {{< math >}}$u \in \mathbb{R}^n${{< /math >}}: unspliced counts
- {{< math >}}$1 \in \mathbb{R}^n${{< /math >}}: vector of ones

We define the design matrix {{< math >}}$X \in \mathbb{R}^{n \times 2}${{< /math >}}:

{{< math >}}
$$X = \begin{bmatrix} s_1 & 1 \\ s_2 & 1 \\ \vdots & \vdots \\ s_n & 1 \end{bmatrix} = [s \quad 1]$$
{{< /math >}}

We write the model:

{{< math >}}
$$u = X\theta + \varepsilon, \text{ with } \theta = \begin{bmatrix} \gamma_0 \\ o \end{bmatrix}$$
{{< /math >}}

We want to minimize the residual sum of squares:

{{< math >}}
$$\min_\theta \|u - X\theta\|^2$$
{{< /math >}}

##### 3.2 Normal Equations

We derive the optimal {{< math >}}$\theta${{< /math >}} by solving the normal equations:

{{< math >}}
$$X^\top X\theta = X^\top u$$
{{< /math >}}

So,

{{< math >}}
$$\theta = (X^\top X)^{-1} X^\top u$$
{{< /math >}}

Let's compute each term step by step.

##### 3.3 Compute $X^\top X$

{{< math >}}
$$X^\top X = \begin{bmatrix} \sum s_i^2 & \sum s_i \\ \sum s_i & n \end{bmatrix}$$
{{< /math >}}

Let's denote:

- {{< math >}}$S = \sum_{i=1}^n s_i${{< /math >}}
- {{< math >}}$S_2 = \sum_{i=1}^n s_i^2${{< /math >}}
- {{< math >}}$n${{< /math >}}: number of cells

So:

{{< math >}}
$$X^\top X = \begin{bmatrix} S_2 & S \\ S & n \end{bmatrix}$$
{{< /math >}}

##### 3.4 Compute $X^\top u$

{{< math >}}
$$X^\top u = \begin{bmatrix} \sum s_i u_i \\ \sum u_i \end{bmatrix}$$
{{< /math >}}

Let's denote:

- {{< math >}}$U = \sum_{i=1}^n u_i${{< /math >}}
- {{< math >}}$SU = \sum_{i=1}^n s_i u_i${{< /math >}}

So:

{{< math >}}
$$X^\top u = \begin{bmatrix} SU \\ U \end{bmatrix}$$
{{< /math >}}

##### 3.5 Solve $\theta = (X^\top X)^{-1} X^\top u$

The inverse of a {{< math >}}$2\times2${{< /math >}} matrix:

{{< math >}}
$$\begin{bmatrix} a & b \\ c & d \end{bmatrix}^{-1} = \frac{1}{ad-bc}\begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$$
{{< /math >}}

So:

{{< math >}}
$$(X^\top X)^{-1} = \frac{1}{S_2n - S^2}\begin{bmatrix} n & -S \\ -S & S_2 \end{bmatrix}$$
{{< /math >}}

Now compute:

{{< math >}}
$$\theta = \frac{1}{S_2n - S^2}\begin{bmatrix} n & -S \\ -S & S_2 \end{bmatrix}\begin{bmatrix} SU \\ U \end{bmatrix}$$
{{< /math >}}

First row:
{{< math >}}
$$n\cdot SU - S\cdot U$$
{{< /math >}}

Second row:
{{< math >}}
$$-S\cdot SU + S_2\cdot U$$
{{< /math >}}

##### 3.6 Final Expression

So we get:

{{< math >}}
$$\gamma_0 = \frac{n\cdot \sum s_i u_i - \sum s_i \cdot \sum u_i}{n\cdot \sum s_i^2 - (\sum s_i)^2} = \frac{\text{Cov}(u,s)}{\text{Var}(s)}$$
{{< /math >}}

{{< math >}}
$$o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}

This completes the derivation of the least squares solution with intercept, showing how to arrive at the covariance and variance expressions.

### Why Normal Equations in 3.2 in Section D give the optimal solution

Let's prove that the normal equations

{{< math >}}
$$X^\top X\theta = X^\top u$$
{{< /math >}}

give us the optimal solution in least squares regression.

#### Goal of Least Squares

We want to minimize the sum of squared errors (residuals):

{{< math >}}
$$L(\theta) = \|u-X\theta\|^2 = (u-X\theta)^\top(u-X\theta)$$
{{< /math >}}

This is a quadratic loss function. To minimize it, we take its gradient with respect to {{< math >}}$\theta${{< /math >}}, set it to zero, and solve.

#### Step-by-step Derivation

Start with:

{{< math >}}
$$L(\theta) = (u-X\theta)^\top(u-X\theta)$$
{{< /math >}}

Expand the product:

{{< math >}}
$$L(\theta) = u^\top u - 2\theta^\top X^\top u + \theta^\top X^\top X\theta$$
{{< /math >}}

Now take the gradient with respect to {{< math >}}$\theta${{< /math >}}:

{{< math >}}
$$\nabla_\theta L = -2X^\top u + 2X^\top X\theta \tag{gradient} $$
{{< /math >}}

Set the gradient to zero:

{{< math >}}
$$-2X^\top u + 2X^\top X\theta = 0$$
{{< /math >}}

Divide both sides by 2:

{{< math >}}
$$X^\top X\theta = X^\top u$$
{{< /math >}}

This is the normal equation. When {{< math >}}$X^\top X${{< /math >}} is invertible (which it is in our case since {{< math >}}$X${{< /math >}} has full column rank), we can solve for {{< math >}}$\theta${{< /math >}}:

{{< math >}}
$$\theta = (X^\top X)^{-1}X^\top u$$
{{< /math >}}

This solution minimizes the squared error because:
1. The loss function is convex (it's quadratic and {{< math >}}$X^\top X${{< /math >}} is positive definite)
2. We found where its gradient is zero
3. The second derivative (Hessian) {{< math >}}$2X^\top X${{< /math >}} is positive definite

Therefore, this solution gives us the global minimum of the least squares problem.

#### Derive the Gradient Equation above for the Least Squares Loss Function

Let's walk through how we get:

{{< math >}}
$$\nabla_\theta L = -2X^\top u + 2X^\top X\theta$$
{{< /math >}}

##### Step 1: Write the Loss Function

The least squares loss is:

{{< math >}}
$$L(\theta) = \|u-X\theta\|^2 = (u-X\theta)^\top(u-X\theta)$$
{{< /math >}}

##### Step 2: Expand the Quadratic Form

Let's expand this expression:

{{< math >}}
$$L(\theta) = u^\top u - 2\theta^\top X^\top u + \theta^\top X^\top X\theta$$
{{< /math >}}

Here's how we get this expansion:

{{< math >}}
$$(u-X\theta)^\top(u-X\theta) = u^\top u - 2\theta^\top X^\top u + \theta^\top X^\top X\theta$$
{{< /math >}}

This follows from two key matrix properties:

1. {{< math >}}$(AB)^T = B^T A^T${{< /math >}}
2. {{< math >}}$u^\top X\theta = \theta^\top X^\top u${{< /math >}} (scalar equality)

##### Step 3: Take the Gradient with Respect to $\theta$

We differentiate term by term:

1. First term:
   {{< math >}}
   $$\nabla_\theta(u^\top u) = 0$$
   {{< /math >}}
   because it doesn't depend on {{< math >}}$\theta${{< /math >}}

2. Second term:
   {{< math >}}
   $$\nabla_\theta(-2\theta^\top X^\top u) = -2X^\top u$$
   {{< /math >}}
   because it's linear in {{< math >}}$\theta${{< /math >}}

3. Third term:
   {{< math >}}
   $$\nabla_\theta(\theta^\top X^\top X\theta) = 2X^\top X\theta \tag{quadratic grad.}$$
   {{< /math >}}
   This is a quadratic form, and its derivative is a standard result from matrix calculus.

Adding all parts together:

{{< math >}}
$$\nabla_\theta L = -2X^\top u + 2X^\top X\theta$$
{{< /math >}}

This gradient expression gives us the direction of steepest ascent of the loss function at any point {{< math >}}$\theta${{< /math >}}. Setting it to zero leads to the normal equations, which give us the optimal solution.

##### Note: Vector Dot Product Commutativity

Please note the identity

{{< math >}}
$$a^\top b = b^\top a$$
{{< /math >}}

holds true. This is a foundational property in linear algebra, and the key point is:

Both {{< math >}}$a^\top b${{< /math >}} and {{< math >}}$b^\top a${{< /math >}} are scalars (i.e., single numbers), and they are equal because they compute the same dot product. Let's break it down:

Let {{< math >}}$a = [a_1, a_2, \ldots, a_n]^\top${{< /math >}}, and {{< math >}}$b = [b_1, b_2, \ldots, b_n]^\top${{< /math >}}.

Then:

{{< math >}}
$$a^\top b = \sum_{i=1}^n a_i b_i$$
{{< /math >}}

{{< math >}}
$$b^\top a = \sum_{i=1}^n b_i a_i$$
{{< /math >}}

But since real number multiplication is commutative (i.e., {{< math >}}$a_i b_i = b_i a_i${{< /math >}}), the two sums are identical:

{{< math >}}
$$a^\top b = b^\top a$$
{{< /math >}}

##### Note on shapes and intuition:

- {{< math >}}$a^\top b${{< /math >}} is a {{< math >}}$1\times1${{< /math >}} scalar (row vector × column vector)
- {{< math >}}$b^\top a${{< /math >}} is also a {{< math >}}$1\times1${{< /math >}} scalar (same logic)
- This is just the dot product: it doesn't matter which order you take the dot product in. It's symmetric. This property is crucial in many derivations in linear algebra and optimization, including our least squares derivation where we used {{< math >}}$u^\top X\theta = \theta^\top X^\top u${{< /math >}}.

Since they're both just numbers, and equal by commutativity, the identity holds.


### Derive Quadratic Form Gradient in Equation (quadratic grad.)

Let's prove the key gradient identity with Matrix Calculus:

{{< math >}}
$$\nabla_\theta(\theta^\top A\theta) = (A + A^\top)\theta$$
{{< /math >}}

And its special case when {{< math >}}$A${{< /math >}} is symmetric:

{{< math >}}
$$\nabla_\theta(\theta^\top A\theta) = 2A\theta$$
{{< /math >}}

#### Setup

Given:
- {{< math >}}$\theta \in \mathbb{R}^p${{< /math >}} is a column vector
- {{< math >}}$A \in \mathbb{R}^{p \times p}${{< /math >}} is a constant matrix
- The scalar function is {{< math >}}$f(\theta) = \theta^\top A\theta = \sum_{i=1}^p \sum_{j=1}^p \theta_i A_{ij} \theta_j${{< /math >}}

#### Overview of Steps

**Step 1**: Write in Summation Form

{{< math >}}
$$f(\theta) = \sum_{i=1}^p \sum_{j=1}^p \theta_i A_{ij} \theta_j$$
{{< /math >}}

**Step 2**: Compute Partial Derivatives

The gradient's k-th component is:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \frac{\partial}{\partial \theta_k} \sum_{i=1}^p \sum_{j=1}^p \theta_i A_{ij} \theta_j$$
{{< /math >}}

**Step 3**: Use Product Rule

When differentiating:

{{< math >}}
$$\frac{\partial}{\partial \theta_k}(\theta_i A_{ij} \theta_j) = A_{ij}(\delta_{ik}\theta_j + \theta_i\delta_{jk})$$
{{< /math >}}

where {{< math >}}$\delta_{ik}${{< /math >}} is the Kronecker delta:

{{< math >}}
$$\delta_{ik} = \begin{cases} 1 & \text{if } i=k \\ 0 & \text{otherwise} \end{cases}$$
{{< /math >}}

**Step 4**: Sum Terms

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \sum_{j=1}^p A_{kj}\theta_j + \sum_{i=1}^p \theta_i A_{ik}$$
{{< /math >}}

**Step 5**: Express in Vector Form

This gives us:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = (A\theta)_k + (A^\top\theta)_k$$
{{< /math >}}

Therefore:

{{< math >}}
$$\nabla_\theta f = A\theta + A^\top\theta = (A + A^\top)\theta$$
{{< /math >}}

**Symmetric Matrix**: a special case 

When {{< math >}}$A = A^\top${{< /math >}}, we get:

{{< math >}}
$$\nabla_\theta(\theta^\top A\theta) = 2A\theta$$
{{< /math >}}

This is the form we use in least squares optimization where {{< math >}}$A = X^\top X${{< /math >}} is symmetric.


#### Detail Proof of Key Steps above

Let's prove that:

{{< math >}}
$$\frac{\partial}{\partial \theta_k} \left(\sum_{i=1}^p \sum_{j=1}^p \theta_i A_{ij} \theta_j\right) = \sum_{j=1}^p A_{kj}\theta_j + \sum_{i=1}^p \theta_i A_{ik}$$
{{< /math >}}

##### **Step-by-step Derivation**

Start with the function:

{{< math >}}
$$f(\theta) = \sum_{i=1}^p \sum_{j=1}^p \theta_i A_{ij} \theta_j$$
{{< /math >}}

Take partial derivative with respect to {{< math >}}$\theta_k${{< /math >}}, using product rule for all {{< math >}}$\theta_i${{< /math >}} and {{< math >}}$\theta_j${{< /math >}} terms:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \sum_{i=1}^p \sum_{j=1}^p A_{ij} \cdot \frac{\partial}{\partial \theta_k}(\theta_i \theta_j)$$
{{< /math >}}

Now:

{{< math >}}
$$\frac{\partial}{\partial \theta_k}(\theta_i \theta_j) = \delta_{ik}\theta_j + \theta_i\delta_{jk}$$
{{< /math >}}

where {{< math >}}$\delta_{ik}${{< /math >}} is the Kronecker delta:

{{< math >}}
$$\delta_{ik} = \begin{cases} 1 & \text{if } i=k \\ 0 & \text{otherwise} \end{cases}$$
{{< /math >}}

Therefore:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \sum_{i=1}^p \sum_{j=1}^p A_{ij}(\delta_{ik}\theta_j + \theta_i\delta_{jk})$$
{{< /math >}}

Split into two sums:

First term:
{{< math >}}
$$\sum_{i=1}^p \sum_{j=1}^p A_{ij}\delta_{ik}\theta_j = \sum_{j=1}^p A_{kj}\theta_j$$
{{< /math >}}

(since {{< math >}}$\delta_{ik}=1${{< /math >}} only when {{< math >}}$i=k${{< /math >}})

Second term:
{{< math >}}
$$\sum_{i=1}^p \sum_{j=1}^p A_{ij}\theta_i\delta_{jk} = \sum_{i=1}^p A_{ik}\theta_i$$
{{< /math >}}

(since {{< math >}}$\delta_{jk}=1${{< /math >}} only when {{< math >}}$j=k${{< /math >}})

##### **Interim Result**

Combining both terms:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \sum_{j=1}^p A_{kj}\theta_j + \sum_{i=1}^p A_{ik}\theta_i$$
{{< /math >}}

In vector form, this means:

{{< math >}}
$$\nabla_\theta(\theta^\top A\theta) = A\theta + A^\top\theta = (A + A^\top)\theta$$
{{< /math >}}

##### **From Component-wise to Vector Form of Gradient**

For each coordinate {{< math >}}$k \in \{1,\ldots,p\}${{< /math >}}, we showed:

{{< math >}}
$$\frac{\partial f}{\partial \theta_k} = \underbrace{\sum_{j=1}^p A_{kj}\theta_j}_{(A\theta)_k} + \underbrace{\sum_{i=1}^p A_{ik}\theta_i}_{(A^\top\theta)_k}$$
{{< /math >}}

So each element of the gradient vector is:

{{< math >}}
$$(\nabla_\theta f)_k = (A\theta)_k + (A^\top\theta)_k$$
{{< /math >}}

##### **Stacking into Vector Form**

The full gradient vector is:

{{< math >}}
$$\nabla_\theta f = \begin{bmatrix} 
(A\theta)_1 + (A^\top\theta)_1 \\
(A\theta)_2 + (A^\top\theta)_2 \\
\vdots \\
(A\theta)_p + (A^\top\theta)_p
\end{bmatrix} = A\theta + A^\top\theta$$
{{< /math >}}
The notation {{< math >}}$(A\theta)_k${{< /math >}} which refers to the k-th component of the matrix-vector product {{< math >}}$A\theta${{< /math >}}. This follows because matrix-vector multiplication simply stacks the components:

{{< math >}}
$$(A\theta)_k = \sum_j A_{kj}\theta_j$$
{{< /math >}}

{{< math >}}
$$A\theta = \begin{bmatrix} 
\sum_{j=1}^p A_{1j}\theta_j \\
\sum_{j=1}^p A_{2j}\theta_j \\
\vdots \\
\sum_{j=1}^p A_{pj}\theta_j
\end{bmatrix}$$
{{< /math >}}
Similarly, for the transpose:
{{< math >}}
$$(A^\top\theta)_k = \sum_i A_{ik}\theta_i$$
{{< /math >}}

##### **Final Result**

Therefore:

{{< math >}}
$$\nabla_\theta(\theta^\top A\theta) = A\theta + A^\top\theta = (A + A^\top)\theta$$
{{< /math >}}



### From Normal Equations to Standard Regression Form (3.6 in Section D)

Let's derive how the normal equations solution:

{{< math >}}
$$\theta = \frac{1}{nS_2-S_1^2}\begin{bmatrix} n & -S_1 \\ -S_1 & S_2 \end{bmatrix}\begin{bmatrix} SU \\ U \end{bmatrix}$$
{{< /math >}}

relates to the standard regression form:

{{< math >}}
$$\gamma_0 = \frac{\text{Cov}(u,s)}{\text{Var}(s)} = \frac{n\sum s_i u_i - \sum s_i \sum u_i}{n\sum s_i^2 - (\sum s_i)^2}$$
{{< /math >}}

with offset:

{{< math >}}
$$o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}

#### Given Matrices and Vectors

Let {{< math >}}$s=[s_1,\ldots,s_n]^T${{< /math >}}, {{< math >}}$u=[u_1,\ldots,u_n]^T${{< /math >}}

Model:
{{< math >}}
$$u_i = \gamma_0 s_i + o + \varepsilon_i \text{ or } u = \gamma_0 s + o\cdot 1 + \varepsilon$$
{{< /math >}}

#### Matrix Form

Design matrix:
{{< math >}}
$$X = [s \quad 1] \in \mathbb{R}^{n\times 2}$$
{{< /math >}}

Parameters:
{{< math >}}
$$\theta = \begin{bmatrix} \gamma_0 \\ o \end{bmatrix}$$
{{< /math >}}

#### Computing Products

{{< math >}}
$$X^T X = \begin{bmatrix} \sum s_i^2 & \sum s_i \\ \sum s_i & n \end{bmatrix} = \begin{bmatrix} S_2 & S_1 \\ S_1 & n \end{bmatrix}$$
{{< /math >}}

{{< math >}}
$$X^T u = \begin{bmatrix} \sum s_i u_i \\ \sum u_i \end{bmatrix} = \begin{bmatrix} SU \\ U \end{bmatrix}$$
{{< /math >}}

#### Solution Components

First row (slope):
{{< math >}}
$$\gamma_0 = \frac{1}{nS_2-S_1^2}(n\cdot SU - S_1\cdot U)$$
{{< /math >}}

Second row (intercept):
{{< math >}}
$$o = \frac{1}{nS_2-S_1^2}(-S_1\cdot SU + S_2\cdot U)$$
{{< /math >}}

#### Statistical Interpretation

Define:
{{< math >}}
$$\bar{s} = \frac{1}{n}\sum s_i \quad \bar{u} = \frac{1}{n}\sum u_i$$
{{< /math >}}

{{< math >}}
$$\text{Var}(s) = \frac{1}{n}\sum(s_i-\bar{s})^2 = \frac{1}{n}S_2 - \bar{s}^2$$
{{< /math >}}

{{< math >}}
$$\text{Cov}(u,s) = \frac{1}{n}\sum(u_i-\bar{u})(s_i-\bar{s}) = \frac{1}{n}SU - \bar{u}\bar{s}$$
{{< /math >}}

#### Final Result

This gives us:
{{< math >}}
$$\gamma_0 = \frac{\text{Cov}(u,s)}{\text{Var}(s)} \text{ and } o = \bar{u} - \gamma_0\bar{s}$$
{{< /math >}}