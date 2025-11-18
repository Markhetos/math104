---
title: "Vector Spaces & Subspaces"
permalink: /notes/vector-spaces/
---

# Math 104: Applied Matrix Theory
## Vector Spaces & Subspaces
**Date:** Nov 17, 2025  
**Instructor:** Dr. Emmanuel Candès  
**References:** Lecture 4 | Laub Ch. 2 | HW1

---

## 1. PRELIMINARY REVIEW

### 1.1 Vectors and Matrices

**Notation:**
- Let $\mathbb{R}^n$ be the set of $n$-tuples of real numbers represented as column vectors
- If $\mathbf{x} \in \mathbb{R}^n$, then $\mathbf{x}$ is a **vertical vector** (column vector)

$$\mathbf{x} = \begin{bmatrix} 
x_1 \\ 
x_2 \\ 
x_3 \\
\vdots \\ 
x_n 
\end{bmatrix}$$

- $\mathbf{x}^T$ is a **row vector** (transpose): $\mathbf{x}^T = (x_1, x_2, \ldots, x_n)$

**Key dimensional results:**
- $\mathbf{x}^T\mathbf{y}$ is a **scalar** (inner product)
- $\mathbf{xy}^T$ is an $n \times n$ **matrix** (outer product)

**Fields:**
- Real numbers: $\mathbf{x} \in \mathbb{R}^n$
- Complex numbers: $\mathbf{x} \in \mathbb{C}^n$

### 1.2 Matrix Shapes

**Common matrix types for** $\mathbf{A} = [a_{ij}] \in \mathbb{R}^{m \times n}$:

1. **Diagonal matrix:** $a_{ij} = 0$ for $i \neq j$ (all off-diagonal elements are zero)

$$\mathbf{D} = \begin{bmatrix} 
d_1 & 0 & 0 & \cdots & 0 \\ 
0 & d_2 & 0 & \cdots & 0 \\ 
0 & 0 & d_3 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\ 
0 & 0 & 0 & \cdots & d_n 
\end{bmatrix}$$

2. **Upper triangular matrix:** $a_{ij} = 0$ for $i > j$ (all elements below diagonal are zero)

3. **Lower triangular matrix:** $a_{ij} = 0$ for $i < j$ (all elements above diagonal are zero)

**Matrix elements:**
- **Main diagonal:** elements $a_{ii}$ where $i = j$
- **Anti-diagonal:** elements $a_{ij}$ where $i + j = n + 1$
- **Blocks:** submatrices within a larger matrix

### 1.3 Determinant

**Definition:** For $\mathbf{A} \in \mathbb{R}^{n \times n}$, we write $\det(\mathbf{A})$ for the determinant of $\mathbf{A}$.

**Geometric interpretation:** The determinant represents the signed volume scaling factor of the linear transformation.

**Key properties:**

1. If $\mathbf{A}$ has a zero row or two equal rows, then $\det(\mathbf{A}) = 0$
2. If $\mathbf{A}$ has a zero column or two equal columns, then $\det(\mathbf{A}) = 0$
3. Interchanging two rows (or columns) changes the sign: $\det(\mathbf{A}) \to -\det(\mathbf{A})$
4. For scalar $\alpha \in \mathbb{R}$: $\det(\alpha\mathbf{A}) = \alpha^n \det(\mathbf{A})$ (for $n \times n$ matrix)
5. $\det(\mathbf{A}^T) = \det(\mathbf{A})$
6. If $\mathbf{A}$ is diagonal or triangular, then $\det(\mathbf{A}) = \prod_{i=1}^n a_{ii}$ (product of diagonal elements)
7. If $\mathbf{A}, \mathbf{B} \in \mathbb{R}^{n \times n}$, then $\det(\mathbf{AB}) = \det(\mathbf{A})\det(\mathbf{B})$
8. $\mathbf{A}$ is invertible if and only if $\det(\mathbf{A}) \neq 0$, and then $\det(\mathbf{A}^{-1}) = \frac{1}{\det(\mathbf{A})}$

### 1.4 Inner Products and Orthogonality

**Inner product:** For $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$,

$$\langle \mathbf{x}, \mathbf{y} \rangle := \mathbf{x}^T\mathbf{y} = \sum_{i=1}^n x_i y_i$$

**Orthogonality:** Vectors $\mathbf{v}$ and $\mathbf{w}$ are **orthogonal** (written $\mathbf{v} \perp \mathbf{w}$) if and only if

$$\mathbf{v} \perp \mathbf{w} \iff \mathbf{v} \cdot \mathbf{w} = 0 \iff \langle \mathbf{v}, \mathbf{w} \rangle = 0$$

---

## 2. VECTOR SPACES

### 2.1 Definition

**Definition:** A **vector space** over a field $\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$ is a set $V$ with two operations:
- Addition: $+ : V \times V \to V$
- Scalar multiplication: $\times : \mathbb{F} \times V \to V$

satisfying the following **axioms** for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and $\alpha, \beta \in \mathbb{F}$:

1. **Commutativity:** $\mathbf{v} + \mathbf{w} = \mathbf{w} + \mathbf{v}$
2. **Associativity (addition):** $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$
3. **Zero vector:** There exists $\mathbf{0} \in V$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$
4. **Additive inverse:** For each $\mathbf{v} \in V$, there exists $-\mathbf{v} \in V$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$
5. **Associativity (scalar mult.):** $(\alpha \cdot \beta)\mathbf{v} = \alpha(\beta \cdot \mathbf{v})$
6. **Distributivity (scalar):** $(\alpha + \beta)\mathbf{v} = \alpha\mathbf{v} + \beta\mathbf{v}$
7. **Distributivity (vector):** $\alpha(\mathbf{v} + \mathbf{w}) = \alpha\mathbf{v} + \alpha\mathbf{w}$
8. **Identity:** $1\mathbf{v} = \mathbf{v}$

**Example 2:** $\mathbb{R}^{m \times n}$ (the set of all $m \times n$ matrices)

- **Addition:**
  ```math
  (\mathbf{A} + \mathbf{B})_{ij} = A_{ij} + B_{ij}
  ```
- **Scalar multiplication:**
  ```math
  (\alpha \mathbf{A})_{ij} = \alpha A_{ij}
  ```

**Example 3:** $P_n$ = polynomials of degree at most $n$

$$P_n = \{p : p(t) = c_0 + c_1 t + c_2 t^2 + \cdots + c_n t^n\}$$

- Addition: $(p + q)(t) = p(t) + q(t)$
- Scalar multiplication: $(\alpha p)(t) = \alpha p(t)$

**Representation via coefficients:**

$$p(t) = c_0 + c_1 t + c_2 t^2 + \cdots + c_n t^n \longleftrightarrow \begin{bmatrix} 
c_0 \\ 
c_1 \\ 
c_2 \\
\vdots \\ 
c_n 
\end{bmatrix} \in \mathbb{R}^{n+1}$$

This shows $P_n$ is **isomorphic** to $\mathbb{R}^{n+1}$.

---

## 3. SUBSPACES

### 3.1 Definition

**Definition:** $W \subset V$ is a **subspace** of vector space $V$ if and only if $\alpha\mathbf{w}_1 + \beta\mathbf{w}_2 \in W$ for all $\alpha, \beta \in \mathbb{F}$ and $\mathbf{w}_1, \mathbf{w}_2 \in W$.

**Equivalently:** A subset $W \subseteq V$ is a subspace if:
1. $\mathbf{0} \in W$ (contains zero vector)
2. $W$ is **closed under addition**: if $\mathbf{u}, \mathbf{v} \in W$, then $\mathbf{u} + \mathbf{v} \in W$
3. $W$ is **closed under scalar multiplication**: if $\mathbf{v} \in W$ and $\alpha \in \mathbb{F}$, then $\alpha\mathbf{v} \in W$

### 3.2 Examples

**Example 1:** Any line through the origin in $\mathbb{R}^2$ is a subspace

**Example 2:** Any plane through the origin in $\mathbb{R}^3$ is a subspace

**Example 3:** Let $V = \mathbb{R}^{n \times n}$ and $W = \{\mathbf{A} \in \mathbb{R}^{n \times n} : \mathbf{A} \text{ is symmetric}\}$. Then $W$ is a subspace of $V$.

**Proof:** 
- $\mathbf{0}$ is symmetric, so $\mathbf{0} \in W$
- If $\mathbf{A}, \mathbf{B}$ symmetric, then $(\mathbf{A} + \mathbf{B})^T = \mathbf{A}^T + \mathbf{B}^T = \mathbf{A} + \mathbf{B}$, so $\mathbf{A} + \mathbf{B} \in W$
- If $\mathbf{A}$ symmetric, then $(\alpha\mathbf{A})^T = \alpha\mathbf{A}^T = \alpha\mathbf{A}$, so $\alpha\mathbf{A} \in W$ ∎

**Example 4 (Non-subspace):** A line NOT through the origin is NOT a subspace (doesn't contain $\mathbf{0}$)

### 3.3 Equality of Subspaces

**Proposition:** If $R$ and $S$ are vector spaces (or subspaces), then $R = S$ if and only if $R \subseteq S$ **and** $S \subseteq R$.

**Proof strategy:** To prove $R = S$:
1. Take arbitrary $\mathbf{r} \in R$ and show $\mathbf{r} \in S$ (proves $R \subseteq S$)
2. Take arbitrary $\mathbf{s} \in S$ and show $\mathbf{s} \in R$ (proves $S \subseteq R$)

---

## 4. LINEAR INDEPENDENCE AND SPAN

### 4.1 Span

**Definition:** The **span** of vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ in $\mathbb{R}^n$ is the set of all linear combinations:

$$\text{span}(\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k) = \{ \mathbf{x} : \mathbf{x} = \sum_{i=1}^k \alpha_i\mathbf{v}_i, \ \alpha_i \in \mathbb{R} \}$$

**Key fact:** The span of any collection of vectors is a subspace.

**Example:** 

$$
\text{span}\({
\begin{bmatrix}
1 \\
0 \\
0 
\end{bmatrix},
\begin{bmatrix}
0 \\
1 \\
0 \end{bmatrix},
\begin{bmatrix}
0 \\
0 \\
1
\end{bmatrix}
\}) = \mathbb{R}^3
$$

### 4.2 Linear Independence

**Definition:** Vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ are **linearly dependent** if there exist scalars $\alpha_1, \alpha_2, \ldots, \alpha_k \in \mathbb{F}$ **not all equal to zero** such that

$$\sum_{i=1}^k \alpha_i\mathbf{v}_i = \mathbf{0}$$

**Definition:** Vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ are **linearly independent** if the only way to have

$$\sum_{i=1}^k \alpha_i\mathbf{v}_i = \mathbf{0}$$

is when $\alpha_1 = \alpha_2 = \cdots = \alpha_k = 0$.

**Geometric intuition:**
- **Linearly dependent:** At least one vector is redundant (can be expressed as a combination of others)
- **Linearly independent:** All vectors are "essential" (none can be expressed as a combination of others)

**Key facts:**
1. Any set containing the zero vector is linearly dependent
2. A set with duplicate vectors is linearly dependent
3. Eigenvectors corresponding to distinct eigenvalues are linearly independent

**Matrix interpretation:** For $\mathbf{V} = [\mathbf{v}_1 \ \mathbf{v}_2 \ \cdots \ \mathbf{v}_k] \in \mathbb{R}^{m \times k}$,
- Columns are linearly **dependent** if and only if there exists $\mathbf{x} \neq \mathbf{0}$ such that $\mathbf{Vx} = \mathbf{0}$
- Columns are linearly **independent** if and only if $\mathbf{Vx} = \mathbf{0}$ implies $\mathbf{x} = \mathbf{0}$

### 4.3 Basis and Dimension

**Definition:** A set of vectors $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n\}$ is a **basis** for $V$ if and only if:
1. The vectors are **linearly independent**
2. $\text{span}(\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n) = V$

**Examples:**

1. **Standard basis for** $\mathbb{R}^n$:

$$\mathbf{e}_1 = \begin{bmatrix} 
1 \\ 
0 \\ 
0 \\
\vdots \\ 
0 
\end{bmatrix}, \quad \mathbf{e}_2 = \begin{bmatrix} 
0 \\ 
1 \\ 
0 \\
\vdots \\ 
0 
\end{bmatrix}, \quad \ldots, \quad \mathbf{e}_n = \begin{bmatrix} 
0 \\ 
0 \\ 
0 \\
\vdots \\ 
1 
\end{bmatrix}$$

2. **Basis for** $P_4$ (polynomials of degree $\leq 4$): $\{1, t, t^2, t^3, t^4\}$

**Uniqueness of representation:** If $\{\mathbf{b}_1, \mathbf{b}_2, \ldots, \mathbf{b}_n\}$ is a basis for $V$, then every $\mathbf{v} \in V$ has a **unique** representation:

$$\mathbf{v} = \sum_{i=1}^n x_i\mathbf{b}_i = \mathbf{Bx}$$

where $\mathbf{B} = [\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_n]$ and $\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$.

**Why unique?** Suppose $\mathbf{v} = \sum x_i\mathbf{b}_i = \sum y_i\mathbf{b}_i$. Then $\sum(x_i - y_i)\mathbf{b}_i = \mathbf{0}$. By linear independence, $x_i - y_i = 0$ for all $i$, so $x_i = y_i$.

### 4.4 Dimension

**Theorem:** All bases of a subspace have the same number of elements.

**Definition:** If a basis for $V$ has $n$ elements, we say $V$ is **n-dimensional** and write $\dim(V) = n$.

**Examples:**

| Vector Space V | Dimension |
|----------------|-----------|
| $\mathbb{R}^n$ | $n$ |
| $\mathbb{R}^{m \times n}$ | $mn$ |
| Symmetric matrices in $\mathbb{R}^{n \times n}$ | $n(n+1)/2$ |
| Upper triangular matrices in $\mathbb{R}^{n \times n}$ | $n(n+1)/2$ |
| $P_n$ (polynomials of degree $\leq n$) | $n + 1$ |

---

## 5. FOUR FUNDAMENTAL SUBSPACES

For matrix $\mathbf{A} \in \mathbb{R}^{m \times n}$, there are four important associated subspaces:

### 5.1 Column Space

**Definition:** The **column space** (or **range**) of $\mathbf{A}$, denoted $\mathcal{C}(\mathbf{A})$ or $\mathcal{R}(\mathbf{A})$, is the span of the column vectors:

$$\mathcal{C}(\mathbf{A}) = \text{span}(\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n)$$

where $\mathbf{a}_j$ is the $j$-th column of $\mathbf{A}$.

**Interpretation:** $\mathcal{C}(\mathbf{A})$ represents **all possible outputs** of the transformation $\mathbf{Ax}$.

**Key fact:** $\mathcal{C}(\mathbf{A}) \subseteq \mathbb{R}^m$

### 5.2 Row Space

**Definition:** The **row space** of $\mathbf{A}$ is the span of the row vectors, which equals the column space of $\mathbf{A}^T$:

$$\mathcal{R}(\mathbf{A}) = \mathcal{C}(\mathbf{A}^T) = \text{span}(\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_m)$$

where $\mathbf{w}_i$ are the columns of $\mathbf{A}^T$ (i.e., rows of $\mathbf{A}$).

**Key fact:** $\mathcal{R}(\mathbf{A}) \subseteq \mathbb{R}^n$

### 5.3 Null Space

**Definition:** The **null space** (or **kernel**) of $\mathbf{A}$, denoted $\mathcal{N}(\mathbf{A})$, is the set of all vectors $\mathbf{x}$ such that $\mathbf{Ax} = \mathbf{0}$:

$$\mathcal{N}(\mathbf{A}) = \{\mathbf{x} \in \mathbb{R}^n : \mathbf{Ax} = \mathbf{0}\}$$

**Interpretation:** $\mathcal{N}(\mathbf{A})$ represents **inputs that are mapped to the zero vector**.

**Key fact:** $\mathcal{N}(\mathbf{A}) \subseteq \mathbb{R}^n$

### 5.4 Left Null Space

**Definition:** The **left null space** of $\mathbf{A}$ is the null space of $\mathbf{A}^T$:

$$\mathcal{N}(\mathbf{A}^T) = \{\mathbf{y} \in \mathbb{R}^m : \mathbf{A}^T\mathbf{y} = \mathbf{0}\}$$

Equivalently, $\mathbf{y}^T\mathbf{A} = \mathbf{0}^T$ (hence "left" null space).

**Key fact:** $\mathcal{N}(\mathbf{A}^T) \subseteq \mathbb{R}^m$

### 5.5 Orthogonality Relations

**Fundamental Theorem:**
1. $\mathcal{N}(\mathbf{A}) \perp \mathcal{R}(\mathbf{A})$ (null space perpendicular to row space)
2. $\mathcal{N}(\mathbf{A}^T) \perp \mathcal{C}(\mathbf{A})$ (left null space perpendicular to column space)

**Direct sum decompositions:**
1. $\mathbb{R}^n = \mathcal{N}(\mathbf{A}) \oplus \mathcal{R}(\mathbf{A})$ (null space direct sum row space)
2. $\mathbb{R}^m = \mathcal{N}(\mathbf{A}^T) \oplus \mathcal{C}(\mathbf{A})$ (left null space direct sum column space)

**Geometric picture:**
```
          ℝⁿ                           ℝᵐ
    ┌──────────────┐             ┌──────────────┐
    │              │             │              │
    │  Row Space   │    A        │ Column Space │
    │     R(A)     │ ────────>   │     C(A      │
    │              │             │              │
    ├──────────────┤             ├──────────────┤
    │  Null Space  │             │ Left Null    │
    │     N(A)     │             │ Space N(Aᵀ)  │
    └──────────────┘             └──────────────┘
         ⊥                             ⊥
```

---

## 6. RANK-NULLITY THEOREM

**Theorem (Rank-Nullity):** For $\mathbf{A} \in \mathbb{R}^{m \times n}$,

$$\text{rank}(\mathbf{A}) + \text{nullity}(\mathbf{A}) = n$$

or equivalently,

$$\dim(\mathcal{C}(\mathbf{A})) + \dim(\mathcal{N}(\mathbf{A})) = \text{number of columns of } \mathbf{A}$$

**Corollary:** 

$$\dim(\mathcal{C}(\mathbf{A})) = \dim(\mathcal{R}(\mathbf{A})) = \text{rank}(\mathbf{A})$$

That is, the **column space** and **row space** have the **same dimension**! (This common dimension is called the **rank**.)

**Proof sketch:**
1. Let $\{\mathbf{v}_1, \ldots, \mathbf{v}_r\}$ be a basis for the row space
2. Let $\{\mathbf{v}_{r+1}, \ldots, \mathbf{v}_n\}$ be a basis for the null space
3. Column space is $\text{span}\{\mathbf{Av}_1, \ldots, \mathbf{Av}_r\}$
4. Show $\{\mathbf{Av}_1, \ldots, \mathbf{Av}_r\}$ are linearly independent
5. Therefore $\dim(\mathcal{C}(\mathbf{A})) = r = \dim(\mathcal{R}(\mathbf{A}))$ ∎

---

## 7. KEY TAKEAWAYS

**Vector Spaces:**
- Abstract sets with addition and scalar multiplication satisfying 8 axioms
- Examples: $\mathbb{R}^n$, matrices, polynomials, function spaces

**Subspaces:**
- Subsets closed under addition and scalar multiplication
- Must contain zero vector
- Geometric: lines/planes through origin

**Linear Independence:**
- Vectors are independent if none is a combination of others
- Matrix columns independent if and only if $\mathbf{Ax} = \mathbf{0}$ implies $\mathbf{x} = \mathbf{0}$

**Basis & Dimension:**
- Basis = linearly independent spanning set
- All bases have same size = dimension
- Representation in basis is unique

**Four Fundamental Subspaces:**
- $\mathcal{C}(\mathbf{A}) \subseteq \mathbb{R}^m$ (column space)
- $\mathcal{R}(\mathbf{A}) \subseteq \mathbb{R}^n$ (row space)
- $\mathcal{N}(\mathbf{A}) \subseteq \mathbb{R}^n$ (null space)
- $\mathcal{N}(\mathbf{A}^T) \subseteq \mathbb{R}^m$ (left null space)

**Orthogonality:**
- $\mathcal{N}(\mathbf{A}) \perp \mathcal{R}(\mathbf{A})$
- $\mathcal{N}(\mathbf{A}^T) \perp \mathcal{C}(\mathbf{A})$

**Rank-Nullity:**
- $\text{rank}(\mathbf{A}) + \text{nullity}(\mathbf{A}) = n$
- Column rank = row rank

---

## 8. COMMON MISTAKES TO AVOID

1. ❌ **Confusing row space with column space locations**
   - ✅ Column space lives in $\mathbb{R}^m$, row space lives in $\mathbb{R}^n$

2. ❌ **Forgetting zero vector in subspace**
   - ✅ Every subspace must contain $\mathbf{0}$

3. ❌ **Mixing up "linearly dependent" condition**
   - ✅ Need "not all $\alpha_i = 0$" for dependence, not "$\alpha_i \neq 0$"

4. ❌ **Thinking any n vectors in** $\mathbb{R}^n$ **form a basis**
   - ✅ Need linear independence, not just correct count

5. ❌ **Confusing dimension with number of elements**
   - ✅ Dimension = number of basis vectors, not total elements in space

---

## 9. PRACTICE PROBLEMS

**Problem 1:** Determine if $W = \{(x, y, z) \in \mathbb{R}^3 : x + 2y - z = 0\}$ is a subspace of $\mathbb{R}^3$.

**Problem 2:** Find a basis for the column space and null space of

$$\mathbf{A} = \begin{bmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{bmatrix}$$

**Problem 3:** Prove that $\dim(\mathbb{R}^{n \times n}_{\text{symmetric}}) = \frac{n(n+1)}{2}$.

**Problem 4:** If $\mathbf{A}$ is $5 \times 8$ with rank 3, what is $\dim(\mathcal{N}(\mathbf{A}))$?

---

**END OF NOTES**
