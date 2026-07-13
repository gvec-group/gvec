# Fast implementation of 2D DFT

The two-dimensional direct fourier transform is used in GVEC, and we use dimension-by-dimension operators storing 1D matrices only for a fast transform, and making use of matrix-matrix multiplication.

## Representation of variables in Fourier modes

A variable in the poloidal and toroidal angle periodic directions $(\thet, \zeta)$ is represented as a truncated real-valued Fourier series of either only cosine, only sine, or both. The mode number in poloidal direction is $0 \leq m \leq \mmax$ and in toroidal direction $-\nmax \leq n \leq \nmax$. It is possible to have a periodicity of the toroidal direction, the 'number of field periods' $\nfp$, such that $\zeta \in [0, 2\pi / \nfp]$.

The full 2D Fourier representation is split into the following parts:

$$
\begin{align}
x(\thet, \zeta) = x^c_{00}
  &+ \sum_{n=1}^{\nmax} x^c_{0n} \cos(0 \cdot \thet - n \nfp \zeta)
   + \sum_{m=1}^{\mmax} \sum_{n=-\nmax}^{\nmax} x^c_{mn} \cos(m\thet - n \nfp \zeta) \\
  &+ \sum_{n=1}^{\nmax} x^s_{0n} \sin(0 \cdot \thet - n \zeta)
   + \sum_{m=1}^{\mmax} \sum_{n=-\nmax}^{\nmax} x^s_{mn} \sin(m\thet - n \nfp \zeta)
\end{align}
$$

with the coefficients of cosine $x^c_{mn}$ and sine $x^s_{mn}$. One can choose to use the full series, with or without the zero mode, or only the cosine series with or without the zero mode, or only the sine series.

Depending on this choice, the 2D matrix of $(m, n)$ coefficients would not be filled. Therefore, only the non-zero coefficients are stored as a contiguous one-dimensional array, along with an array of same size containing the associated $m$ and $n \nfp$ mode numbers. Since only one single 1D array is used independent of the different choices, the start and end positions of the sine and cosine coefficients are also stored.

## Transform to real space

A sampling in real space is necessary, since we need to take integrals of nonlinear products of these variables (and of their derivatives in $\thet$ and $\zeta$). For periodic functions, the trapezoidal rule on an equidistant set of points is known to be spectrally accurate. We need to choose the same point positions for evaluating cosine and sine series,

$$
\begin{align}
\thet_i &= 2\pi \left(i - \tfrac{1}{2}\right) / \mnyq, \quad i = 1, \ldots, \mnyq \\
\zeta_j &= \frac{2\pi}{\nfp} \left(j - \tfrac{1}{2}\right) / \nnyq, \quad j = 1, \ldots, \nnyq
\end{align}
$$

with the number of points $(\mnyq, \nnyq)$. The number of integration points is typically chosen to be a factor 4 of the maximum mode number, due to the nonlinearity of the integral. (To guarantee the orthogonality of the product of two Fourier basis functions with maximum mode number $N$, the minimum number of integration points is $2N+1$.)

The transform then reads

$$
\begin{align}
x(\thet_i, \zeta_j) &= \sum_m \sum_n x^c_{mn} \cos(m\thet_i - n \nfp \zeta_j)
+ \sum_m \sum_n x^s_{mn} \sin(m\thet_i - n \nfp \zeta_j)
\end{align}
$$

where the sum over mode indices is not specified yet. This operation can be made more efficient by splitting the sum into the two directions. First,

$$
\begin{align}
x_{ij} = x(\thet_i, \zeta_j) =
  &\sum_m \sum_n \cos(m\thet_i)\, x^c_{mn} \cos(n \nfp \zeta_j) + \sin(m\thet_i)\, x^c_{mn} \sin(n \nfp \zeta_j) \\
+ &\sum_m \sum_n \sin(m\thet_i)\, x^s_{mn} \cos(n \nfp \zeta_j) - \cos(m\thet_i)\, x^s_{mn} \sin(n \nfp \zeta_j)
\end{align}
$$

which can be grouped as

$$
\begin{align}
x_{ij}
&= \sum_n \left(
    \sum_m \cos(m\thet_i)\, x^c_{mn} + \sin(m\thet_i)\, x^s_{mn}
  \right) \cos(n \nfp \zeta_j) \\
&\quad + \left(
    \sum_m \sin(m\thet_i)\, x^c_{mn} - \cos(m\thet_i)\, x^s_{mn}
  \right) \sin(n \nfp \zeta_j) \\
&= \sum_n \left(
    \sum_m
    \begin{pmatrix} \cos(m\thet_i) & \sin(m\thet_i) \\ \sin(m\thet_i) & -\cos(m\thet_i) \end{pmatrix}
    \cdot
    \begin{pmatrix} x^c_{mn} \\ x^s_{mn} \end{pmatrix}
  \right)
  \cdot
  \begin{pmatrix} \cos(n \nfp \zeta_j) \\ \sin(n \nfp \zeta_j) \end{pmatrix}^\top
\end{align}
$$

Since the result for all points $(i,j)$ is a matrix of size $[\mnyq \times \nnyq]$, one can write the operation in two consecutive matrix-matrix multiplications:

1. The Fourier coefficients $x^c_{mn}$ and $x^s_{mn}$ **must be copied into a matrix** of size $[\mtotal \times \ntotal]$, with $\mtotal = f(\mmax+1)$, $\ntotal = (2\nmax+1)$, where $f=1$ if sine or cosine series only, and $f=2$ if the full series. The few non-existing coefficients are set to zero.
2. The first matrix-matrix product has a left matrix of size $[2\mnyq \times \mtotal]$. This matrix is precomputed.
3. The result is a matrix of size $[2\mnyq \times \ntotal]$. It can be recast into a matrix of size $[\mnyq \times 2\ntotal]$, without any additional copy. (The result can be seen as a 3D array of size $[\mnyq \times 2 \times \ntotal]$, allowing both matrix shapes.)
4. The second matrix-matrix product has a right matrix of size $[2\ntotal \times \nnyq]$. This matrix is also precomputed.
5. The result is a 2D matrix of size $[\mnyq \times \nnyq]$.

Note that for the computation of a derivative in $\thet$ or $\zeta$ direction, the corresponding derivative matrices are precomputed and used instead:

$$
\frac{\partial}{\partial\thet}
\begin{pmatrix} \cos(m\thet) & \sin(m\thet) \\ \sin(m\thet) & -\cos(m\thet) \end{pmatrix}_{\!\thet_i}
\qquad
\frac{\partial}{\partial\zeta}
\begin{pmatrix} \cos(n \nfp \zeta) \\ \sin(n \nfp \zeta) \end{pmatrix}_{\!\zeta_j}
$$

The total number of multiplications for the two matrix-matrix products is $(2\mnyq\, \mtotal\, \ntotal) + (\mnyq\, 2\ntotal\, \nnyq) = 2\mnyq\, \ntotal(\mtotal + \nnyq)$.

## Transform to Fourier space

The integrand $g(\thet, \zeta)$ is projected onto the Fourier coefficients via integration. Using the orthogonality of the Fourier basis, one can write

$$
\begin{align}
g^c_{mn} &= \int_0^{2\pi} \int_0^{2\pi} g(\thet, \zeta)\, \cos(m\thet - n\zeta)\, d\thet\, d\zeta, \\
g^s_{mn} &= \int_0^{2\pi} \int_0^{2\pi} g(\thet, \zeta)\, \sin(m\thet - n\zeta)\, d\thet\, d\zeta
\end{align}
$$

Note that the normalization factor is not included here.

The integral is approximated by the sum over the integration points with a constant integration weight (not included here). With the integrand evaluated at the points $g_{ij} = g(\thet_i, \zeta_j)$, we get

$$
\begin{align}
g^c_{mn}
  &= \sum_{i=1}^{\mnyq} \sum_{j=1}^{\nnyq} g_{ij}\, \cos(m\thet_i - n\zeta_j) \\
  &= \sum_{i=1}^{\mnyq} \sum_{j=1}^{\nnyq}
     \cos(m\thet_i)\, g_{ij}\, \cos(n\zeta_j) + \sin(m\thet_i)\, g_{ij}\, \sin(n\zeta_j) \\
  &= \sum_{i=1}^{\mnyq}
     \begin{pmatrix} \cos(m\thet_i) \\ \sin(m\thet_i) \end{pmatrix}^\top
     \cdot \left(
       \sum_{j=1}^{\nnyq} g_{ij}
       \begin{pmatrix} \cos(n\zeta_j) \\ \sin(n\zeta_j) \end{pmatrix}
     \right)
\end{align}
$$

or for the generalized full series

$$
\begin{pmatrix} g^c_{mn} \\ g^s_{mn} \end{pmatrix}
= \sum_{i=1}^{\mnyq}
  \begin{pmatrix} \cos(m\thet_i) & \sin(m\thet_i) \\ \sin(m\thet_i) & -\cos(m\thet_i) \end{pmatrix}^\top
  \cdot \left(
    \sum_{j=1}^{\nnyq} g_{ij}
    \begin{pmatrix} \cos(n\zeta_j) \\ \sin(n\zeta_j) \end{pmatrix}
  \right)
$$

Using the abbreviations $\mtotal = f(\mmax+1)$, $\ntotal = (2\nmax+1)$, where $f=1$ if sine or cosine series only, and $f=2$ if the full series, the two matrix-matrix products are:

1. For the first matrix-matrix product, the left matrix holds $g_{ij}$ and has size $[\mnyq \times \nnyq]$; the right matrix is of size $[\nnyq \times 2\ntotal]$. This matrix is precomputed.
2. The result is a matrix of size $[\mnyq \times 2\ntotal]$, which can be recast into a matrix of size $[2\mnyq \times \ntotal]$.
3. The second matrix-matrix product has a left matrix of size $[\mtotal \times 2\mnyq]$. This matrix is precomputed.
4. The result is the Fourier mode coefficients in a matrix of size $[\mtotal \times \ntotal]$, which must be **copied** to the 1D array data structure introduced in the first section.

Note that the precomputed matrices are the transposed versions of those used in the previous section. The projection can also be used with the $\thet$, $\zeta$ derivatives of the Fourier basis, using again precomputed matrices.

The total number of multiplications for the two matrix-matrix products is

$$
(\mnyq\, \nnyq\, 2\ntotal) + (\mtotal\, 2\mnyq\, \ntotal)
= 2\mnyq\, \ntotal(\mtotal + \nnyq)
$$
