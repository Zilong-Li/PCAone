# Mathematics of `--winem`

This note describes the working-tree implementation reviewed on 2026-10-24,
including the EMU Newton solve, PCAngsd expectation-change test, and default
over-relaxation cap of 8. The implementation is in
[`Halko.cpp`](../src/Halko.cpp), [`Data.cpp`](../src/Data.cpp), and
[`LinAlg.cpp`](../src/LinAlg.cpp). The equations describe current behavior;
they do not assert global convergence or a unique solution.

## 1. State and scales

There are \(N\) samples and \(M\) retained SNPs. Write

\[
X=[x_1,\ldots,x_M]\in\mathbb R^{N\times M},\qquad
U\in\mathbb R^{N\times r},\quad U^\top U=I_r,
\]

where \(r=\texttt{--em-k}\), defaulting to \(K=\texttt{-k}\).
The model basis \(U\) is `Uem`; the final output basis can have a different
rank \(K\). The sketch has \(\ell=r+s\) columns, where \(s\) is the effective
EM oversampling count, and
\(\Omega\in\mathbb R^{N\times\ell}\).

Population allele frequencies \(f_j\) are estimated before winEM and held
fixed during its passes. The fitting matrix is centered but not genetically
standardized:

| Mode | Matrix entry during fitting | Reconstruction-to-frequency map |
|---|---|---|
| EMU | \(x_{ij}=g_{ij}/P-f_j\), for observed dosage \(g_{ij}\) and ploidy \(P\) | \(f_j+u_i^\top c_j\) |
| PCAngsd | \(x_{ij}=\mathbb E[g_{ij}\mid D_{ij},\pi_{ij}]-2f_j\), with diploid \(g_{ij}\in\{0,1,2\}\) | \(f_j+\tfrac12u_i^\top c_j\) |

EMU initially puts missing entries at zero, the centered population mean.
PCAngsd initially uses genotype posteriors under \(\pi_{ij}=f_j\).
Before the first window basis exists, in-core winEM uses these initial
entries. Out-of-core EMU receives missing entries as NaN and replaces them
with zero when the basis is empty.

## 2. EMU: solve each SNP conditional on the current basis

For one SNP, let \(O\) and \(A\) be its observed and missing sample sets.
Use \(U_O,U_A\) for the corresponding rows of \(U\), and define

\[
a=-f,\qquad b=1-f,\qquad
\beta=U_O^\top x_O,\qquad
C(v)=\operatorname{clip}(v,a,b).
\]

The unregularized self-consistency equations are

\[
x_A=C(U_Ac),\qquad
c=U^\top x=\beta+U_A^\top C(U_Ac).
\tag{1}
\]

These are the stationary imputation equations of the ordinary EMU update
\(x_A\leftarrow C(U_AU^\top x)\) when \(U\) is held fixed. winEM solves
them further within each SNP update rather than taking only one imputation
step. Observed calls remain unchanged.

### Regression when no bound is active

If \(U_Ac\) lies inside the box, (1) becomes

\[
(I-U_A^\top U_A)c=\beta,
\qquad I-U_A^\top U_A=U_O^\top U_O.
\tag{2}
\]

Thus the coefficients regress the observed calls on the observed rows of
the sample basis. This explains why the E-step needs no old imputations.
It also explains possible rank deficiency when few samples observe a SNP.

### The implemented regularized objective

The implementation adds \(\lambda=10^{-10}\) and solves

\[
R_\lambda(c)=(1+\lambda)c-\beta-U_A^\top C(U_Ac)=0.
\tag{3}
\]

Define an antiderivative of clipping,

\[
\psi(t)=
\begin{cases}
at-\tfrac12a^2,&t<a,\\
\tfrac12t^2,&a\le t\le b,\\
bt-\tfrac12b^2,&t>b.
\end{cases}
\]

Then (3) is the gradient equation of

\[
\Phi_\lambda(c)=
\frac{1+\lambda}{2}\|c\|_2^2-\beta^\top c
-\sum_{i\in A}\psi(u_i^\top c).
\tag{4}
\]

Away from clipping boundaries its Hessian is

\[
J_\lambda(c)=(1+\lambda)I-
\sum_{\substack{i\in A\\a<u_i^\top c<b}}u_iu_i^\top
\succeq\lambda I.
\tag{5}
\]

Here \(u_i\) is the column vector formed from row \(i\) of \(U\).
Orthonormality of \(U\) gives the inequality. Consequently (4) is strongly
convex and has a unique minimizer, even for a rank-deficient observed design.
At a clipping boundary the code treats the call as free, choosing a valid
generalized Hessian.

The initial point solves the all-free system

\[
\bigl(I-U_A^\top U_A+\lambda I\bigr)c_0=\beta.
\]

At a subsequent point \(c\), let \(S\subseteq A\) be the calls outside the
bounds. The Newton target \(\widehat c\) solves

\[
\left(I-U_A^\top U_A+\lambda I+
\sum_{i\in S}u_iu_i^\top\right)\widehat c
=\beta+\sum_{i\in S}u_iC(u_i^\top c).
\tag{6}
\]

Set \(d=\widehat c-c\) and backtrack \(\alpha=1,1/2,1/4,\ldots\), testing

\[
\Phi_\lambda(c+\alpha d)
\le\Phi_\lambda(c)+10^{-4}\alpha R_\lambda(c)^\top d.
\tag{7}
\]

The code allows 30 line-search trials and 50 outer iterations. It declares
the SNP solve converged when

\[
\|R_\lambda(c)\|_2\le10^{-8}(1+\|c\|_2),
\tag{8}
\]

then writes \(x_A=C(U_Ac)\). A solve that reaches the iteration cap is counted
as unconverged and still writes its current clipped prediction.

The ridge slightly changes (1): at an exact solution of (3), the residual
of the unregularized coefficient equation is \(-\lambda c\).
Thus the implemented solve is a regularized approximation to EMU's
unregularized conditional fixed point, not an algebraically identical solve.

## 3. PCAngsd: take up to three posterior fixed-point steps

Let \(L_i(g)\) be the genotype likelihood for one SNP and sample. For an
individual allele frequency \(\pi_i\), form unnormalized posterior weights

\[
q_{i0}=L_i(0)(1-\pi_i)^2,\quad
q_{i1}=2L_i(1)\pi_i(1-\pi_i),\quad
q_{i2}=L_i(2)\pi_i^2,\quad Z_i=\sum_{g=0}^2q_{ig}.
\]

The posterior mean and variance are

\[
\mu_i(\pi_i)=\frac{q_{i1}+2q_{i2}}{Z_i},\qquad
v_i(\pi_i)=\frac{q_{i1}+4q_{i2}}{Z_i}-\mu_i(\pi_i)^2.
\tag{9}
\]

Starting from the stored centered expectations \(x^{(0)}\), each inner step is

\[
c^{(t)}=U^\top x^{(t)},\qquad
\pi_i^{(t)}=\operatorname{clip}
\left(f+\tfrac12u_i^\top c^{(t)},\epsilon,1-\epsilon\right),
\quad\epsilon=10^{-4},
\]

\[
x_i^{(t+1)}=\mu_i(\pi_i^{(t)})-2f.
\tag{10}
\]

The current implementation takes at most three steps, stopping sooner if
\(\max_i|x_i^{(t+1)}-x_i^{(t)}|<10^{-6}\).
It retains the last expectations for the next pass. Unlike the EMU solve,
this is a stateful partial inner solve.

Holding \(U\) fixed, its stationary equation is

\[
c=U^\top\left[\mu\!\left(
\operatorname{clip}(f\mathbf1+\tfrac12Uc,\epsilon,1-\epsilon)
\right)-2f\mathbf1\right].
\tag{11}
\]

For a certain genotype, the posterior mean does not depend on the prior.
For a flat likelihood it is \(2\pi_i\). Consequently certain-or-flat
likelihoods recover EMU-like equations after accounting for the factor of
two and the different clipping bounds. Reusing the EMU routine unchanged
would use the wrong scale and bounds.

### Missing information and the local inner rate

Differentiating the posterior gives

\[
\frac{d\mu_i}{d\pi_i}=
\frac{v_i}{\pi_i(1-\pi_i)}.
\]

For calls whose frequency is not clipped, define

\[
\rho_i=\frac{v_i}{2\pi_i(1-\pi_i)}.
\tag{12}
\]

The coefficient map in (11) has local Jacobian
\(U^\top\operatorname{diag}(\rho_i)U\), with the corresponding derivative
set to zero for frequencies strictly beyond a clipping bound. Local
contraction requires its spectral radius to be below one. A flat likelihood
has \(\rho_i=1\), a certain call has \(\rho_i=0\), and an individual
likelihood can give \(\rho_i>1\).

For over-relaxation statistics the code uses
\(\operatorname{clip}(\rho_i,0,1)\) from each call's last inner step,
counted once. This bounded statistic is a heuristic; it is not the exact
Jacobian when clipping is active or a raw \(\rho_i\) exceeds one.

## 4. Window products and basis updates

For a block of SNPs \(J\), first update its entries using the newest model
basis. Then compute

\[
T_J=X_J^\top\Omega\in\mathbb R^{|J|\times\ell},\qquad
H_h\leftarrow H_h+X_JT_J,
\tag{13}
\]

where \(h\in\{1,2\}\) chooses one of two window accumulators.
At a scheduled basis update,

\[
H=H_1+H_2=QR,\qquad R=L\Sigma W^\top,
\]

\[
\Omega\leftarrow Q,\qquad
U_{\rm raw}\leftarrow QL_{[:,1:r]}.
\tag{14}
\]

Column signs of \(\Omega\) are aligned heuristically with its previous
version; this does not change the subspace.
Equation (14) takes the leading left singular vectors of \(H\).
It does not compute a full PCA decomposition of \(X\) at each update.

If all contributions use the same matrix and sketch, then
\(H=XX^\top\Omega\): this is a block power update. During the early window
phase, the accumulators can contain contributions from different sketches
and from imputations made under different model bases. Equation (13), rather
than a single global product \(XX^\top\Omega\), describes that phase exactly.

### The band schedule

Let \(B=\texttt{-w}\), a power of two at least 4, and index passes from
\(p=0\). In-core, the current band contains

\[
b_p=\min(2^{p+1},B)
\]

nominal SNP blocks. The code assigns blocks to the first or second half of
the current band and uses both full-band updates and selected half-band or
adjacent-window updates. More precisely, for zero-based block index \(j\)
and one-based position \(i\) within the current band, define

\[
A_{pj}=[p>0]\,[j+1=2^{p-1}]\,[2^p<B].
\]

An update occurs exactly when

\[
\bigl(j+1\ge b_p\ \text{or}\ A_{pj}\bigr)
\quad\text{and}\quad
\bigl(i=b_p\ \text{or}\ i=b_p/2\ \text{or}\ A_{pj}\bigr).
\tag{15}
\]

At \(i=b_p\), the code clears \(H_1\) and restarts the position counter;
at other updates it clears \(H_2\). For \(2^p\ge B\), both accumulators
are also cleared at the start of the pass. In this full-band phase, (15)
allows only the update at the end of the pass: there is one basis update
per data pass, not an additional halfway update.

Writing \(L=\log_2B\), the in-core band first reaches \(B\) on one-based
pass \(L\). The first pass that clears both accumulators before reading
the full band is \(L+1\). That pass initializes the relaxation reference;
the first possible extrapolation is on pass \(L+2\).
For \(B=64\), these are passes 6, 7, and 8 respectively.

Out-of-core EMU uses physical I/O blocks: its band size doubles from
`bandFactor` and is capped at `nblocks`; the adjacent-window boundary is
multiplied by `bandFactor`. Rounding the physical block count can change
the early schedule. The per-SNP solve is still stateless. Out-of-core
PCAngsd winEM is currently unsupported.

## 5. Per-sample over-relaxation

In full-band passes, estimate each sample's average missing information as

\[
\bar\rho_i=
\begin{cases}
M^{-1}\sum_j\mathbf1\{g_{ij}\text{ missing}\},&\text{EMU},\\
M^{-1}\sum_j\operatorname{clip}(\rho_{ij},0,1),&\text{PCAngsd}.
\end{cases}
\tag{16}
\]

EMU counts its missing calls once. PCAngsd resets and recomputes its sums
each pass. The extrapolation weights are

\[
\omega_i=\min\left\{
\frac{1}{1-\min(\bar\rho_i,1-10^{-6})},\ \omega_{\max}
\right\},\qquad
\omega_{\max}=\texttt{--winem-relax},
\tag{17}
\]

with default \(\omega_{\max}=8\); a cap of 1 disables extrapolation.

Let \(U_-\) be the previous full-band reference basis. First remove the
arbitrary rotation between bases using orthogonal Procrustes alignment:

\[
U_{\rm raw}^\top U_-=A\Lambda B^\top,\qquad
U_a=U_{\rm raw}AB^\top.
\]

Then extrapolate the rows and reorthonormalize:

\[
\widetilde U=U_-+
\operatorname{diag}(\omega_i)(U_a-U_-),\qquad
U_+=\operatorname{qf}(\widetilde U).
\tag{18}
\]

The code replaces the leading \(r\) columns of the sketch by \(U_+\),
reorthonormalizes the entire sketch, and saves \(U_+\) as the next reference.

The motivation is a scalar local error model
\(e_{t+1}\approx\bar\rho_i e_t\). Extrapolating its increment by
\(1/(1-\bar\rho_i)\) would cancel that scalar error. The actual coupled
matrix update is more complicated; (17) is not its exact inverse Jacobian.
At an aligned fixed subspace the increment in (18) is zero, so the
extrapolation preserves that fixed subspace. This does not establish global
stability or guarantee the same solution from different initializations.

If the basis-change statistic increases during eligible full-band passes,
the code halves the cap, with a lower bound of 1. This changes future steps;
it does not roll back the current step or test an objective before accepting
the extrapolation.

## 6. Stopping criteria

The implemented `mev` is the mean projection length, not the mean squared
projection length. For two orthonormal \(N\times r\) bases,

\[
\operatorname{mev}(U,V)=\frac1r\sum_{a=1}^r\|U^\top v_a\|_2,
\qquad d_U=1-\operatorname{mev}(U_p,U_{p-1}).
\tag{19}
\]

It equals one for identical subspaces. For nearby subspaces,
\(d_U\approx(2r)^{-1}\sum_a\sin^2\theta_a\), where \(\theta_a\) are
their principal angles. Thus a small value measures a small squared
subspace increment, not the remaining error to a limiting solution.
Floating-point roundoff can make the computed \(d_U\) slightly negative.

PCAngsd also computes

\[
d_E=\frac{\sum_j\|x_{j,\rm after}-x_{j,\rm before}\|_2^2}
{\sum_j\|x_{j,\rm after}\|_2^2}.
\tag{20}
\]

The snapshots bracket that SNP's inner update in the pass. In clean
full-band passes this is the squared relative Frobenius change of the
expectation matrix. If the denominator is zero, the code sets \(d_E=0\).
EMU sets \(d_E=0\) unconditionally because it solves its SNP problems
directly rather than propagating stored imputations between passes.

Stopping requires all three conditions:

\[
\text{passes}\ge\lceil\log_2 B\rceil+3,\qquad
d_U<\tau,\qquad d_E<\tau,
\quad\tau=\texttt{--tol-em}.
\tag{21}
\]

The default \(\tau\) for winEM is \(10^{-7}\); `--maxiter` caps data passes.
For example, (20) below \(10^{-7}\) means a relative expectation increment
below approximately \(3.16\times10^{-4}\), not a relative increment below
\(10^{-7}\). Under a scalar linear contraction with rate \(q<1\), remaining
error can be roughly \(q/(1-q)\) times the latest increment; neither (19)
nor (20) estimates that rate.

## 7. Final decomposition and PCAngsd covariance

After fitting, freeze \(U_\star=U_{\rm em}\). In-core, EMU solves every SNP
once more under \(U_\star\); PCAngsd takes another up-to-three inner steps
under \(U_\star\). PCAngsd's finalization is not an unlimited inner solve.

With default genetic scaling, EMU forms

\[
Y_{ij}=\frac{\sqrt P\,X_{ij}}{\sqrt{f_j(1-f_j)}}.
\tag{22}
\]

PCAngsd first decomposes its unstandardized centered expectations,
\(Y=X\). Nondefault scaling options can change or bypass the default
transformation; (22) describes `--scale -9`. The implementation avoids
division by a near-zero standard deviation.

The final solver is warm-started and runs full bands. For a fixed final
matrix, one such pass has

\[
T=Y^\top\Omega=Q_T R_T,\qquad H=YT,
\qquad B_T=R_T^{-\top}H^\top=Q_T^\top Y^\top.
\]

If \(B_T=L_T\Sigma_TW_T^\top\), the reported rank-\(K\) decomposition is

\[
U_{\rm out}=W_{T,[:,1:K]},\qquad
S_{\rm out}=\Sigma_{T,1:K},\qquad
V_{\rm out}=Q_TL_{T,[:,1:K]}.
\tag{23}
\]

The inverse notation assumes full sketch rank; the code uses a pivoted QR
solve when constructing \(B_T\). Further power passes refine this randomized
decomposition to `--tol-rsvd`, subject to `--maxp`. Out-of-core EMU recreates
the same frozen-basis imputations and scaling as each block is read in each
final pass. Changing \(K\) changes the requested final decomposition rank;
the fitted model rank remains \(r\).
The `.sigvals` are the entries of \(S_{\rm out}\), and the `.eigvals` are
\(S_{{\rm out},a}^2/M\), the eigenvalues of \(YY^\top/M\) approximated by
this decomposition.

For BEAGLE input and default scaling, covariance construction is a separate
step after (23). Starting from the final unstandardized expectations
\(X_{\rm fin}\), form

\[
\Pi=\operatorname{clip}\left(
\mathbf1 f^\top+\tfrac12U_\star U_\star^\top X_{\rm fin},
\epsilon,1-\epsilon\right).
\]

Recompute the posterior once at \(\Pi\), giving means \(\mu_{ij}\) and
variances \(v_{ij}\), and define

\[
Z_{ij}=\frac{\mu_{ij}-2f_j}{\sqrt{2f_j(1-f_j)}}.
\]

For polymorphic sites, the written covariance satisfies

\[
C_{ih}=\frac1M\sum_j Z_{ij}Z_{hj}\quad(i\ne h),
\]

\[
C_{ii}=\frac1M\sum_j
\frac{v_{ij}+(\mu_{ij}-2f_j)^2}{2f_j(1-f_j)}.
\tag{24}
\]

Thus the diagonal includes posterior genotype uncertainty. The `.eigvecs2`
of a BEAGLE run contains the top \(K\) eigenvectors of this corrected
covariance; these need not equal the sample singular vectors in `.eigvecs`
from the unstandardized expectation matrix.

## 8. Relation to ordinary EM-PCA

Ordinary EM-PCA alternates a rank-\(r\) PCA decomposition and an imputation
or posterior update. winEM interleaves those operations across SNP windows:

\[
\boxed{\text{update a SNP block from the newest }U
\;\longrightarrow\;\text{accumulate }X_JX_J^\top\Omega
\;\longrightarrow\;\text{update }U,\Omega.}
\]

In the clean full-band phase, a stationary solution must satisfy the
conditional SNP equations together with a stationary sketch/basis update.
If the sketch has converged to the relevant leading invariant subspace,
the basis condition is the usual leading-PC condition. The EMU ridge,
finite inner and outer tolerances, sketch approximation, and possible
multiple stationary solutions qualify any claim that the methods produce
exactly the same answer.

The convex objective (4) belongs to an individual EMU SNP solve with a fixed
basis. It is not a demonstrated global objective that every window update
or extrapolation decreases. Similarly, solving PCAngsd's equation (11)
does not establish maximum likelihood for a separately specified
low-rank genotype-likelihood model.

The equations assume a fitted basis exists. The current implementation still
accepts `--maxiter 0` for winEM, leaving that basis empty; the reviewed BEAGLE
path can crash during finalization. That is an implementation edge case,
not a defined zero-pass version of the mathematics above.
