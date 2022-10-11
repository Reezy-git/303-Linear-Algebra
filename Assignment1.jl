### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 22220222-fdce-41e9-96b0-50b8109487b2
using LinearAlgebra

# ╔═╡ a783ca5c-e85f-483e-b3bb-88e25aa56a83
begin
	using Plots
	
	function a(λ)
	    return a = (λ^3+5λ-11λ^2+14)/(λ^2-11λ+10)
	end

	plot(
	    -4:0.1:14, a, 
	    xlims=(-1.5,13), ylims=(-1,11), 
	    xlabel="λ", ylabel="a",
	    title="Solutions to Charactaristic Polynomial",
	    label="Solns",
		gridalpha=0.3, minorgrid=true, thickness_scaling=1.3, showaxis=true
	)
	scatter!([0], [1.4], label="(0, 1.4)")
	vline!([0], label="", linecolor=RGB(0,0,0), lw=1.6)
	hline!([0], label="", linecolor=RGB(0,0,0), lw=1.6)
end

# ╔═╡ d5305a50-2338-11ed-1517-1defdfa17b48
md"""
# **MATH303-22S2** - _Assignment 1_

`Due Date`: **Monday 19th September 5pm**

Submission by: Richard Hodges, Henry Hastings

Student IDs: 11139318, 35556669

Code is written in Julia and presented in a Pluto.jl notebook.
"""

# ╔═╡ 505f4ad4-6f69-4101-8745-ac1e75f367a4
md"""
### *Question One*
"""

# ╔═╡ 910ed377-7895-4da6-ad1d-d2fa5a78d8ca
md"""
$A=\left[\begin{array}{cc}
3 & 4\\
0 & -2\\
0 & 1\\
0 & 2
\end{array}\right]$
"""

# ╔═╡ a01b5da6-9d2d-429b-ac55-93b8a3858a99
md"""
###### a) Householder Reflectors:
"""

# ╔═╡ 64a5e043-7b71-4559-86cc-a316bfc3d668
md"""
Since column 1 is correctly formatted we only need to perform one Householder reflection to diagonalize this matrix. We take the portion of the matrix we will be acting upon and add padding. This results in:

$\vec{m}=\left[
	\begin{array}{cc}
	0\\
	-2\\
	1\\
	2
	\end{array}
\right]$
"""

# ╔═╡ 36bf99dc-9374-4c8e-9c7f-91657411f06f
md"""
We must find $θ=\lVert\vec{m}\rVert₂$ and subtract $θ\times e₂$ from our $\vec{m}$:

$θ=\lVert\vec{m}\rVert₂=3$
"""

# ╔═╡ 609409e6-84a2-4c3d-9d6e-4ebc64a08e5a
md"""
$\vec{v}=\vec{m}-θ*e₂=\left[\begin{array}{c}
0\\
-5\\
1\\
2
\end{array}\right]$
"""

# ╔═╡ 1140f421-8205-4f58-bd91-061257737717
md"""
We then find 

$β=\frac{2}{\vec{v}^{T}\vec{v}}$
"""

# ╔═╡ 35fe45e1-476e-4dca-a496-b27abac02aef
md"""
We use the following formula to find $H$:

$H = I - β \vec{v}\vec{v}^{T} = \frac{1}{3}\left[\begin{array}{cc}
3 & 0 & 0 & 0\\
0 & -2 & 1 & 2\\
0 & 1 & 2.8 & -0.4\\
0 & 2 & -0.4 & 2.2
\end{array}\right]$
"""

# ╔═╡ 087e7c14-7f4f-4c26-809d-bf0da209ae24
md"""
Thus, noting there is only one householder matrix used in this transformation, we find:

$Q= H = \frac{1}{3}\left[\begin{array}{cc}
3 & 0 & 0 & 0\\
0 & -2 & 1 & 2\\
0 & 1 & 2.8 & -0.4\\
0 & 2 & -0.4 & 2.2
\end{array}\right]
\quad R =QA=\left[\begin{array}{cc}
3 & 4\\
0 & 3\\
0 & 0\\
0 & 0
\end{array}\right]$

-----
"""

# ╔═╡ 90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
md"""
###### b) Givens Rotations:
"""

# ╔═╡ 6937d7e4-cb7e-4d3c-8247-6b832e267038
md"""
We first must find $c$ and $s$ for row 3 column 2 of A using the given rule:

$c =\frac{a_{22}}{\sqrt{a_{22}^2+a_{32}^2}} = \frac{-2}{\sqrt{5}},\quad s = \frac{a_{32}}{\sqrt{a_{22}^2+a_{32}^2}} = \frac{1}{\sqrt{5}}$
"""

# ╔═╡ ed0c0cd2-2f34-4c65-810b-e711d51dd50f
md"""
We create our Givens rotation matrix like so:

$G_{32}= \left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & c & s & 0\\
0 & -s & c & 0\\
0&0&0&1
\end{array}\right]$
"""

# ╔═╡ 64c214dd-6424-40d2-a0c5-6fa459384fde
md"""
This results in our first Givens matrix:

$G_{32}= \left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & \frac{-2}{\sqrt{5}} & \frac{1}{\sqrt{5}} & 0\\
0 & \frac{-1}{\sqrt{5}} & \frac{-2}{\sqrt{5}} & 0\\
0&0&0&1
\end{array}\right]
\quad\Rightarrow\quad
G_{32}A = \left[\begin{array}{c}
3 & 4\\
0 & \sqrt{5}\\
0 & 0\\
0 & 2
\end{array}\right]$
"""

# ╔═╡ c4d119cb-5c97-4cbb-a0fe-671b3235e439
md"""
We repeat the process with row 4 column 2, using $\,G_{32}A\,$ as our new  "$A$ ":

$c =\frac{a_{22}}{\sqrt{a_{22}^2+a_{42}^2}} = \frac{\sqrt{5}}{3},\quad s = \frac{a_{42}}{\sqrt{a_{22}^2+a_{42}^2}} = \frac{2}{3}$


$\Rightarrow\quad G_{42}=\left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & \frac{\sqrt{5}}{3} & 0 & \frac{2}{3}\\
0 & 0 & 1 & 0\\
0 & -\frac{2}{3} & 0 & \frac{\sqrt{5}}{3}
\end{array}\right]
\quad\Rightarrow\quad
G_{42}G_{32}A=
\left[\begin{array}{cc}
3 & 4\\
0 & 3\\
0 & 0\\
0 & 0
\end{array}\right]=R$
"""

# ╔═╡ a15aeec1-01ff-4b06-966c-148f4da3d220
md"""
We obtain Q in the following way:

$Q=G_{32}^TG_{42}^T=\left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & -\frac{2}{3} & \frac{-1}{\sqrt{5}} & \frac{4}{3\sqrt{5}}\\
0 & \frac{1}{3} & \frac{-2}{\sqrt{5}} & \frac{-2}{3\sqrt{5}}\\
0 & \frac{2}{3} & 0 & \frac{\sqrt{5}}{3}
\end{array}\right]$
"""

# ╔═╡ 4bb8c855-3b9e-47cb-b581-c303e27ae925
md"""
### *Question Two*
As $A\vec{x}=\vec{b}$ is an underdetermined system, following the notes in section 4.1.3 we must find the economy QR factorisation of $A^T$, then solve $\hat{R}^T\vec{u}=\vec{b}$ for $\vec{u}$ and then finally set $\vec{x}=Y\vec{u}$ to find the LSS for $A\vec{x}=\vec{b}$. 
"""

# ╔═╡ 825e4aca-45d4-42af-8c50-1ab5ee670bad
md"""

Given the QR factors of Aᵀ we can find its economy QR factorisation by trimming off the row of zeros in $R$ and then removing the necessary column in $Q$:

$A^T=Y\hat{R}=\frac{1}{3}\left[\begin{array}{cc}
2 & -2\\
1 & 2\\
2 & 1
\end{array}\right]\left[\begin{array}{cc}
1 & 2\\
0 & 1
\end{array}\right]=\frac{1}{3}\left[\begin{array}{cc}
2 & 2\\
1 & 4\\
2 & 5
\end{array}\right]$
"""

# ╔═╡ 798513af-4235-4c1d-8470-d0a4744af373
md"""
This gives us the following expression to evaluate, following our steps

$\hat{R}^T\vec{u}=\left[\begin{array}{cc}
1\\1
\end{array}\right]$
"""

# ╔═╡ 74a3bb34-448e-457e-85cb-b2689152361e
md"""
We can solve this system using forward substitution

$\left[\begin{array}{cc}
1 & 0 & : & 1\\
2 & 1 & : & 1
\end{array}\right]\quad\Rightarrow\quad u_{1}=1,\quad u_{2}=-1 \quad\Rightarrow\quad\vec{u}=\left[\begin{array}{cc}
1\\-1
\end{array}\right]$
"""

# ╔═╡ db055454-a28f-4d72-9db5-26ee5b8f7656
md"""
We use this $\,\vec{u}\,$ to find $\,\vec{x}\,$ continuing our steps:

$Y\vec{u}=\vec{x} \quad\Rightarrow\quad \frac{1}{3}\left[\begin{array}{cc}
2 & -2\\
1 & 2\\
2 & 1
\end{array}\right]
\left[\begin{array}{cc}
1\\
-1
\end{array}\right]=\vec{x}\quad\Rightarrow\quad\vec{x}=\frac{1}{3}\left[\begin{array}{cc}
4\\
-1\\
1
\end{array}\right]$

-----
"""

# ╔═╡ 303b1fe8-0fb0-4807-80c0-3cea5f2532e3
md"""
b) We can find the general solution to the equation $A\vec{x}=\vec{b}$ by doing Gaussian elimination on the system, noting we already know $A$ from our expression of $A^T$

$\frac{1}{3}\left[\begin{array}{cc}
2 & 1 & 2\\
2 & 4 & 5
\end{array}\right]
\left[\begin{array}{cc}x_1\\x_2\\x_3\end{array}\right] =
\left[\begin{array}{cc} 1\\1 \end{array}\right]
\quad\Leftrightarrow\quad
\left[\begin{array}{cc} 
2 & 1 & 2 & : & 3\\
2 & 4 & 5 & : & 3
\end{array}\right]\quad\Rightarrow\quad
\left[\begin{array}{cc} 
2 & 1 & 2 & : & 3\\
0 & 3 & 3 & : & 0
\end{array}\right]$

$\quad\Rightarrow\quad
x_1=\frac{3-t}{2},\quad x_2=-t, \quad x_3=t$
"""

# ╔═╡ 83bc2406-76b8-4410-a9ea-41d413920dd3
md"""
A particular solution to this equation could be

$\left[\begin{array}{cc}x_1\\x_2\\x_3\end{array}\right] =
\left[\begin{array}{cc}\frac{3}{2}\\0\\0\end{array}\right]$
"""

# ╔═╡ 5d220621-c193-4e90-8aff-2610abbc70d6
md"""

### *Question 3*

We are given:

$B=A+\vec{u}\vec{v}^T$
$B^{-1}=A^{-1} - αA^{-1}\vec{u}\vec{v}^TA^{-1}$
"""

# ╔═╡ ae920fb5-c334-4ce1-abc5-8ded3a365dad
md"""
Thus we must show:

$I=(A+\vec{u}\vec{v}^T)(A^{-1} - αA^{-1}\vec{u}\vec{v}^TA^{-1})$
"""

# ╔═╡ 0d90647e-0e0b-409d-b7c7-2714874586a3
md"""
Giving the following:
 
$I=I+\vec{u}\vec{v}^TA^{-1} - α\vec{u}\vec{v}^TA^{-1} - α\vec{u}\vec{v}^TA^{-1}\vec{u}\vec{v}^TA^{-1}$

$α\vec{u}\vec{v}^TA^{-1} + α\vec{u}\vec{v}^TA^{-1}\vec{u}\vec{v}^TA^{-1}=\vec{u}\vec{v}^TA^{-1}$

$α\vec{u}\vec{v}^T+α\vec{u}\vec{v}^TA^{-1}\vec{u}\vec{v}^T = \vec{u}\vec{v}^T$

$α\vec{u}\vec{v}^T+α\vec{u}(\vec{v}^TA^{-1}\vec{u})\vec{v}^T = \vec{u}\vec{v}^T$

$α\vec{u}\vec{v}^T+α\vec{u}\vec{v}^T(\vec{v}^TA^{-1}\vec{u}) = \vec{u}\vec{v}^T$

$\vec{u}\vec{v}^T(α+α\vec{v}^TA^{-1}\vec{u}) = \vec{u}\vec{v}^T$
"""

# ╔═╡ 9018d133-bdc8-4f2f-9292-81a91f951bbe
md"""
We note that $\,\vec{v}^TA^{-1}\vec{u}\,$, and thus $\,α+α\vec{v}^TA^{-1}\vec{u}\,$ is a scalar meaning

$α+α\vec{v}^TA^{-1}\vec{u}=1$

"""

# ╔═╡ 39f0015a-3bd2-4743-9732-ccbdd5a81dc6
md"""
Rearranging for $\,α\,$:

$α=\frac{1}{1+\vec{v}^TA^{-1}\vec{u}}$

Therefore we have a value for $\,α\,$ which satisfies $\,BB^{-1}=I\,$ making the statement for $\,B^{-1}\,$ true.

Therefore:

$B^{-1}=A^{-1} - αA^{-1}\vec{u}\vec{v}^TA^{-1}$

With:

$α=\frac{1}{1+\vec{v}^TA^{-1}\vec{u}}$
"""

# ╔═╡ 5aac61a0-1037-474b-99c2-54341d4fea51
md"""

### *Question 4*


a) Going with the eigenvalue route 

$det(A-λ I) \quad=\quad det\left[\begin{array}{cc}
1-λ & 1 & 0\\
1 & a-λ & -2\\
0 & -2 & 10-λ
\end{array}\right]$
"""

# ╔═╡ 5dc95a33-5e81-4a41-b35e-738daafe0f39
md"""
Gives us the characteristic polynomial

$P(λ)=(1-λ)((a-λ)(10-λ)-4)-(10-λ)$

$P(λ)=-λ^3+11λ^2+aλ^2-11λa-5λ+10a-14$
"""

# ╔═╡ 5b2add48-45dd-434d-962d-1cd12aa26700
md"""
We set $P(λ) = 0$ to find λ the eigenvalues of $A$

$-λ^3+11λ^2+aλ^2-11λa-5λ+10a-14=0$
"""


# ╔═╡ 7a7864f6-150d-4549-bebe-d7b98b7f75da
md"""
We can set $λ$ to $0$, then solve for $a$

$10a-14=0 \quad\Rightarrow\quad a=\frac{7}{5}$
"""

# ╔═╡ 9cae94f4-917b-473e-a1e8-fd93ba58a82c
md"""
This means there is only one solution for $\,a\,$ ($\frac{7}{5}$) that gives us an eigen value of $0$. The cubic $P(λ)$ for any $a\,\in\,\mathbb{R}\,$ is smooth and continuous. This tells us that the roots to the characteristic polynomial can change sign, but it can only happen when $a=\frac{7}{5}$. 
"""

# ╔═╡ ea73172f-e64f-47ad-a484-c1f6b92e925f
md"""
We can solve again when $a∈\lbrace 1,\, 2 \rbrace$

$a=1 \;\Rightarrow\; \lambda\approx-0.21,\,1.76,\,10.43$
$a=2 \;\Rightarrow\; \lambda\approx0.252,\,2.27,\,10.48$

So we can see when $a>\frac{7}{5}$ we get positive eigenvalues and thus a positive definite matrix. $\,\square$

"""

# ╔═╡ 1a18de8f-f637-4832-9594-fe8089f57f4c
md"""
Supplementarily we can make $a$ the subject of the system and plot the solution of for $a$ with respect to $\,λ$. A little algebra yields:

$a=\frac{λ^3-11λ^2+5λ+14}{(λ-1)(λ-10)}$
"""

# ╔═╡ b4ad26cf-8330-43ce-92c8-946a07fa8392
md"""
This has vertical asymptotes at $λ=1$ and $λ=10$. We can plot this as well. For every value of $a$ we can draw a horizontal line and get the 3 roots of our characteristic polynomial, the eigenvalues of A. By observation it is clear the only time all 3 eigenvalues of A are positive is when $a>\frac{7}{5}$."""

# ╔═╡ dfb880f8-0795-45f9-a497-0508e007c573
md"""
-----"""

# ╔═╡ 285fd45b-3c0d-41e4-a3e2-e79530a279ef
md"""
b) A is positive semi-definite when $\,a\ge\frac{7}{5}\,$, as this is when $λ≥0$ for all three eigenvalues of $A$.

-----
"""

# ╔═╡ d9f1b7e8-b6f9-4489-a3b7-8b7021b3854b
md"""
c) When $a=5$ we get the matrix

$A=\left[\begin{array}{cc}
1 & 1 & 0\\
1 & 5 & -2\\
0 & -2 & 10
\end{array}\right]$
"""

# ╔═╡ 861a46b3-9752-43ac-a2cb-72e6051be1d4
md"""
We make this matrix upper diagonal through row reduction

$A\xrightarrow[]{R₂=R₂-R₁}\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 4 & -2\\
0 & -2 & 10
\end{array}\right]
\quad\xrightarrow[]{R₃=R₃+\frac{1}{2}R₂}\quad
\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 4 & -2\\
0 & 0 & 9
\end{array}\right]=U$
"""

# ╔═╡ 75329d4e-1cc7-41f2-8f57-3c6b135412c6
md"""
We can find the inverse of $U$ by using row reduction to make $U$ $I$ and doing the same operations on $I$

$\left[\begin{array}{cc}
1 & 1 & 0 & : & 1 & 0 & 0\\
0 & 4 & -2 & : & 0 & 1 & 0\\
0 & 0 & 9 & : & 0 & 0 & 1
\end{array}\right] \quad\Rightarrow\quad \left[\begin{array}{cc}
1 & 1 & 0 & : & 1 & 0 & 0\\
0 & 4 & -2 & : & 0 & 1 & 0\\
0 & 0 & 1 & : & 0 & 0 & \frac{1}{9}
\end{array}\right]$
$\Rightarrow\quad \left[\begin{array}{cc}
1 & 1 & 0 & : & 1 & 0 & 0\\
0 & 1 & 0 & : & 0 & \frac{1}{4} & \frac{1}{18}\\
0 & 0 & 1 & : & 0 & 0 & \frac{1}{9}
\end{array}\right]\quad\Rightarrow\quad \left[\begin{array}{cc}
1 & 0 & 0 & : & 1 & \frac{-1}{4} & \frac{-1}{18}\\
0 & 1 & 0 & : & 0 & \frac{1}{4} & \frac{1}{18}\\
0 & 0 & 1 & : & 0 & 0 & \frac{1}{9}
\end{array}\right]$
"""

# ╔═╡ 06b10a24-c820-4488-8373-2fe557e323b8
md"""
Thus:

$U^{-1}=\left[\begin{array}{cc}
1 & \frac{-1}{4} & \frac{-1}{18}\\
0 & \frac{1}{4} & \frac{1}{18}\\
0 & 0 & \frac{1}{9}
\end{array}\right]$

"""

# ╔═╡ 01885659-043a-4494-b9d7-2d5e043da4a0
md"""
Multiplying $U^{-1}$ by $A$ provides:

$L=AU^{-1}=\left[\begin{array}{cc}
1 & 0 & 0\\
1 & 1 & 0\\
0 & -\frac{1}{2} & 1
\end{array}\right]$
"""

# ╔═╡ bdd7d69d-b975-4645-9241-850606732402
md"""
This gives us:

$A=LDL^T=
\left[\begin{array}{cc}
1 & 0 & 0\\
1 & 1 & 0\\
0 & -\frac{1}{2} & 1
\end{array}\right]
\left[\begin{array}{cc}
1 & 0 & 0\\
0 & 4 & 0\\
0 & 0 & 9
\end{array}\right]
\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 1 & -\frac{1}{2}\\
0 & 0 & 1
\end{array}\right]$

Where $\,D\,$ is a diagonal matrix with the diagonal values of $\,U$; the modified Cholesky factorisation of $\,A$

--------
"""

# ╔═╡ 4dcc313b-8580-44d4-8349-5edbc34e2d99
md"""
d) We find a matrix $\,R\,$ such that:

$R=D^{\frac{1}{2}}L^T=\left[\begin{array}{cc}
\sqrt{1} & 0 & 0\\
0 & \sqrt{4} & 0\\
0 & 0 & \sqrt{9}
\end{array}\right]
\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 1 & -\frac{1}{2}\\
0 & 0 & 1
\end{array}\right] = 
\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 2 & -1\\
0 & 0 & 3
\end{array}\right]$
"""

# ╔═╡ 869cc088-da24-4af4-82b1-7d05d2061cc2
md"""
We observe that if $\,A\vec{x}=\vec{b}\,$ then $\,R^TR\vec{x}=\vec{b}\,$. Hence, if we let $\,R\vec{x}=\vec{y}\,$ we can solve $\,R^T\vec{y}=\vec{b}\,$ and then solve $\,R\vec{x}=\vec{y}\,$ afterwards, which is significantly easier due to $\,R\,$ already being in upper-triangular form.

$R^T\vec{y}=\vec{b}\;\Leftrightarrow\;\left[\begin{array}{cc}
1 & 0 & 0 & : & 1\\
1 & 2 & 0 & : & 1\\
0 & -1 & 3 & : & 1
\end{array}\right]$
"""

# ╔═╡ ed80a134-d420-4bfd-a81f-4a12039d6424
md"""
This gives us 

$\vec{y}=\left[\begin{array}{cc} 1 \\ 0\\ \frac{1}{3}\end{array}\right]\quad\Rightarrow\quad R\vec{x}=\left[\begin{array}{cc} 1 \\ 0\\ \frac{1}{3}\end{array}\right]$
"""

# ╔═╡ d9bf4eb8-74c2-4e51-8270-833c5acad20f
md"""
$R\vec{x}=\vec{y}\;\Leftrightarrow\;
\left[\begin{array}{cc}
1 & 1 & 0 & : & 1\\
0 & 2 & -1 & : & 0\\
0 & 0 & 3 & : & \frac{1}{3}
\end{array}\right]
\quad\Rightarrow\quad
\vec{x} = \left[\begin{array}{cc}
\frac{17}{18}\\
\frac{1}{18}\\
\frac{1}{9}
\end{array}\right]$
"""

# ╔═╡ 0e864931-9bb1-4c49-a43d-5c8f31eb0670
md"""
### *Question 5*
"""

# ╔═╡ 49a127ce-12c2-4a4b-8f1b-9647bfe4ac8b
results = [
	1948 10.3; 1952 10.4; 1956 10.5; 1960 10.2
	1964 10.0; 1968 9.95; 1972 10.14; 1976 10.06
	1980 10.25; 1984 9.99; 1988 9.92; 1992 9.96
	1996 9.84; 2000 9.87; 2004 9.85; 2008 9.69
	2012 9.63; 2016 9.81; 2020 9.80
]

# ╔═╡ c1c9383a-ca3c-4d6d-ab99-c9c494de67a8
trace1 = scatter(results[:,1], results[:,2], title="Olympics sprint times", label="Winning time", xlabel="Year", ylabel="Time (sec)")


# ╔═╡ 277b39bd-640f-4e35-8655-1a90663d1d62
md"""
a) 

 $A_1$ is the matrix that represents the linear approximation of $\,T$, that is $\,T_1= c_0 +c_1Y\,$ where $\,Y\,$ is the vector of years, being the first column of the results matrix, and the $\,c\,$'s are scalars (note: $\,Y\,$ is the leftmost column of $\,A_1\,$). We aim to find the LSS of $\,A_1 \vec{x} = \vec{b}\,$, where $\,\vec{b}\,$ is the vector of times.
"""

# ╔═╡ de1cffa9-90a8-4da1-81de-8eae3a759467
A₁=hcat(results[:,1], ones(19,1))

# ╔═╡ e9bae10b-581c-46aa-985b-201b9f187ab0
b=results[:,2]

# ╔═╡ 22e83600-a041-4b4f-b851-ad21ff77cfe9
md"""
We then use the qr( ) command to form the QR factorisation of $\,A_1$
"""

# ╔═╡ aefef568-8858-4374-88c1-227b6a84b7ac
Q₁, R₁ = qr(A₁)

# ╔═╡ d3fb2cf5-d2ac-430a-a195-ad449445442c
md"""
$A\vec{x}=\vec{b}\quad\Rightarrow\quad Q_1R_1\vec{x}=\vec{b}$

Thus we let $\,R_1\vec{x}=\vec{y_1}\,$, and multiply both sides on the left by $\,Q^T\,$ to get $\,Q_1^T\vec{b}=\vec{y_1}$

"""

# ╔═╡ 17f69850-60fb-46e1-a594-15c0684fd4a1
md"""
We find $Q_1^T\vec{b}=\vec{y_1}$, and then solve $R_1\vec{x}=\vec{y_1}\,$ by back substitution.
"""

# ╔═╡ 7af09d5f-a2b6-4b88-a1b8-b4167b9773fd
y₁ = Q₁'b

# ╔═╡ e960c4aa-8b3e-4b02-85a4-c4bd8f6bb44c
c0₁=y₁[2]/R₁[2,2]

# ╔═╡ 389b3833-37c8-4c78-acbd-e6a486e1506a
c1₁=(y₁[1]-R₁[1,2]*c0₁)/R₁[1,1]

# ╔═╡ 3aad13b3-167b-4418-bd80-e696521867ec
x₁=[c1₁; c0₁]

# ╔═╡ 9121e6fa-bf3e-49d5-8aab-eba223592643
md"""
Using the backslash command $\,\vec{x}=A_1$\ $\vec{b}\,$ yields the same result:
"""

# ╔═╡ 77c702eb-0df3-46d8-a342-3d0669919266
A₁\b

# ╔═╡ dc6c9ae5-2c76-4c70-90d6-e0d4995c8c59
md"""
In both cases, this yields a function 

$\,t= 28.169 - 0.00915351y\,$ 

to predict a time, $\,t\,$, in a specific year, $\,y$.
"""

# ╔═╡ 7d666b9f-65ff-45f7-a20b-71e1abd101d1
md"""
 $T_1$ represents the outputs of this function, that is the predicted results from the degree one fitting polynomial.

"""

# ╔═╡ d46c9512-bb51-4fb7-b09b-0b5a99475b37
T₁=A₁*x₁

# ╔═╡ 7bef8142-53a2-4a25-8be1-560a09a622b5
begin
	p₁ = plot(trace1)
	plot!(p₁, results[:,1], T₁, label="Linear Approximation")
end

# ╔═╡ 6e48740c-b073-4724-a28c-cec97515548e
md"""
The error of the liner approximation was calculated by the 2-norm of $\,\vec{b}-T_{1}$

"""

# ╔═╡ 7aa23b7f-bf62-414e-b1a6-2745fc91c910

error₁=norm((b-T₁),2)

# ╔═╡ 55900b5e-ec55-471b-a07c-6d396ca4b5f2
md"""
-----
"""

# ╔═╡ cf7fcde0-bc12-4650-8b29-68d79a499583
md"""

b)

 $A_2$ is the matrix that represents the quadratic approximation of $\,T$, that is $\,T_2= c_0 +c_1Y +c_2Y^2\,$ where $\,Y\,$ is the vector of years, and the $\,c\,$'s are scalars (note: $\,Y^2\,$ is the leftmost column of $\,A_2\,$). We aim to find the LSS of $\,A_2 \vec{x} = \vec{b}\,$, where $\,\vec{b}\,$ is the vector of times.
"""

# ╔═╡ 1f574ae9-0233-4636-9769-ee32bfbb09c2
A₂= hcat(results[:,1].^2, A₁)

# ╔═╡ c892c55c-5049-451f-bd45-45f96c88ebae
md"""
We then use the qr( ) command to form the QR factorisation of $\,A_2$
"""

# ╔═╡ a8f93713-7f0a-407d-98de-e041b2fc56fe
Q₂, R₂ = qr(A₂)

# ╔═╡ 0d19e731-6d97-41db-916b-958eadfba5a2
md"""
Using the same logic as in 5a) we find $\,Q_2^T\vec{b}=\vec{y_2}\,$, and then solve $\,R_2\vec{x}=\vec{y_2}$
"""

# ╔═╡ 43a07f5c-2ff8-431d-b924-d2dff0120232
y₂=Q₂'b

# ╔═╡ b88853cc-0e8d-4991-b925-1c521b3be7a5
c0₂ = y₂[3] / R₂[3,3]

# ╔═╡ 1b53d348-944c-441c-bf0b-a3fe1a413e9a
c1₂ = (y₂[2] - R₂[2,3]*c0₂) / R₂[2,2]

# ╔═╡ 3761ec47-5e29-4c06-9723-82aeba84e8ab
c2₂ = (y₂[1] - R₂[1,3]*c0₂ - R₂[1,2]*c1₂) / R₂[1,1]

# ╔═╡ 1567a906-d428-4b10-9cc4-ce459e5e5568
x₂=[c2₂, c1₂, c0₂]

# ╔═╡ 62b5ccd5-1b52-48b0-a464-99dd1da420f6
md"""
Using the backslash command $\,\vec{x}=A_2$\ $\vec{b}\,$ yields the same result:
"""

# ╔═╡ d7971b50-5dee-42de-a614-39b0dfbf1510
A₂\b

# ╔═╡ 1e0c9bc7-2ad4-4d2e-84c1-b60ea2992dd2
md"""
In both cases, this yields a function 

$\,t= 256.458 -0.239311y + 5.80035*10^{-5}y^2\,$ 

to predict a time, $\,t$, in a specific year, $\,y$.
"""

# ╔═╡ dce9b7c8-9f62-4755-a177-0414b7020d30
md"""
 $T_2\,$ represents the outputs of this function, that is the predicted results from the degree two fitting polynomial.

"""

# ╔═╡ 5e424167-fe2e-4098-9b5b-d5df661a0162
T₂= A₂*x₂

# ╔═╡ 441c552a-6e63-4093-979d-cbd703907324
begin
	p₂ = plot(trace1)
	plot!(p₂, results[:,1], T₂, label="Quadratic Approximation")
end

# ╔═╡ e82140e0-0fb2-4d23-a741-35dd24254115
md"""
The error of the quadratic approximation was calculated by the 2-norm of $\,\vec{b}-T_{2}$

"""

# ╔═╡ 95afd6a0-0b81-4055-8e0e-07070a757596
error₂=norm((b-T₂),2)

# ╔═╡ f1960673-6eaa-4c89-ab22-47f2ea54a354
md"""
-----
"""

# ╔═╡ b9c671c8-365e-445a-b817-debccac08e9f
md"""

c)

 $A_3\,$ is the matrix that represents the linear approximation of $\,T\,$, that is $\,T_3= c_0 +c_1Y +c_2Y^2 +c_3Y^3\,$ where $\,Y\,$ is the vector of years, and the $\,c\,$'s are scalars (note: $\,Y^3\,$ is the leftmost column of $\,A_3\,$). We aim to find the LSS of $\,A_3 \vec{x} = \vec{b}\,$, where $\,\vec{b}\,$ is the vector of times.
"""

# ╔═╡ 19016495-ad30-4678-b644-c3966debd03f
A₃=hcat(results[:,1].^3, A₂)

# ╔═╡ 6d50fdef-edaf-4c4b-99c3-0f62f1edb5fb
md"""
We then use the qr( ) command to form the QR factorisation of $A_3$
"""

# ╔═╡ f8ad08bc-1cd8-4749-8c18-d8e8cdb4745f
Q₃, R₃ = qr(A₃)

# ╔═╡ 2d588d91-702d-4f1b-97f2-17ffbed7c342
md"""
Using the same logic as in 5a) and 5b) we find $\,Q_3^T\vec{b}=\vec{y_3}\,$, and then solve $\,R_3\vec{x}=\vec{y_3}$
"""

# ╔═╡ fb02188a-1649-4bee-b900-ec3125a776a3
y₃=Q₃'b

# ╔═╡ 69329cb1-a1b6-40b8-a5c1-f17703413573
c0₃ = y₃[4] / R₃[4,4]

# ╔═╡ bcadbafd-ba0c-411d-9488-c941ead6edb4
c1₃= (y₃[3] - R₃[3,4]*c0₃) / R₃[3,3]

# ╔═╡ 1908087b-7c35-46d2-9867-e9c120a1a7b5
c2₃ = (y₃[2] - R₃[2,4]*c0₃ - R₃[2,3]*c1₃) / R₃[2,2]

# ╔═╡ 359ab131-4458-4c81-968e-6e9c7dc1ca3a
c3₃ = (y₃[1] - R₃[1,4]*c0₃ - R₃[1,3]*c1₃ - R₃[1,2]*c2₃) / R₃[1,1]

# ╔═╡ 303a6960-ae39-427e-937f-ad9cbe1f9699
x₃ = [c3₃; c2₃; c1₃; c0₃]

# ╔═╡ 535952e1-7a21-4536-a8e0-d305266ced80
md"""
The QR factorisation yields a function 

$t= -5528.27 +8.50904y -0.00435177y^2 + 7.40889*10^{-7}y^3$ 

to predict a time, $\,t\,$, in a specific year, $\,y$.
"""

# ╔═╡ 449d3005-c163-44a2-966d-7ec7807685ef
md"""
 $T_3\,$ represents the outputs of this function, that is the predicted results from the degree three fitting polynomial with the QR factorization.

"""

# ╔═╡ 4a591b8d-606e-4464-bd1e-12b97683aafd
T₃ = A₃*x₃

# ╔═╡ 917c7a8e-ffa4-4f45-842b-e22177d0f017
begin
	p₃ = plot(trace1)
	plot!(p₃, results[:,1], T₃, label="QR Cubic Approximation")
end

# ╔═╡ 8c47ba4d-8379-4de2-a45a-2d55031950dc
md"""
The error of this cubic approximation was calculated by the 2-norm of $\,\vec{b}-T_{3}$

"""

# ╔═╡ c5cc2c9f-cf9e-4c3f-af5a-fa2affa4439e
error₃ = norm((b-T₃),2)

# ╔═╡ 00a16c77-f3ee-47e9-97f4-3edea6ffe3a3
md"""
In this case, using the backslash command yielded a different LSS. This is because we have violated the assumption of the QR factorized matrix having full rank (part 4.1 in the lecture notes), with $A_3$ having rank 3 and thus a nullity of 1. As a result, the QR factorization does not occur correctly, providing an incorrect LSS. 

"""

# ╔═╡ dd08f116-7c1f-4955-adc5-067e671ad785
md"""
The rank of $A_3$ is 3, as seen by the Julia rank command.

"""

# ╔═╡ 16737963-a3c6-416a-a3c0-d229f38aa296
r = rank(A₃)

# ╔═╡ d7bf2fe8-d735-4e54-b18c-96df3e009018
x3 = A₃\b

# ╔═╡ 3181b6ab-a987-488a-8ec2-b7a61e0e0964
md"""
The backslash command yields a function 

$t= 0.000225011 +0.14879y -0.000137756y^2 + 3.291131\times10^{-8}y^3$

to predict a time, $\,t\,$, in a specific year, $\,y$
"""

# ╔═╡ b6ad9dd5-86e0-4150-acf0-6fe2a5df57a1
md"""
 $T3\,$, not to be confused with $\,T_3\,$ represents the outputs of this function, that is the predicted results from the degree three fitting polynomial using the backslash command.

"""

# ╔═╡ a5d807d2-32cf-440b-9fb7-da3fd3d1f422
T3 = A₃*x3

# ╔═╡ ee946d95-1173-46a1-9818-c45d2db586a9
begin
	p3 = plot(trace1)
	plot!(p3, results[:,1], T3, label="Backslash Cubic Approximation")
end

# ╔═╡ 3a00c878-e61e-4422-b9d6-61db67a939d5
md"""
The error of this cubic approximation was calculated by the 2-norm of $\,\vec{b}-T3\,$

"""

# ╔═╡ 8120b7bc-f1c6-490a-8f62-716c2c3c1d5d
error3 = norm((T3-b),2)

# ╔═╡ d50d9181-6d95-4c14-b891-05ae2f29d043
md"""
d)
"""

# ╔═╡ 4f08e254-a354-4931-bf6a-d021d5cbcb3a


begin
	p_all = plot(trace1)
	plot!(p_all, results[:,1], T3, label="Backslash Cubic Approximation")
	plot!(p_all, results[:,1], T₃, label="QR Cubic Approximation")
	plot!(p_all, results[:,1], T₁, label="Linear Approximation")
	plot!(p_all, results[:,1], T₂, label="Quadratic Approximation")
end

# ╔═╡ 8677aa1c-2026-4e38-a8e8-b87707f9c34d
md"""Note that the $\,x^3\,$ coefficient of both our cubic models is virtually zero. Because of this, the cubic approximations are virtually the same as the quadratic approximation. From this we can conclude there is very little benefit to actually adding a cubic approximation. 

-------
"""

# ╔═╡ 415d4096-4fae-4cb3-9a60-c710c3c519a3
md"""e) 2024 predictions are simply found by setting $\,x\,$ to 2024 in the specific functions, which is the same as assessing the dot product of the vector 

$\,\vec{Y_i}=[2024^i \quad 2024^{i-1}\quad ...\quad 2024^1 \quad 1]^T\,$
"""

# ╔═╡ 12ec190d-5cdb-4771-b300-18c3d7f3e484
md"""
and the specific solution vector $\,\vec{x_i}\,$ for each degree $\,i\,$ approximation, up to degree three. The difference between the quadratic and cubic approximations, especially between the quadratic and backslash cubic approximations, demonstrates how the approximations give the same result to a hundredth, or even less, of a second, further highlighting the similarities we see in their graphic representations.
"""

# ╔═╡ e6d0ce9b-3a4d-484e-8523-e55b32c4f32b
md""" Linear approximation in seconds:"""

# ╔═╡ 9589aee5-ffb4-4095-8b77-8bdf0bb20f9b
[2024 1] * x₁  # Linear

# ╔═╡ e7318975-9cd5-4185-a950-3bcb05d71ddd
md""" Quadratic approximation in seconds:"""

# ╔═╡ 7ced8a88-4167-4450-880f-81fc9e349526
[2024^2 2024 1]*x₂  # Quadratic

# ╔═╡ 238d62cd-fa00-42b0-805c-645826b2736d
md""" Cubic approximation in seconds, using QR:"""

# ╔═╡ d2d85219-56b6-4829-89d7-710b0ef09896
[2024^3 2024^2 2024 1]*x₃  # Cubic QR

# ╔═╡ 1d4be66e-3ce3-41a7-82ea-ffb55bcf44df
md""" Cubic approximation in seconds, using the backslash command:"""

# ╔═╡ fff1c6ce-1230-469e-891a-4870b969c9e4
[2024^3 2024^2 2024 1]*x3  # Cubic backslash

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Plots = "~1.33.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "dc4405cee4b2fe9e1108caec2d760b7ea758eca2"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.5"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "5856d3031cdb1f3b2b6340dfdc66b6d9a149a374"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.2.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cf0a9940f250dc3cb6cc6c6821b4bf8a4286cf9c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "3697c23d09d5ec6f2088faa68f0d926b6889b5be"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.67.0+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "59ba44e0aa49b87a8c7a8920ec76f8afe87ed502"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.3.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "6872f9594ff273da6d13c7c1a1545d5a8c7d0c1c"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "6062b3b25ad3c58e817df0747fc51518b9110e5f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.33.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "e7eac76a958f8664f2718508435d058168c7953d"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─d5305a50-2338-11ed-1517-1defdfa17b48
# ╠═22220222-fdce-41e9-96b0-50b8109487b2
# ╟─505f4ad4-6f69-4101-8745-ac1e75f367a4
# ╟─910ed377-7895-4da6-ad1d-d2fa5a78d8ca
# ╟─a01b5da6-9d2d-429b-ac55-93b8a3858a99
# ╟─64a5e043-7b71-4559-86cc-a316bfc3d668
# ╟─36bf99dc-9374-4c8e-9c7f-91657411f06f
# ╟─609409e6-84a2-4c3d-9d6e-4ebc64a08e5a
# ╟─1140f421-8205-4f58-bd91-061257737717
# ╟─35fe45e1-476e-4dca-a496-b27abac02aef
# ╟─087e7c14-7f4f-4c26-809d-bf0da209ae24
# ╟─90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
# ╟─6937d7e4-cb7e-4d3c-8247-6b832e267038
# ╟─ed0c0cd2-2f34-4c65-810b-e711d51dd50f
# ╟─64c214dd-6424-40d2-a0c5-6fa459384fde
# ╟─c4d119cb-5c97-4cbb-a0fe-671b3235e439
# ╟─a15aeec1-01ff-4b06-966c-148f4da3d220
# ╟─4bb8c855-3b9e-47cb-b581-c303e27ae925
# ╟─825e4aca-45d4-42af-8c50-1ab5ee670bad
# ╟─798513af-4235-4c1d-8470-d0a4744af373
# ╟─74a3bb34-448e-457e-85cb-b2689152361e
# ╟─db055454-a28f-4d72-9db5-26ee5b8f7656
# ╟─303b1fe8-0fb0-4807-80c0-3cea5f2532e3
# ╟─83bc2406-76b8-4410-a9ea-41d413920dd3
# ╟─5d220621-c193-4e90-8aff-2610abbc70d6
# ╟─ae920fb5-c334-4ce1-abc5-8ded3a365dad
# ╟─0d90647e-0e0b-409d-b7c7-2714874586a3
# ╟─9018d133-bdc8-4f2f-9292-81a91f951bbe
# ╟─39f0015a-3bd2-4743-9732-ccbdd5a81dc6
# ╟─5aac61a0-1037-474b-99c2-54341d4fea51
# ╟─5dc95a33-5e81-4a41-b35e-738daafe0f39
# ╟─5b2add48-45dd-434d-962d-1cd12aa26700
# ╟─7a7864f6-150d-4549-bebe-d7b98b7f75da
# ╟─9cae94f4-917b-473e-a1e8-fd93ba58a82c
# ╟─ea73172f-e64f-47ad-a484-c1f6b92e925f
# ╟─1a18de8f-f637-4832-9594-fe8089f57f4c
# ╟─b4ad26cf-8330-43ce-92c8-946a07fa8392
# ╟─a783ca5c-e85f-483e-b3bb-88e25aa56a83
# ╟─dfb880f8-0795-45f9-a497-0508e007c573
# ╟─285fd45b-3c0d-41e4-a3e2-e79530a279ef
# ╟─d9f1b7e8-b6f9-4489-a3b7-8b7021b3854b
# ╟─861a46b3-9752-43ac-a2cb-72e6051be1d4
# ╟─75329d4e-1cc7-41f2-8f57-3c6b135412c6
# ╟─06b10a24-c820-4488-8373-2fe557e323b8
# ╟─01885659-043a-4494-b9d7-2d5e043da4a0
# ╟─bdd7d69d-b975-4645-9241-850606732402
# ╟─4dcc313b-8580-44d4-8349-5edbc34e2d99
# ╟─869cc088-da24-4af4-82b1-7d05d2061cc2
# ╟─ed80a134-d420-4bfd-a81f-4a12039d6424
# ╟─d9bf4eb8-74c2-4e51-8270-833c5acad20f
# ╟─0e864931-9bb1-4c49-a43d-5c8f31eb0670
# ╟─49a127ce-12c2-4a4b-8f1b-9647bfe4ac8b
# ╟─c1c9383a-ca3c-4d6d-ab99-c9c494de67a8
# ╟─277b39bd-640f-4e35-8655-1a90663d1d62
# ╟─de1cffa9-90a8-4da1-81de-8eae3a759467
# ╟─e9bae10b-581c-46aa-985b-201b9f187ab0
# ╟─22e83600-a041-4b4f-b851-ad21ff77cfe9
# ╠═aefef568-8858-4374-88c1-227b6a84b7ac
# ╟─d3fb2cf5-d2ac-430a-a195-ad449445442c
# ╟─17f69850-60fb-46e1-a594-15c0684fd4a1
# ╠═7af09d5f-a2b6-4b88-a1b8-b4167b9773fd
# ╠═e960c4aa-8b3e-4b02-85a4-c4bd8f6bb44c
# ╠═389b3833-37c8-4c78-acbd-e6a486e1506a
# ╠═3aad13b3-167b-4418-bd80-e696521867ec
# ╟─9121e6fa-bf3e-49d5-8aab-eba223592643
# ╠═77c702eb-0df3-46d8-a342-3d0669919266
# ╟─dc6c9ae5-2c76-4c70-90d6-e0d4995c8c59
# ╟─7d666b9f-65ff-45f7-a20b-71e1abd101d1
# ╟─d46c9512-bb51-4fb7-b09b-0b5a99475b37
# ╟─7bef8142-53a2-4a25-8be1-560a09a622b5
# ╟─6e48740c-b073-4724-a28c-cec97515548e
# ╠═7aa23b7f-bf62-414e-b1a6-2745fc91c910
# ╟─55900b5e-ec55-471b-a07c-6d396ca4b5f2
# ╟─cf7fcde0-bc12-4650-8b29-68d79a499583
# ╟─1f574ae9-0233-4636-9769-ee32bfbb09c2
# ╟─c892c55c-5049-451f-bd45-45f96c88ebae
# ╠═a8f93713-7f0a-407d-98de-e041b2fc56fe
# ╟─0d19e731-6d97-41db-916b-958eadfba5a2
# ╠═43a07f5c-2ff8-431d-b924-d2dff0120232
# ╠═b88853cc-0e8d-4991-b925-1c521b3be7a5
# ╠═1b53d348-944c-441c-bf0b-a3fe1a413e9a
# ╠═3761ec47-5e29-4c06-9723-82aeba84e8ab
# ╠═1567a906-d428-4b10-9cc4-ce459e5e5568
# ╟─62b5ccd5-1b52-48b0-a464-99dd1da420f6
# ╠═d7971b50-5dee-42de-a614-39b0dfbf1510
# ╟─1e0c9bc7-2ad4-4d2e-84c1-b60ea2992dd2
# ╟─dce9b7c8-9f62-4755-a177-0414b7020d30
# ╠═5e424167-fe2e-4098-9b5b-d5df661a0162
# ╟─441c552a-6e63-4093-979d-cbd703907324
# ╟─e82140e0-0fb2-4d23-a741-35dd24254115
# ╠═95afd6a0-0b81-4055-8e0e-07070a757596
# ╟─f1960673-6eaa-4c89-ab22-47f2ea54a354
# ╟─b9c671c8-365e-445a-b817-debccac08e9f
# ╠═19016495-ad30-4678-b644-c3966debd03f
# ╟─6d50fdef-edaf-4c4b-99c3-0f62f1edb5fb
# ╠═f8ad08bc-1cd8-4749-8c18-d8e8cdb4745f
# ╟─2d588d91-702d-4f1b-97f2-17ffbed7c342
# ╟─fb02188a-1649-4bee-b900-ec3125a776a3
# ╠═69329cb1-a1b6-40b8-a5c1-f17703413573
# ╠═bcadbafd-ba0c-411d-9488-c941ead6edb4
# ╠═1908087b-7c35-46d2-9867-e9c120a1a7b5
# ╠═359ab131-4458-4c81-968e-6e9c7dc1ca3a
# ╠═303a6960-ae39-427e-937f-ad9cbe1f9699
# ╟─535952e1-7a21-4536-a8e0-d305266ced80
# ╟─449d3005-c163-44a2-966d-7ec7807685ef
# ╟─4a591b8d-606e-4464-bd1e-12b97683aafd
# ╟─917c7a8e-ffa4-4f45-842b-e22177d0f017
# ╟─8c47ba4d-8379-4de2-a45a-2d55031950dc
# ╠═c5cc2c9f-cf9e-4c3f-af5a-fa2affa4439e
# ╟─00a16c77-f3ee-47e9-97f4-3edea6ffe3a3
# ╟─dd08f116-7c1f-4955-adc5-067e671ad785
# ╠═16737963-a3c6-416a-a3c0-d229f38aa296
# ╠═d7bf2fe8-d735-4e54-b18c-96df3e009018
# ╟─3181b6ab-a987-488a-8ec2-b7a61e0e0964
# ╟─b6ad9dd5-86e0-4150-acf0-6fe2a5df57a1
# ╠═a5d807d2-32cf-440b-9fb7-da3fd3d1f422
# ╟─ee946d95-1173-46a1-9818-c45d2db586a9
# ╟─3a00c878-e61e-4422-b9d6-61db67a939d5
# ╠═8120b7bc-f1c6-490a-8f62-716c2c3c1d5d
# ╟─d50d9181-6d95-4c14-b891-05ae2f29d043
# ╟─4f08e254-a354-4931-bf6a-d021d5cbcb3a
# ╟─8677aa1c-2026-4e38-a8e8-b87707f9c34d
# ╟─415d4096-4fae-4cb3-9a60-c710c3c519a3
# ╟─12ec190d-5cdb-4771-b300-18c3d7f3e484
# ╟─e6d0ce9b-3a4d-484e-8523-e55b32c4f32b
# ╟─9589aee5-ffb4-4095-8b77-8bdf0bb20f9b
# ╟─e7318975-9cd5-4185-a950-3bcb05d71ddd
# ╟─7ced8a88-4167-4450-880f-81fc9e349526
# ╟─238d62cd-fa00-42b0-805c-645826b2736d
# ╟─d2d85219-56b6-4829-89d7-710b0ef09896
# ╟─1d4be66e-3ce3-41a7-82ea-ffb55bcf44df
# ╟─fff1c6ce-1230-469e-891a-4870b969c9e4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
