### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 22220222-fdce-41e9-96b0-50b8109487b2
using LinearAlgebra # importing the LinearAlgebra package

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
a) Householder Reflectors:
"""

# ╔═╡ 9cdf0092-5d3e-4521-8d4a-07c962a0772d
begin
	M = [
		3 4
		0 -2
		0 1
		0 2
	]
	H = Matrix(I, 4, 4)
	m_=[0; -2; 1; 2]  # being lazy and not slicing from M then padding
	θ=norm(m_)
	e₂=[0; 1; 0; 0]  # in my actual code implementation I slice from identity matrix
	v_ = m_ - θ*e₂
	β = 2 / norm(v_)^2
	H = H - β*v_*transpose(v_)
end;

# ╔═╡ 64a5e043-7b71-4559-86cc-a316bfc3d668
md"""
Since column 1 is correctly formatted we only need to perform one Householder reflection to diagonalize this matrix. We take the portion of the matrix we will be acting apon and add padding. This results in:

$M=\left[
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
We must find `θ=norm(M)` this is valid since we have the special case where M is a column vector and subtract $θ*e₂$ from our $M$

$θ=norm(M)=3$
"""

# ╔═╡ 609409e6-84a2-4c3d-9d6e-4ebc64a08e5a
md"""
$\vec{v}=M-θ*e₂=\left[\begin{array}{c}
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
Using this to find H

$H = I - β \vec{v}\vec{v}^{T} = \frac{1}{3}\left[\begin{array}{cc}
3 & 0 & 0 & 0\\
0 & -2 & 1 & 2\\
0 & 1 & 2.8 & -0.4\\
0 & 2 & -0.4 & 2.2
\end{array}\right]$
"""

# ╔═╡ 087e7c14-7f4f-4c26-809d-bf0da209ae24
md"""
Thus we find 

$Q = \frac{1}{3}\left[\begin{array}{cc}
3 & 0 & 0 & 0\\
0 & -2 & 1 & 2\\
0 & 1 & 2.8 & -0.4\\
0 & 2 & -0.4 & 2.2
\end{array}\right]
\quad R =\left[\begin{array}{cc}
3 & 4\\
0 & 3\\
0 & 0\\
0 & 0
\end{array}\right]$

-----
"""

# ╔═╡ 90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
md"""
b) Givens Rotations:
"""

# ╔═╡ 50b5be53-e0a1-46f9-978d-6512f993eae8
function givens_helper(A, row , col)
	"""Takes input, matrix, and row/column target and returns values for Givens rotations"""
	c, s = A[col,col], A[row,col]
	r² = c^2 + s^2
    return c, s, r²
end

# ╔═╡ 6937d7e4-cb7e-4d3c-8247-6b832e267038
md"""
We first must find c, s & r using the algorithm defined by the code for the function above: givens_helper with A, 3, 2

$c = -2,\quad s = 1,\quad r=√5$
"""

# ╔═╡ ed0c0cd2-2f34-4c65-810b-e711d51dd50f
md"""
We put create our Givens rotation matrix like so:

$G_{32}= \left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & \frac{c}{r} & \frac{s}{r} & 0\\
0 & -\frac{s}{r} & \frac{c}{r} & 0\\
0&0&0&1
\end{array}\right]$
"""

# ╔═╡ 64c214dd-6424-40d2-a0c5-6fa459384fde
md"""
We use this to create our first givens matrix

$G_{32}= \left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & \frac{-2}{\sqrt{5}} & \frac{1}{\sqrt{5}} & 0\\
0 & \frac{-1}{\sqrt{5}} & \frac{-2}{\sqrt{5}} & 0\\
0&0&0&1
\end{array}\right] \quad\quad G_{32}A = \left[\begin{array}{c}
3 & 4\\
0 & \sqrt{5}\\
0 & 0\\
0 & 2
\end{array}\right]$
"""

# ╔═╡ c4d119cb-5c97-4cbb-a0fe-671b3235e439
md"""
We repeat the process with row 4 column 2.

$G_{42}=\left[\begin{array}{cc}
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
\end{array}\right]$
"""

# ╔═╡ a15aeec1-01ff-4b06-966c-148f4da3d220
md"""
We obtain Q in the following way.

$Q=G_{32}^TG_{42}^T=\left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & -\frac{2}{3} & -\frac{1}{\sqrt{5}} & \frac{4}{3\sqrt{5}}\\
0 & \frac{1}{3} & -\frac{2}{\sqrt{5}} & -\frac{2}{3\sqrt{5}}\\
0 & \frac{2}{3} & 0 & \frac{\sqrt{5}}{3}
\end{array}\right]$
"""

# ╔═╡ aafcf37b-8b71-4bb2-95d6-b6860db77e3c
begin
	G32=[
		1 0 0 0
		0 -2/√5 1/√5 0
		0 -1/√5 -2/√5 0
		0 0 0 1
	]
	G32*M
	c, s, r² = givens_helper(G32*M, 4, 2)
	r=sqrt(r²)
	G42 = [
		1 0 0 0
		0 c/r 0 s/r
		0 0 1 0
		0 -s/r 0 c/r
	]
	G42*G32*M
	G32'*G42'

end;

# ╔═╡ 825e4aca-45d4-42af-8c50-1ab5ee670bad
md"""
### *Question Two*
Given the QR factors of Aᵀ we can find the economy QR

$A^T=Y\hat{R}=\frac{1}{3}\left[\begin{array}{cc}
2 & -2\\
1 & 2\\
2 & 1
\end{array}\right]\left[\begin{array}{cc}
1 & 2\\
0 & 1
\end{array}\right]$
"""

# ╔═╡ 798513af-4235-4c1d-8470-d0a4744af373
md"""
This gives us

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
\end{array}\right]\quad\Rightarrow\quad x_{1}=1,\quad x_{2}=-1$
"""

# ╔═╡ db055454-a28f-4d72-9db5-26ee5b8f7656
md"""
We use this $\vec{x}$ to find $\vec{b}$ where

$Y\vec{x}=\vec{b} \quad\Rightarrow\quad \frac{1}{3}\left[\begin{array}{cc}
2 & -2\\
1 & 2\\
2 & 1
\end{array}\right]
\left[\begin{array}{cc}
1\\
-1
\end{array}\right]=\vec{b}\quad\Rightarrow\quad\vec{b}=\frac{1}{3}\left[\begin{array}{cc}
4\\
-1\\
1
\end{array}\right]$

-----
"""

# ╔═╡ 303b1fe8-0fb0-4807-80c0-3cea5f2532e3
md"""
b) We can find the general solution to the equation $A\vec{x}=b$ by doing Gaussian elimination on the system

$\frac{1}{3}\left[\begin{array}{cc}
2 & 1 & 2\\
2 & 4 & 5
\end{array}\right]
\left[\begin{array}{cc}x_1\\x_2\\x_3\end{array}\right] =
\left[\begin{array}{cc} 1\\1 \end{array}\right]
\quad\Rightarrow\quad
\left[\begin{array}{cc} 
2 & 1 & 2 & : & 3\\
2 & 4 & 5 & : & 3
\end{array}\right]$

$\Rightarrow\quad
\left[\begin{array}{cc} 
2 & 1 & 2 & : & 3\\
0 & 3 & 3 & : & 0
\end{array}\right]
\quad\Rightarrow\quad
x_1=3/2\quad x_2=-t \quad x_3=-t$
"""

# ╔═╡ 83bc2406-76b8-4410-a9ea-41d413920dd3
md"""
A particular solution to this equation could be

$\left[\begin{array}{cc}x_1\\x_2\\x_3\end{array}\right] =
\left[\begin{array}{cc}\frac{3}{2}\\0\\0\end{array}\right]$
"""

# ╔═╡ d725f5f6-d6a0-455d-bcbb-178553d7f22e
md"""
### *Question 3*
"""

# ╔═╡ c28a11e2-0563-4c6b-a05c-119baa9a8e9b
#place holder

# ╔═╡ 24591fea-a98f-40d0-bfd5-369c2988e0ae
md"""
### *Question 4*
"""

# ╔═╡ 5aac61a0-1037-474b-99c2-54341d4fea51
md"""
a) Going with the eigen-value route 

$det(A-λ I) \quad\Rightarrow\quad det\left[\begin{array}{cc}
1-λ & 1 & 0\\
1 & a-λ & -2\\
0 & -2 & 10-λ
\end{array}\right]$
"""

# ╔═╡ 5dc95a33-5e81-4a41-b35e-738daafe0f39
md"""
Gives us the characteristic polynomial

$(1-λ)((a-λ)(10-λ)-4)-(10-λ)=0$

$-λ^3+11λ^2+aλ^2-11λa-5λ+10a-14=0$
"""

# ╔═╡ 7a7864f6-150d-4549-bebe-d7b98b7f75da
md"""
We can think of this as a function of $a$ and set $λ$ to $0$, then solve for $a$

$10a-14=0 \quad\Rightarrow\quad a=\frac{7}{5}$
"""

# ╔═╡ 9cae94f4-917b-473e-a1e8-fd93ba58a82c
md"""
That is to say there is only one solution for a that gives us an eigen value of 0. i.e. this equation can only cross 0 at one point. This tells us that only one of the roots to the charactaristic polynomial can change sign, and it happens when $a=\frac{7}{5}$. 
"""

# ╔═╡ ea73172f-e64f-47ad-a484-c1f6b92e925f
md"""
We can solve again when $a∈\lbrace 1,\, 2 \rbrace$

$a=1 \;\Rightarrow\; \lambda\approx-0.21,\,1.76,\,10.43$
$a=2 \;\Rightarrow\; \lambda\approx0.252,\,2.27,\,10.48$

So we can see when $a>\frac{7}{5}$ 

we get positive eigen values and thus a positive definite matrix. $\square$

-----
"""

# ╔═╡ 285fd45b-3c0d-41e4-a3e2-e79530a279ef
md"""
b) A is positive semi-definite when $a=\frac{7}{5}$ as above.

-----
"""

# ╔═╡ d9f1b7e8-b6f9-4489-a3b7-8b7021b3854b
md"""
c) When $a=5$ we get the matrix

$\left[\begin{array}{cc}
1 & 1 & 0\\
1 & 5 & -2\\
0 & -2 & 10
\end{array}\right]$
"""

# ╔═╡ 861a46b3-9752-43ac-a2cb-72e6051be1d4
md"""
We make this matrix upper diagonal through row reduction

$\xrightarrow[]{R₂=R₂-R₁}\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 4 & -2\\
0 & -2 & 10
\end{array}\right]
\quad\xrightarrow[]{R₃=R₃-\frac{1}{2}R₂}\quad
\left[\begin{array}{cc}
1 & 1 & 0\\
0 & 4 & -2\\
0 & 0 & 9
\end{array}\right]$
"""

# ╔═╡ 1e61ec72-cab3-41e7-9f0d-f88b9975add8
begin
	A=[
		1 1 0
		1 5 -2
		0 -2 10
	]
	U=[
		1 1 0
		0 4 -2
		0 0 9
	]
	L=A*inv(U)
end

# ╔═╡ 01885659-043a-4494-b9d7-2d5e043da4a0
md"""
We find 

$L=AU^{-1}=\left[\begin{array}{cc}
1 & 0 & 0\\
1 & 1 & 0\\
0 & -\frac{1}{2} & 1
\end{array}\right]$
"""

# ╔═╡ bdd7d69d-b975-4645-9241-850606732402
md"""
Giving us 

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

--------
"""

# ╔═╡ 4dcc313b-8580-44d4-8349-5edbc34e2d99
md"""
d) We need R such that

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
We observe that if $R\vec{x}=\vec{y}$ then $R^{T}\vec{y}=\vec{b}$. Solving 

$R^T\vec{x}=\vec{b}\;\Leftrightarrow\;\left[\begin{array}{cc}
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


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised
julia_version = "1.7.3"
manifest_format = "2.0"
[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╟─d5305a50-2338-11ed-1517-1defdfa17b48
# ╠═22220222-fdce-41e9-96b0-50b8109487b2
# ╟─505f4ad4-6f69-4101-8745-ac1e75f367a4
# ╟─910ed377-7895-4da6-ad1d-d2fa5a78d8ca
# ╟─a01b5da6-9d2d-429b-ac55-93b8a3858a99
# ╟─9cdf0092-5d3e-4521-8d4a-07c962a0772d
# ╟─64a5e043-7b71-4559-86cc-a316bfc3d668
# ╟─36bf99dc-9374-4c8e-9c7f-91657411f06f
# ╟─609409e6-84a2-4c3d-9d6e-4ebc64a08e5a
# ╟─1140f421-8205-4f58-bd91-061257737717
# ╟─35fe45e1-476e-4dca-a496-b27abac02aef
# ╟─087e7c14-7f4f-4c26-809d-bf0da209ae24
# ╟─90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
# ╠═50b5be53-e0a1-46f9-978d-6512f993eae8
# ╟─6937d7e4-cb7e-4d3c-8247-6b832e267038
# ╟─ed0c0cd2-2f34-4c65-810b-e711d51dd50f
# ╟─64c214dd-6424-40d2-a0c5-6fa459384fde
# ╟─c4d119cb-5c97-4cbb-a0fe-671b3235e439
# ╟─a15aeec1-01ff-4b06-966c-148f4da3d220
# ╟─aafcf37b-8b71-4bb2-95d6-b6860db77e3c
# ╟─825e4aca-45d4-42af-8c50-1ab5ee670bad
# ╟─798513af-4235-4c1d-8470-d0a4744af373
# ╟─74a3bb34-448e-457e-85cb-b2689152361e
# ╟─db055454-a28f-4d72-9db5-26ee5b8f7656
# ╟─303b1fe8-0fb0-4807-80c0-3cea5f2532e3
# ╟─83bc2406-76b8-4410-a9ea-41d413920dd3
# ╟─d725f5f6-d6a0-455d-bcbb-178553d7f22e
# ╠═c28a11e2-0563-4c6b-a05c-119baa9a8e9b
# ╟─24591fea-a98f-40d0-bfd5-369c2988e0ae
# ╟─5aac61a0-1037-474b-99c2-54341d4fea51
# ╟─5dc95a33-5e81-4a41-b35e-738daafe0f39
# ╟─7a7864f6-150d-4549-bebe-d7b98b7f75da
# ╟─9cae94f4-917b-473e-a1e8-fd93ba58a82c
# ╟─ea73172f-e64f-47ad-a484-c1f6b92e925f
# ╟─285fd45b-3c0d-41e4-a3e2-e79530a279ef
# ╟─d9f1b7e8-b6f9-4489-a3b7-8b7021b3854b
# ╟─861a46b3-9752-43ac-a2cb-72e6051be1d4
# ╠═1e61ec72-cab3-41e7-9f0d-f88b9975add8
# ╟─01885659-043a-4494-b9d7-2d5e043da4a0
# ╟─bdd7d69d-b975-4645-9241-850606732402
# ╟─4dcc313b-8580-44d4-8349-5edbc34e2d99
# ╟─869cc088-da24-4af4-82b1-7d05d2061cc2
# ╟─ed80a134-d420-4bfd-a81f-4a12039d6424
# ╟─d9bf4eb8-74c2-4e51-8270-833c5acad20f
# ╟─0e864931-9bb1-4c49-a43d-5c8f31eb0670
# ╠═49a127ce-12c2-4a4b-8f1b-9647bfe4ac8b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
