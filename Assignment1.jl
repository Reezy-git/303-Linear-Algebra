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
Householder Reflectors:
"""

# ╔═╡ 9cdf0092-5d3e-4521-8d4a-07c962a0772d
begin
	M = [
		3 4
		0 -2
		0 1
		0 2
	];
	H = Matrix(I, 4, 4)
	m_=[
		0
		-2
		1
		2
	]
	θ=norm(m_)
	e₂=[
		0
		1
		0
		0
	]
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

We must find `θ=norm(M)` this is valid since we have the special case where M is a column vector and subtract `θ*e₂` from our `M`

$θ=norm(M)=3$

$\vec{v}=M-θ*e₂=\left[\begin{array}{c}
0\\
-5\\
1\\
2
\end{array}\right]$

We then find 

$β=\frac{2}{\vec{v}^{T}\vec{v}}$

Using this to find H

$H = I - β \vec{v}\vec{v}^{T} = \frac{1}{3}\left[\begin{array}{cc}
3 & 0 & 0 & 0\\
0 & -2 & 1 & 2\\
0 & 1 & 2.8 & -0.4\\
0 & 2 & -0.4 & 2.2
\end{array}\right]$

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
"""

# ╔═╡ 90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
md"""
Givens Rotations:
"""

# ╔═╡ 50b5be53-e0a1-46f9-978d-6512f993eae8
function givens_helper(a, b)
        if b == 0
            c = sign(a)
            if c == 0
                c = 1
            end
        elseif a == 0
            c = 0
            s = sign(b)
            r = abs(b)
        elseif abs(a) > abs(b)
            t = b/a
            u = sign(a) * √(1 + t^2)
            c = 1 / u
            s = c * t
            r = a * u
        else
            t = a / b
            u = sign(b) * √(1 + t^2)
            s = 1 / u
            c = s * t
            r = b * u
        end
        return c, s, r
end;

# ╔═╡ 63ef1671-08de-4556-aede-c87165243993
md"""
re-do givens helper from class notes
"""

# ╔═╡ f1e8717a-b80d-452c-a808-851d5735aa7a
givens_helper(4, 2)

# ╔═╡ 6937d7e4-cb7e-4d3c-8247-6b832e267038
md"""
We first must find c, s & r using the algorithm defined by the code for the function above: givens_helper with a = 

$c = -2,\quad s = 1,\quad r=√5$

We use this to create our first givens matrix

$G_{32}= \left[\begin{array}{cc}
1 & 0 & 0 & 0\\
0 & \frac{-2}{\sqrt{5}} & \frac{1}{\sqrt{5}} & 0\\
0 & \frac{-1}{\sqrt{5}} & \frac{-2}{\sqrt{5}} & 0\\
0&0&0&1
\end{array}\right] \quad\quad G_{32}A = \left[\begin{array}{c}
3 & 4\\
0 & \frac{5}{\sqrt{5}}\\
0 & 0\\
0 & 2
\end{array}\right]$
"""

# ╔═╡ 825e4aca-45d4-42af-8c50-1ab5ee670bad
md"""
### *Question Two*
"""

# ╔═╡ 41fc5b67-869c-4b1c-bb5e-cbcd6d2217c1
begin
	Q = 1/3 * [
		2 -2 -1
		1 2 -2
		2 1 2
	]

	R =[
		1 2
		0 1
		0 0
	]
	b =[
		1
		1
	]
	Aᵀ=Q*R
end;

# ╔═╡ 470ae356-7ebe-4d9a-b582-b3796e36d7c1
md"""
Given the QR factors of Aᵀ we can find the economy QR

$A^T=Y\hat{R}=\frac{1}{3}\left[\begin{array}{cc}
2 & -2\\
1 & 2\\
2 & 1
\end{array}\right]\left[\begin{array}{cc}
1 & 2\\
0 & 1
\end{array}\right]$

This gives us

$\hat{R}^T\vec{u}=\left[\begin{array}{cc}
1\\1
\end{array}\right]$
We can solve this system using forward substitution

$\left[\begin{array}{cc}
1 & 0 & : & 1\\
2 & 1 & : & 1
\end{array}\right]\quad\Rightarrow\quad x_{1}=1,\quad x_{2}=-1$

We use this $\vec{x}$ to 
"""

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
# ╠═9cdf0092-5d3e-4521-8d4a-07c962a0772d
# ╟─64a5e043-7b71-4559-86cc-a316bfc3d668
# ╠═90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
# ╠═50b5be53-e0a1-46f9-978d-6512f993eae8
# ╟─63ef1671-08de-4556-aede-c87165243993
# ╠═f1e8717a-b80d-452c-a808-851d5735aa7a
# ╠═6937d7e4-cb7e-4d3c-8247-6b832e267038
# ╟─825e4aca-45d4-42af-8c50-1ab5ee670bad
# ╟─41fc5b67-869c-4b1c-bb5e-cbcd6d2217c1
# ╠═470ae356-7ebe-4d9a-b582-b3796e36d7c1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
