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

Submission by: Richard Hodges, Henry

Student IDs: 11139318, 

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

Since column 1 is correctly formatted we only need to perform one Householder reflection to diagonalize this matrix. We take the portion of the matrix we will be acting apon and add padding. This results in:
"""

# ╔═╡ 3930a39b-c26b-4a42-a3da-8550e2fce0f9
md"""
$M=\left[
	\begin{array}{cc}
	0\\
	-2\\
	1\\
	2
	\end{array}
\right]$
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
We must find `θ=norm(M)` this is valid since we have the special case where M is a column vector and subtract `θ*e₂` from our `M`
"""

# ╔═╡ 882a0b1c-85d8-475f-be5e-1828faa10196
md"""
$θ=norm(M)=3$
"""

# ╔═╡ 7a949596-5c6e-4faa-8971-154774836859
md"""
$\vec{v}=M-θ*e₂=\left[\begin{array}{c}
0\\
-5\\
1\\
2
\end{array}\right]$
"""

# ╔═╡ ace59a76-703d-49eb-bea3-479c245ac022
H

# ╔═╡ fe69bf18-119a-48ea-bee5-1ba5faa0b678
md"""
The outterproduct of our vector *v* gives us the householder matrix H

$H=\frac{1}{3}\left[\begin{array}{c}
1 & 0 & 0 & 0\\
0 & -2 & 1 & 2
\end{array}\right]$
"""

# ╔═╡ 90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
md"""
Givens Rotations:
"""

# ╔═╡ e0121f81-04db-4c35-af63-867fc073cb16
begin
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
    end

    A = [
        3 4
        0 -2
        0 1
        0 2
    ]

    function givens(A, i, j)
        c, s, r = givens_helper(A[i-1,j], -A[i,j])
        size = maximum((length(A[:,1]), length(A[1,:])))
        G=Matrix{Float64}(I, size, size)
        G[i,j] = s / r
        G[j,j] = c/ r
        G[i,i] = c/ r
        G[j,i] = -s/ r
        return G, G*A
    end
end;

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
Given the QR factors of Aᵀ we can find

$Aᵀ=\frac{1}{3}\left[\begin{array}{cc}
2 & 2\\
1 & 4\\
2 & 5\end{array}\right]$
"""

# ╔═╡ b4469009-a045-4932-bc26-0b8943d8ff83


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
# ╟─3930a39b-c26b-4a42-a3da-8550e2fce0f9
# ╠═9cdf0092-5d3e-4521-8d4a-07c962a0772d
# ╟─64a5e043-7b71-4559-86cc-a316bfc3d668
# ╟─882a0b1c-85d8-475f-be5e-1828faa10196
# ╟─7a949596-5c6e-4faa-8971-154774836859
# ╟─ace59a76-703d-49eb-bea3-479c245ac022
# ╟─fe69bf18-119a-48ea-bee5-1ba5faa0b678
# ╟─90ef3c8a-c8d7-47e6-b7e6-0e10a72a255f
# ╟─e0121f81-04db-4c35-af63-867fc073cb16
# ╟─825e4aca-45d4-42af-8c50-1ab5ee670bad
# ╟─41fc5b67-869c-4b1c-bb5e-cbcd6d2217c1
# ╟─470ae356-7ebe-4d9a-b582-b3796e36d7c1
# ╠═b4469009-a045-4932-bc26-0b8943d8ff83
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
