### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ ffda8bd0-48fa-11ed-25df-f39d3c26627c
md"""
# **MATH303-22S2** - _Assignment 2_

`Due Date`: **Monday 17th October 5pm**

Submission by: Richard Hodges, Henry Hastings

Student IDs: 11139318, 35556669

Code is written in Julia and presented in a Pluto.jl notebook.
"""

# ╔═╡ 31a0905d-fd6a-42eb-aa86-6d71d2bc6d2f
md"""
### Question 1

a)

$\begin{array}{}
\text{maximize }   & -x_1& +& x_2\\ \\
\text{subject to } & x_1 & - & x_2 & + & x_3 & & & & & = & 10 \\
                   & -x_1 & - & x_2 & & & + & x_4 & & & =& -5\\
                   & 2x_1 & + & x_2 & & & & & + & x_5 & = & 40
\end{array}$

Where $x_3,\,x_4,\,x_5$ are slack variables.
"""

# ╔═╡ d54b0145-4e8c-4425-b19c-307b72c9dade
md"""
$A=\left[\begin{array}{}
-1 & 1 & 1 & 0 & 0\\
-1 & -1 & 0 & 1 & 0\\
2 & 1 & 0 & 0 & 1
\end{array}\right]\quad \vec{b}=\left[\begin{array}{}
10\\-5\\40
\end{array}\right]$
"""

# ╔═╡ 34fc3279-53d6-4978-a8a6-3880fb68607d
md"""
b) Let $B=\{3,4,5\}$ as $A_B=I\Rightarrow\exists \;A_B^{-1}$

$\vec{x}_B=A_B^{-1}\vec{b}=\left[\begin{array}{}
10\\ -5 \\ 40 \end{array}\right]$
"""

# ╔═╡ 63abeecf-3e9e-403d-81c1-e0669eea22ae
md"""
Thus $x_4$ leaves $B$ as this corresponds to the most negative element of $\vec{x}_B$
"""

# ╔═╡ 5674fd90-f9e5-460e-b849-e842b9223d9b
md"""
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & \sum\limits_{i=1}^n c_i x_i \\
    \;\;\text{s.t.} & l_j \le \sum\limits_{i=1}^n a_{ij} x_i \le u_j & j = 1 \ldots m \\
    & l_i \le x_i \le u_i & i = 1 \ldots n.
\end{align}
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─ffda8bd0-48fa-11ed-25df-f39d3c26627c
# ╟─31a0905d-fd6a-42eb-aa86-6d71d2bc6d2f
# ╟─d54b0145-4e8c-4425-b19c-307b72c9dade
# ╟─34fc3279-53d6-4978-a8a6-3880fb68607d
# ╠═63abeecf-3e9e-403d-81c1-e0669eea22ae
# ╠═5674fd90-f9e5-460e-b849-e842b9223d9b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
