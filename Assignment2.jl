### A Pluto.jl notebook ###
# v0.19.13

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
**a)**

$\begin{align}{}
\textsf{minimise}   &\quad -x_1 + x_2\\
\textsf{subject to } &\quad x_1-x_2+x_3  &=10 \\
                   &\quad -x_1-x_2+x_4 &=-5\\
                   &\quad 2x_1+x_2+x_5 &=40\\
				   &\quad\underline{x}\geq0
\end{align}$
Where $x_3,\,x_4,\,x_5$ are slack variables.
"""

# ╔═╡ d54b0145-4e8c-4425-b19c-307b72c9dade
md"""
$A=\left[\begin{array}{}
-1 & 1 & 1 & 0 & 0\\
-1 & -1 & 0 & 1 & 0\\
2 & 1 & 0 & 0 & 1
\end{array}\right]\quad \underline{b}=\left[\begin{array}{}
10\\-5\\40
\end{array}\right]$
"""

# ╔═╡ 34fc3279-53d6-4978-a8a6-3880fb68607d
md"""
**b)** Let $B=\{3,4,5\}$ as $A_B=I\Rightarrow\exists \;A_B^{-1}$
$\vec{x}_B=A_B^{-1}\vec{b}=\left[\begin{array}{}
10\\ -5 \\ 40 \end{array}\right]$
"""

# ╔═╡ 63abeecf-3e9e-403d-81c1-e0669eea22ae
md"""
Thus $x_4$ leaves $B$ as this corresponds to the most negative element of $\vec{x}_B$
"""

# ╔═╡ 19d0433a-6092-467f-a56f-b86730ad4631
md"""
$\underline{a}_6=\underline{b}-\left[\begin{array}{} 1 \\ 0 \\ 0 \end{array}\right] -
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] -
\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right] +
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] =
\left[\begin{array}{} 9 \\ -5 \\ 39 \end{array}\right]$
"""

# ╔═╡ c02f780a-c965-4c35-a7ef-3f88fd251d3f
md"""
So $B^\#=\{3,\,5,\,6\}$ is a B.F.S for the phase 1 problem for this L.P.
"""

# ╔═╡ 0d7d1da3-3848-44b9-8d04-0ec52673a825
md"""
##### Phase 1 problem:

$\begin{align}{}
\textsf{minimise} & \quad x_6\\
\textsf{subject to:} &\quad [A\,|\,a_6]\left[\begin{array}{}\vec{x}\\x_6\end{array}\right]=\vec{b}\\
&\quad x_6,\,\underline{x}\geq0
\end{align}$
"""

# ╔═╡ 41b9a441-1338-4551-a0d7-ee6ea3878080
md"""
**c)** $B=\{3,\,5,\,6\}$ by b) thus $N=\{1,\,2,\,4\}$

$\begin{align}
\Rightarrow\vec{r} & = \vec{0} -\left[\begin{array}{} 
-1 & -1 & 2 \\ 1 & -1 & 1 \\ 0 & 1 & 0
\end{array}\right]
\left[\begin{array}{} 
1 & 0 & 0 \\ 1.8 & 7.2 & -0.2\\ 0 & 1 & 0
\end{array}\right]
\left[\begin{array}{} 
0 \\ 0 \\ 1
\end{array}\right]\\
& = \left[\begin{array}{} 
-1 & -1 & 2 \\ 1 & -1 & 1 \\ 0 & 1 & 0
\end{array}\right]
\left[\begin{array}{} 
0 \\ -0.2 \\ 0
\end{array}\right]=
\left[\begin{array}{} 
-0.2 \\ -0.2 \\ 0.2
\end{array}\right]
\end{align}$
"""

# ╔═╡ a96c5a04-5595-4a2a-befb-f3e77c59bd29
md"""
Choose $x_1$ to enter as it has a lower index.

$\hat{\underline{b}}=A_B^{-1}\underline{b}=\left[\begin{array}{} 
1 \\ 1 \\ 1\end{array}\right]\quad\quad \underline{d}=-A_B^{-1}\underline{a}_1=\left[\begin{array}{} 2.8 \\ 5 \\ -0.2 \end{array}\right]$
"""

# ╔═╡ 29512407-ae19-4115-853b-ae3a1c864004
md"""
Therefore $x_6$ must be the leavign variable

$\Rightarrow B=\{1,\,3,\,5\},\quad N=\{2,\,4,\,6,\,\}$
$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & 0 & 0 \\ 2 & 0 & 1
\end{array}\right]
\Rightarrow A_B^{-1}=\left[\begin{array}{} 
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 & 9 \\ -1 & 1 & -5 \\ 1 & 0 & 39
\end{array}\right]$
"""

# ╔═╡ b0c98553-2fcb-4c69-b23c-fa922ada1123
md"""
$\Rightarrow \underline{r}=\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right]
-A_N^TA_B^{-T}\underline{0}=\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right]$
"""

# ╔═╡ ed2283c2-a273-4922-9601-4246a716cddf
md"""
So $B=\{1,\,3,\,5\}$ is an optimal basis for the phase 1 problem and thus a B.F.S. for the original L.P.
"""

# ╔═╡ e0c74398-2c65-4804-8936-02135659d23c
md"""
**d)** 

$\Rightarrow B=\{1,\,3,\,5\},\quad N=\{2,\,4\}$
$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & 0 & 0 \\ 2 & 0 & 1
\end{array}\right]
\quad A_B^{-1}=\left[\begin{array}{} 
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 \\ -1 & 1  \\ 1 & 0 
\end{array}\right]$
"""

# ╔═╡ b7ad32ee-dbc3-4b8f-8af5-fa4f22a7182b
md"""
$\begin{align}\Rightarrow\underline{r}&=\left[\begin{array}{} -1 \\ 0 \end{array}\right]-\left[\begin{array}{} 
1 & -1 & 1 \\ 0 & 1 & 0
\end{array}\right]\left[\begin{array}{} 
0 & 1 & 0 \\ -1 & -1 & 2 \\ 0 & 0 & 1
\end{array}\right]\left[\begin{array}{} -1 \\ 0 \\ 0 \end{array}\right]\\
&=\left[\begin{array}{} -1 \\ 0 \end{array}\right]-\left[\begin{array}{}
1 & -1 & 1 \\ 0 & 1 & 0 \end{array}\right]\left[\begin{array}{}
0 \\ 1 \\ 0 \end{array}\right]\\
& = \left[\begin{array}{} -1 \\ 0 \end{array}\right] - \left[\begin{array}{} 
-1 \\ 1 \end{array}\right]\\
&=\left[\begin{array}{} 0 \\ -1\end{array}\right] 
\end{align}$
"""

# ╔═╡ 826a61ba-3442-4426-b889-2392f0272895
md"""
So $x_4$ is the entering variable
"""

# ╔═╡ 4bfa234f-a4d7-45c4-a1ef-8908f6e09b18
md"""
$\hat{\underline{b}}=\left[\begin{array}{}
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\left[\begin{array}{} 
10 \\ -5 \\ 40
\end{array}\right]= \left[\begin{array}{} 5 \\ 15 \\ 30 \end{array}\right],
\quad\hat{\underline{d}}=\left[\begin{array}{} 
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right],
\quad \underline{a}_4=\left[\begin{array}{} 1 \\ 1 \\ -2 \end{array}\right]$
"""

# ╔═╡ aa19cdab-3f98-46a9-8735-37cac3fa4f29
md"""
Giving us $x_5$ as the leaving variable
"""

# ╔═╡ 4c8dc4d8-6163-432c-8a29-6da48c154343
md"""
$\Rightarrow B=\{1,\,3,\,4\}, \quad N=\{2,\,5\}$

$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & 0 & 1 \\ 2 & 0 & 0
\end{array}\right]
\quad A_B^{-1}=\left[\begin{array}{} 
0 & 0 & 0.5 \\ 1 & 0 & 0.5 \\ 0 & 1 & 0.5
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 \\ -1 & 0  \\ 1 & 1
\end{array}\right]$
"""

# ╔═╡ c15fd4d4-504a-4aed-b7c6-f4ef17720186
md"""
$\begin{align}
\Rightarrow\underline{r}&=\left[\begin{array}{} -1 \\ 0 \end{array}\right] -
\left[\begin{array}{} 
1 & -1 & 1\\ 0 & 0 & 1
\end{array}\right]\left[\begin{array}{} 
0 & 1 & 0 \\ 0 & 0 & 1 \\ 0.5 & 0.5 & 0.5
\end{array}\right]\left[\begin{array}{} -1 \\ 0 \\ 0 \end{array}\right]\\
&=\left[\begin{array}{} -1 \\ 0\end{array}\right] - \left[\begin{array}{} 
-0.5 \\ -0.5 \end{array}\right] = \left[\begin{array}{} 
0.5 \\ 0.5 \end{array}\right]
\end{align}$
"""

# ╔═╡ 8403d3fd-8606-4401-9e63-23182ec8141a
md"""
So $x_2$ is the entering variable"""

# ╔═╡ 31b11798-95ba-4d82-8ac3-e311e2c6e367
md"""
$\hat{\underline{b}}=A_B^{-1}\underline{b}=\left[\begin{array}{} 
0 & 0 & 0.5 \\ 1 & 0 & 0.5 \\ 0 & 1 & 0.5
\end{array}\right]\left[\begin{array}{} 
10 \\ -5 \\ 40
\end{array}\right]=\left[\begin{array}{} 
20 \\ 30 \\ 15
\end{array}\right],\quad\underline{d}=A_B^{-1}\underline{a}_2= \left[\begin{array}{} 
0 & 0 & 0.5 \\ 1 & 0 & 0.5 \\ 0 & 1 & 0.5
\end{array}\right]\left[\begin{array}{} 
1 \\ -1 \\ 1
\end{array}\right] = \left[\begin{array}{} 
-0.5 \\ -1.5 \\ 0.5
\end{array}\right]$
"""

# ╔═╡ ca23707c-1e54-4c2e-be2c-f84f35ac8605
md"""
$\frac{-\hat{b}_1}{d_1}=40\quad\quad\frac{-\hat{b}_2}{d_2}=20$
"""

# ╔═╡ 0117f75d-69b9-4da9-a84f-6738e3a3f1c4
md"""
As $\frac{-\hat{b}_2}{d_2}<\frac{-\hat{b}_1}{d_1}$ we get $a_3$ as the leaving variable.
"""

# ╔═╡ f3213504-148f-4b22-b440-c657b93b0534
md"""
$\Rightarrow B=\{1,\,2,\,4\}, \quad N=\{3,\,5\}$

$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & -1 & 1 \\ 2 & 1 & 0
\end{array}\right]
\quad A_B^{-1}=\left[\begin{array}{} 
\frac{-1}{3} & 0 & \frac{1}{3} \\ \frac{2}{3} & 0 & \frac{1}{3} \\ \frac{1}{3} & 1 & \frac{2}{3}
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 \\ 0 & 0  \\ 0 & 1
\end{array}\right]$
"""

# ╔═╡ efbae989-0099-4aa0-9601-3b97e4d9c83d
md"""
$\begin{align}
\Rightarrow\underline{r} & = \underline{0} -\left[\begin{array}{} 
1 & 0 & 0 \\ 0 & 0 & 1
\end{array}\right]\left[\begin{array}{} 
\frac{-1}{3} & \frac{2}{3} & \frac{1}{3} \\
0 & 0 & 1\\
\frac{1}{3} & \frac{1}{3} & \frac{2}{3}
\end{array}\right]\left[\begin{array}{} 
-1 \\ -1 \\ 0
\end{array}\right]\\
&= \left[\begin{array}{}\frac{1}{3}\\\frac{2}{3} \end{array}\right]
\end{align}$
"""

# ╔═╡ 5823e3d2-5dfb-41c1-a7c7-51062f9ed976
md"""
As $\underline{r}\geq\underline{0}$ we have an optimal basis

$\underline{x}^*_B = A_B^{-1} \underline{b} = \left[\begin{array}{} 10 \\ 20 \\ 25 \end{array}\right] \quad \underline{x}_N^* = \underline{0}$
"""

# ╔═╡ 6ea3d17f-e9ee-438c-ae5b-f5c46eac3a06
md"""
$\underline{x}^*=\left[\begin{array}{} 
10 \\ 20 \\ 0 \\ 25 \\ 0
\end{array}\right]$
"""

# ╔═╡ 2d838108-8b40-481a-98a4-8c1e7e907b43
md"""
This gives minimal value of 

$-x_1-x_2=-30$
"""

# ╔═╡ adf6d7e0-68d1-4e76-ae4b-68129219b272
md"""Maximal value of

$x_1+x_2=30$
"""

# ╔═╡ d8b655a5-79e4-4cc8-82e6-288add3554a5
md"""
The solution ot the L.P. is 

$x_1=10,\;x_2=20;\;x_1+x_2=30$
"""

# ╔═╡ a3acc4ae-ad02-45ca-98ed-f42f7c994ea1
md"""All inverses for this question were found using Symbolab."""

# ╔═╡ 6d83b238-643e-4ad5-b7c5-637d5ebbbeef


# ╔═╡ 5674fd90-f9e5-460e-b849-e842b9223d9b
md"""

\left[\begin{array}{} \end{array}\right]

```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & \sum\limits_{i=1}^n c_i x_i \\
    \;\;\text{s.t.} & l_j \le \sum\limits_{i=1}^n a_{ij} x_i \le u_j & j = 1 \ldots m \\
    & l_i \le x_i \le u_i & i = 1 \ldots n.
\end{align}
```
"""

# ╔═╡ Cell order:
# ╟─ffda8bd0-48fa-11ed-25df-f39d3c26627c
# ╟─31a0905d-fd6a-42eb-aa86-6d71d2bc6d2f
# ╟─d54b0145-4e8c-4425-b19c-307b72c9dade
# ╟─34fc3279-53d6-4978-a8a6-3880fb68607d
# ╟─63abeecf-3e9e-403d-81c1-e0669eea22ae
# ╟─19d0433a-6092-467f-a56f-b86730ad4631
# ╟─c02f780a-c965-4c35-a7ef-3f88fd251d3f
# ╟─0d7d1da3-3848-44b9-8d04-0ec52673a825
# ╟─41b9a441-1338-4551-a0d7-ee6ea3878080
# ╟─a96c5a04-5595-4a2a-befb-f3e77c59bd29
# ╟─29512407-ae19-4115-853b-ae3a1c864004
# ╟─b0c98553-2fcb-4c69-b23c-fa922ada1123
# ╟─ed2283c2-a273-4922-9601-4246a716cddf
# ╟─e0c74398-2c65-4804-8936-02135659d23c
# ╟─b7ad32ee-dbc3-4b8f-8af5-fa4f22a7182b
# ╟─826a61ba-3442-4426-b889-2392f0272895
# ╟─4bfa234f-a4d7-45c4-a1ef-8908f6e09b18
# ╟─aa19cdab-3f98-46a9-8735-37cac3fa4f29
# ╟─4c8dc4d8-6163-432c-8a29-6da48c154343
# ╟─c15fd4d4-504a-4aed-b7c6-f4ef17720186
# ╟─8403d3fd-8606-4401-9e63-23182ec8141a
# ╟─31b11798-95ba-4d82-8ac3-e311e2c6e367
# ╟─ca23707c-1e54-4c2e-be2c-f84f35ac8605
# ╟─0117f75d-69b9-4da9-a84f-6738e3a3f1c4
# ╟─f3213504-148f-4b22-b440-c657b93b0534
# ╟─efbae989-0099-4aa0-9601-3b97e4d9c83d
# ╟─5823e3d2-5dfb-41c1-a7c7-51062f9ed976
# ╟─6ea3d17f-e9ee-438c-ae5b-f5c46eac3a06
# ╟─2d838108-8b40-481a-98a4-8c1e7e907b43
# ╟─adf6d7e0-68d1-4e76-ae4b-68129219b272
# ╟─d8b655a5-79e4-4cc8-82e6-288add3554a5
# ╠═a3acc4ae-ad02-45ca-98ed-f42f7c994ea1
# ╠═6d83b238-643e-4ad5-b7c5-637d5ebbbeef
# ╠═5674fd90-f9e5-460e-b849-e842b9223d9b
