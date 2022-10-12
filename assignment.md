# **MATH303-22S2** - _Assignment 2_

`Due Date`: **Monday 17th October 5pm**

`Submission by`: Richard Hodges, Henry Hastings

`Student IDs`: 11139318, 35556669

### Question 1
**a)**

$$\begin{align*}{}
\textsf{minimise}   &\quad -x_1 + x_2\\
\textsf{subject to } &\quad x_1-x_2+x_3  &=10 \\
                   &\quad -x_1-x_2+x_4 &=-5\\
                   &\quad 2x_1+x_2+x_5 &=40\\
				   &\quad\underline{x}\geq0
\end{align*}$$

Where $x_3,\,x_4,\,x_5$ are slack variables.

$$A=\left[\begin{array}{}
-1 & 1 & 1 & 0 & 0\\
-1 & -1 & 0 & 1 & 0\\
2 & 1 & 0 & 0 & 1
\end{array}\right]\quad \underline{b}=\left[\begin{array}{}
10\\-5\\40
\end{array}\right]$$

**b)** Let $B=\{3,4,5\}$ as $A_B=I\Rightarrow\exists \;A_B^{-1}$

$$\vec{x}_B=A_B^{-1}\vec{b}=\left[\begin{array}{}
10\\ -5 \\ 40 \end{array}\right]$$

Thus $x_4$ leaves $B$ as this corresponds to the most negative element of $\vec{x}_B$

**Explain artificial variable maybe??**

$$\underline{a}_6=\underline{b}-\left[\begin{array}{} 1 \\ 0 \\ 0 \end{array}\right] -
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] -
\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right] +
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] =
\left[\begin{array}{} 9 \\ -5 \\ 39 \end{array}\right]$$

So $B^\#=\{3,\,5,\,6\}$ is a B.F.S for the phase 1 problem for this L.P.

**Phase 1 problem:**

$$\begin{align*}{}
\textsf{minimise} & \quad x_6\\
\textsf{subject to:} &\quad [A\,|\,a_6]\left[\begin{array}{}\vec{x}\\x_6\end{array}\right]=\vec{b}\\
&\quad x_6,\,\underline{x}\geq0
\end{align*}$$

