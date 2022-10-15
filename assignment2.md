# **MATH303-22S2** - _Assignment 2_

`Due Date:` **Monday 17th October 5pm**

`Submission by:` Richard Hodges, Henry Hastings

`Student IDs:` 11139318, 35556669



## **Question 1**
**a)**

$$\begin{align*}{}
\operatornamewithlimits{\textsf{minimise}}_{\underline{x}}   &\quad -x_1 + x_2\\
\textsf{subject to } &\quad x_1-x_2+x_3  &=10 \\
                   &\quad -x_1-x_2+x_4 &=-5\\
                   &\quad 2x_1+x_2+x_5 &=40\\
				   &\quad\underline{x}\geq0
\end{align*}$$

Where $x_3,\,x_4,\,x_5$ are slack variables.


$$
A=\left[\begin{array}{}
-1 & 1 & 1 & 0 & 0\\
-1 & -1 & 0 & 1 & 0\\
2 & 1 & 0 & 0 & 1
\end{array}\right]\quad 
\underline{b}=\left[\begin{array}{}
10\\-5\\40
\end{array}\right]
$$

-------------

**b)** Let $B=\{3,4,5\}$ as $A_B=I\Rightarrow\exists \;A_B^{-1}=I$


$$\underline{x}_B=A_B^{-1}\underline{b}=\left[\begin{array}{}
10\\ -5 \\ 40 \end{array}\right]$$


Thus $x_4$ leaves $B$ as this corresponds to the most negative element of $\underline{x}_B$

We can replace $x_4$ in $B$  with a singular artificial variable as seen in part 8.1 of the lecture notes. The artificial variable, $\underline{a}_6$, which is adjoined to the matrix $A$, is obtained in the following way:




$$\underline{a}_6=\underline{b}-\left[\begin{array}{} 1 \\ 0 \\ 0 \end{array}\right] -
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] -
\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right] +
\left[\begin{array}{} 0 \\ 1 \\ 0\end{array}\right] =
\left[\begin{array}{} 9 \\ -5 \\ 39 \end{array}\right]$$



From part 8.1 of the lecture notes we see that $B^\#=\{3,\,5,\,6\}$ is a B.F.S of the following L.P: 

$$[A\,|\,a_6]\left[\begin{array}{}\vec{x}\\x_6\end{array}\right]=\underline{b}$$


This represents the constraints of our phase 1 problem, as outlined in part 8.1 of the lecture notes. The phase 1 problem requires us to minimise $x_6$ given this constraint, removing it from $B$ which yields a B.F.S of our original L.P. 




**Phase 1 problem:**

$$\begin{align*}{}
\operatornamewithlimits{\textsf{minimise}}_{\underline{x}}  & \quad x_6\\
\textsf{subject to:} & \quad [A\,|\,a_6]\left[\begin{array}{}\vec{x}\\x_6\end{array}\right]=\underline{b}\\
&\quad x_6,\,\underline{x}\geq0
\end{align*}$$


---------------------


**c)** We solve this L.P using the simplex method as described in Algorithm 1 of the week 7 notes.


$B=\{3,\,5,\,6\}$ by **b)** thus $N=\{1,\,2,\,4\}$


$$\Rightarrow A_B= \left[\begin{array}{} 
1 & 0 & 9 \\ 0 & 0 & -5 \\ 0 & 1 & 39
\end{array}\right]
\Rightarrow A_B^{-1}=\left[\begin{array}{} 
1 & 1.8 & 0 \\ 0 & 7.2 & 1 \\ 0 & -0.2 & 0
\end{array}\right],\quad A_N=\left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & -1 & 1 \\ 2 & 1 & 0
\end{array}\right]$$ 


$$\underline{c}_N = \underline{0},\quad\underline{c}_B= \left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right]$$


We define $\underline{r}= \underline{c}_N -A_N^TA_B^{-T}\underline{c}_B$


$$\begin{align*}
\Rightarrow\underline{r} & = \vec{0} -\left[\begin{array}{} 
-1 & -1 & 2 \\ 1 & -1 & 1 \\ 0 & 1 & 0
\end{array}\right]
\left[\begin{array}{} 
1 & 0 & 0 \\ 1.8 & 7.2 & -0.2\\ 0 & 1 & 0
\end{array}\right]
\left[\begin{array}{} 
0 \\ 0 \\ 1
\end{array}\right]=\left[\begin{array}{} 
-0.2 \\ -0.2 \\ 0.2
\end{array}\right]
\end{align*}$$


Choose $x_1$ to enter as it has a lower index of the variables corresponding to the negative values of $\underline{r}$, using Bland's anti-cycling rules.

$$\hat{\underline{b}}=A_B^{-1}\underline{b}=\left[\begin{array}{} 
1 \\ 1 \\ 1\end{array}\right]\quad\quad \underline{d}=-A_B^{-1}\underline{a}_1=\left[\begin{array}{} 2.8 \\ 5 \\ -0.2 \end{array}\right]$$


 $x_6$ must be the leaving variable as it represents the only negative value in $\underline{d}$.


$$\Rightarrow B=\{1,\,3,\,5\},\quad N=\{2,\,4,\,6\}$$


$$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & 0 & 0 \\ 2 & 0 & 1
\end{array}\right]
\Rightarrow A_B^{-1}=\left[\begin{array}{} 
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 & 9 \\ -1 & 1 & -5 \\ 1 & 0 & 39
\end{array}\right]$$


$$\underline{c}_N= \left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right],\quad\underline{c}_B = \underline{0}$$


$$\Rightarrow \underline{r}=\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right]
-A_N^TA_B^{-T}\underline{0}=\left[\begin{array}{} 0 \\ 0 \\ 1\end{array}\right]$$


As $\underline{r}≥0$, $B=\{1,\,3,\,5\}$ is an optimal basis for the phase 1 problem and thus a B.F.S. for the original L.P.

-------

**d)** We solve this L.P using the simplex method as described in Algorithm 1 of the week 7 notes


$$\Rightarrow B=\{1,\,3,\,5\},\quad N=\{2,\,4\}$$


$$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & 0 & 0 \\ 2 & 0 & 1
\end{array}\right]
\quad A_B^{-1}=\left[\begin{array}{} 
0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\quad A_N=\left[\begin{array}{} 
1 & 0 \\ -1 & 1  \\ 1 & 0 
\end{array}\right]$$


$$\underline{c}_N= \left[\begin{array}{} -1 \\ 0 \end{array}\right],\quad\underline{c}_B = \left[\begin{array}{} -1 \\ 0 \\ 0 \end{array}\right]\\$$


$$
\begin{align*}
    \Rightarrow\underline{r} &= \left[\begin{array}{} -1 \\ 0 \end{array}\right] - \left[\begin{array}{} 
    1 & -1 & 1 \\ 0 & 1 & 0
    \end{array}\right]\left[\begin{array}{} 
    0 & 1 & 0 \\ -1 & -1 & 2 \\ 0 & 0 & 1
    \end{array}\right]\left[\begin{array}{} 
    -1 \\ 0 \\ 0 
    \end{array}\right]\\

    &= \left[\begin{array}{} -1 \\ 0 \end{array}\right]-\left[\begin{array}{}
    1 & -1 & 1 \\ 0 & 1 & 0 \end{array}\right]\left[\begin{array}{}
    0 \\ 1 \\ 0 \end{array}\right]\\

    &= \left[\begin{array}{} -1 \\ 0 \end{array}\right] - \left[\begin{array}{} 
    -1 \\ 1 \end{array}\right]\\
    &=\left[\begin{array}{} 0 \\ -1\end{array}\right] 
\end{align*}
$$


So $x_4$ is the entering variable


$$
\hat{\underline{b}}=\left[\begin{array}{}
    0 & -1 & 0 \\ 
    1 & -1 & 0 \\ 
    0 & 2 & 1
\end{array}\right]\left[\begin{array}{} 
    10 \\ -5 \\ 40
\end{array}\right]= \left[\begin{array}{} 
    5 \\ 15 \\ 30 
\end{array}\right],
\quad{\underline{d}}=\left[\begin{array}{} 
    0 & -1 & 0 \\ 1 & -1 & 0 \\ 0 & 2 & 1
\end{array}\right]\underline{a}_4=\left[\begin{array}{} 
    1 \\ 1 \\ -2 
\end{array}\right]$$


Giving us $x_5$ as the leaving variable, as this corresponds to the only negative entry in $\underline{d}$


$$\Rightarrow B=\{1,\,3,\,4\}, \quad N=\{2,\,5\}$$


$$\Rightarrow A_B= \left[\begin{array}{} 
    -1 & 1 & 0 \\ 
    -1 & 0 & 1 \\ 
    2 & 0 & 0
\end{array}\right]
\quad A_B^{-1}=\left[\begin{array}{} 
    0 & 0 & 0.5 \\ 
    1 & 0 & 0.5 \\ 
    0 & 1 & 0.5
\end{array}\right]\quad A_N=\left[\begin{array}{} 
    1 & 0 \\ 
    -1 & 0  \\ 
    1 & 1
\end{array}\right]$$


$$
\underline{c}_N= \left[\begin{array}{} 
    -1 \\ 0 
\end{array}\right],
\quad \underline{c}_B = \left[\begin{array}{} 
    -1 \\ 0 \\ 0 
\end{array}\right]
$$


$$\begin{align*}
\Rightarrow\underline{r}&=\left[\begin{array}{} -1 \\ 0 \end{array}\right] -
\left[\begin{array}{} 
1 & -1 & 1\\ 0 & 0 & 1
\end{array}\right]\left[\begin{array}{} 
0 & 1 & 0 \\ 0 & 0 & 1 \\ 0.5 & 0.5 & 0.5
\end{array}\right]\left[\begin{array}{} -1 \\ 0 \\ 0 \end{array}\right]\\
&=\left[\begin{array}{} -1 \\ 0\end{array}\right] - \left[\begin{array}{} 
-0.5 \\ -0.5 \end{array}\right] = \left[\begin{array}{} 
0.5 \\ 0.5 \end{array}\right]
\end{align*}$$


So $x_2$ is the entering variable


$$
\hat{\underline{b}}=A_B^{-1}\underline{b}=\left[\begin{array}{} 
    0 & 0 & 0.5 \\ 
    1 & 0 & 0.5 \\ 
    0 & 1 & 0.5
\end{array}\right]\left[\begin{array}{} 
    10 \\ -5 \\ 40
\end{array}\right]=\left[\begin{array}{} 
    20 \\ 30 \\ 15
\end{array}\right]$$


$$\underline{d}=A_B^{-1}\underline{a}_2= \left[\begin{array}{} 
    0 & 0 & 0.5 \\ 
    1 & 0 & 0.5 \\ 
    0 & 1 & 0.5
\end{array}\right]\left[\begin{array}{} 
    1 \\ -1 \\ 1
\end{array}\right] = \left[\begin{array}{} 
    -0.5 \\ -1.5 \\ 0.5
\end{array}\right]$$


$$\frac{-\hat{b}_1}{d_1}=40\quad\quad\frac{-\hat{b}_2}{d_2}=20$$


As $\frac{-\hat{b}_2}{d_2}<\frac{-\hat{b}_1}{d_1}$ we get $x_3$ as the leaving variable, which corresponds to $\hat{b}_2$ and ${d_2}$.


$$\Rightarrow B=\{1,\,2,\,4\}, \quad N=\{3,\,5\}$$


$$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 1 & 0 \\ -1 & -1 & 1 \\ 2 & 1 & 0
\end{array}\right]
\Rightarrow\quad A_B^{-1}=\left[\begin{array}{} 
\frac{-1}{3} & 0 & \frac{1}{3} \\ \frac{2}{3} & 0 & \frac{1}{3} \\ \frac{1}{3} & 1 & \frac{2}{3}
\end{array}\right],\quad A_N=\left[\begin{array}{} 
1 & 0 \\ 0 & 0  \\ 0 & 1
\end{array}\right]$$


$$\underline{c}_N= \underline{0},\quad\underline{c}_B = \left[\begin{array}{} -1 \\ -1 \\ 0 \end{array}\right]\\$$


$$\begin{align*}
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
\end{align*}$$


As $\underline{r}\geq\underline{0}$ we have an optimal basis


$$\underline{x}^*_B = A_B^{-1} \underline{b} = \left[\begin{array}{} 10 \\ 20 \\ 25 \end{array}\right] \quad \underline{x}_N^* = \underline{0}$$


$$\underline{x}^*=\left[\begin{array}{} 
10 \\ 20 \\ 0 \\ 25 \\ 0
\end{array}\right]$$


This gives minimal value of 


$$-x_1-x_2=-30$$


and a maximal value of


$$x_1+x_2=30$$


The solution ot the L.P. is 

$x_1=10,\;x_2=20;\;x_1+x_2=30$


All inverses for this question were found using Symbolab.


## **Question 2**

In standard form our constraints are:


$$\begin{align*}
-2x_1 - x_2 + x_3 &= -70\\
x_1 + x_2 + x_4 &= 40\\
-x_1 - 3x_2 + x_5 &= -90\\
\underline{x}\geq 0
\end{align*}$$


Where $x_3,\,x_4,\,x_5$ are slack variables.


$$\Rightarrow A=\left[\begin{array}{} 
-2 & -1 & 1 & 0 & 0\\
1 & 1 & 0 & 1 & 0 \\
-1 & -3 & 0 & 0 & 1
\end{array}\right] \quad\quad \underline{b} = \left[\begin{array}{} 
-70 \\ 40 \\ -90
\end{array}\right]$$


Let $B=\{3,\,4,\,5\}$ as $\exists\;A_B^{-1}$ as $A_B=I_3$


$$\underline{x}_B = \left[\begin{array}{} -70\\ 40 \\ -90 \end{array}\right] = \underline{b}$$


We can replace $x_5$ in $B$ as it corresponds to the most negative entry of $\underline{x}_B$. We replace it with a singular artificial variable as seen in part 8.1 of the lecture notes. The artificial variable, ${x}_6$, which is adjoined to the matrix $A$, is obtained in the following way:


$$\underline{a}_6 = \underline{b} - \underline{a}_3 - \underline{a}_4 + \underline{a}_5 = \left[\begin{array}{} -71 \\ 39 \\ -90 \end{array}\right]$$


This yields a phase 1 problem with $B=\{3,\,4,\,6\}$:


$$\begin{align*}
\operatornamewithlimits{\textsf{minimise}}_{\underline{x}}  \quad & x_6\\
\textsf{subject to} \quad & [A\,|\,\underline{a}_6]\left[\begin{array}{} \underline{x} \\ x_6\end{array}\right] = \underline{b}\\
& \underline{x},\,x_6\geq \underline{0}
\end{align*}$$


We solve this L.P using the simplex method as described in Algorithm 1 of the week 7 notes


$$B=\{3,\,4,\,6\}\quad N=\{1,\,2,\,5\}$$
$$A_B= \left[\begin{array}{} 
1 & 0 & -71 \\  0 & 1 & 39 \\ 0 & 0 & -90
\end{array}\right] \quad A_B^{-1} = \left[\begin{array}{} 
1 & 0 & \frac{-71}{90} \\ 0 & 1 & \frac{39}{90} \\ 0 & 0 & \frac{-1}{90}
\end{array}\right] \quad A_N = \left[\begin{array}{} 
-2 & -1 & 0 \\ 1 & 1 & 0 \\ -1 & -3 & 1
\end{array}\right]$$


$$\underline{c}_N= \underline{0},\quad\underline{c}_B = \left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right]\\$$


$$\begin{align*}
\underline{r} &= \underline{0} -\left[\begin{array}{} 
-2 & -1 & 0 \\ 1 & 1 & 0 \\ -1 & -3 & 1
\end{array}\right] \left[\begin{array}{} 
1 & 0 & 0 \\ 0 & 1 & 0 \\ \frac{-71}{90} & \frac{39}{90} & \frac{-1}{90}
\end{array}\right]\left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right] \\
&=\left[\begin{array}{} \frac{-1}{90} \\ \frac{-3}{90} \\ \frac{1}{90} \end{array}\right]
\end{align*}$$


So $x_2$ is the entering variable.


$$\hat{\underline{b}} = A_B^{-1}\underline{b} = \left[\begin{array}{} 1 \\ 1 \\ 1 \end{array}\right] \quad \underline{d} = -A_B^{-1} \underline{a}_2 =\left[\begin{array}{} \frac{-123}{90} \\ \frac{27}{90} \\ \frac{3}{90} \end{array}\right]$$


$$\frac{-\hat{b}_1}{d_1} =\frac{90}{123}\quad \frac{-\hat{b}_3}{d_3} =\frac{90}{3}$$

So $x_3$ is the leaving variable, as $\frac{-\hat{b}_1}{d_1}<\frac{-\hat{b}_3}{d_3}$ and $x_3$ corresponds to $\hat{b}_1$ and ${d_1}$.


$$\Rightarrow B=\{2,\,4,\,6\}\quad N=\{1,\,3,\,5\}$$

$$\Rightarrow A_B= \left[\begin{array}{} 
-1 & 0 & -71 \\ 1 & 1 & 39 \\ -3 & 0 & -90
\end{array}\right]\Rightarrow \quad A_B^{-1} = \left[\begin{array}{} 
\frac{30}{41} & 0 & \frac{-71}{90} \\ \frac{9}{41} & 1 & \frac{39}{90} \\ \frac{-1}{41} & 0 & \frac{-1}{90}
\end{array}\right], \quad A_N = \left[\begin{array}{} 
-2 & 1 & 0 \\ 1 & 0 & 0 \\ -1 & 0 & 1
\end{array}\right]$$

$$\underline{c}_N= \underline{0},\quad\underline{c}_B = \left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right]\\$$


$$\begin{align*}
\underline{r} &= \underline{0} -\left[\begin{array}{} 
-2 & 1 & 0 \\ 1 & 0 & 0 \\ -1 & 0 & 1
\end{array}\right] \left[\begin{array}{} 
\frac{30}{41} & \frac{9}{41} & \frac{-1}{41} \\ 0 & 1 & 0 \\ \frac{-71}{90} & \frac{39}{90} & \frac{-1}{90}
\end{array}\right]\left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right] \\
&=\left[\begin{array}{} \frac{-5}{123} \\ \frac{1}{41} \\ \frac{-1}{123} \end{array}\right]
\end{align*}$$


Therefore $x_1$ is the entering variable.


$$\hat{\underline{b}} = A_B^{-1}\underline{b} = \left[\begin{array}{} \frac{30}{41} \\ \frac{50}{41} \\ \frac{40}{41} \end{array}\right] \quad \underline{d} = -A_B^{-1} \underline{a}_2 =\left[\begin{array}{} \frac{109}{123} \\ \frac{-37}{123} \\ \frac{-5}{123} \end{array}\right]$$

$$\frac{-\hat{b}_2}{d_2} =\frac{150}{37}\quad \frac{-\hat{b}_3}{d_3} =\frac{120}{5}$$

So $x_4$ is the leaving variable, as $\frac{-\hat{b}_2}{d_2}<\frac{-\hat{b}_3}{d_3}$ and $x_4$ corresponds to $\hat{b}_2$ and ${d_2}$

$$\Rightarrow B=\{1,\,2,\,6\}\quad N=\{3,\,4,\,5\}$$
$$\Rightarrow A_B= \left[\begin{array}{} 
-2 & -1 & -71 \\ 1 & 1 & 39 \\ -1 & -3 & -90
\end{array}\right]\Rightarrow \quad A_B^{-1} = \left[\begin{array}{} 
\frac{27}{37} & \frac{123}{37} & \frac{32}{37} \\
\frac{51}{37} & \frac{109}{37} & \frac{7}{37} \\
\frac{-2}{37} & \frac{-5}{37} & \frac{-1}{37}
\end{array}\right], \quad A_N = \left[\begin{array}{} 
1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1
\end{array}\right]$$


$$\underline{c}_N= \underline{0},\quad\underline{c}_B = \left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right]\\$$


$$\begin{align*}
\underline{r} &= \underline{0} -\left[\begin{array}{} 
\frac{27}{37} & \frac{51}{37} & \frac{-2}{37} \\
\frac{123}{37} & \frac{109}{37} & \frac{-5}{37} \\
\frac{32}{37} & \frac{7}{37} & \frac{-1}{37}
\end{array}\right] \left[\begin{array}{} 0 \\ 0 \\ 1 \end{array}\right] \\
&=\left[\begin{array}{} \frac{2}{37} \\ \frac{5}{37} \\ \frac{1}{37} \end{array}\right] \geq \underline{0}
\end{align*}$$


Thus, the optimal basis is $B=\{2,\,4,\,6\}$

The optimal basis still contains the artificial variable. This shows us that there is no basic feasible solution to our set of constraints without the artificial variable. Thus, this shows there are no vectors satisfying our constraints.

All inverses for this question were found using Symbolab.

---------------

## **Question 3**

**a)**


$$A=\left[\begin{array}{} 
1948 & 1 \\
\vdots & \vdots \\
2020 & 1
\end{array}\right] \quad \underline{b} = \left[\begin{array}{} 
10.3 \\ \vdots \\ 9.80
\end{array}\right] \underline{x} = \left[\begin{array}{} 
x_1 \\ x_2 
\end{array}\right] \quad \underline{t}= \left[\begin{array}{} 
t_1 \\ \vdots \\ t_{19}
\end{array}\right]$$


Using this we can rewrite our constraints in the form


$$\left[\begin{array}{}\begin{array}{} 
A \\ -A \end{array} & \bigg | & \begin{array}{} -I_{19} \\ -I_{19} 
\end{array} \end{array}\right]\left[\begin{array}{} 
\underline{x} \\ \underline{t}
\end{array}\right]\leq \left[\begin{array}{} 
\underline{b} \\ -\underline{b}
\end{array}\right]$$

Thus we add 38 slack variables using $I_{38}$, naming the variables $s_1,\,\dots,\,s_{38}$ which form the vector $\underline{s}$. This results in the standard form LP constraint.

$$\left[\begin{array}{}\begin{array}{} 
A \\ -A \end{array} & \bigg| & \begin{array}{} -I_{19} \\ -I_{19} 
\end{array} & \bigg| & I_{38} \end{array}\right]\left[\begin{array}{} 
\underline{x} \\ \underline{t} \\ \underline{s}
\end{array}\right] \left[\begin{array}{} 
\underline{b} \\ -\underline{b}
\end{array}\right]$$


The LP in standard form is:


$$\begin{align*}
\operatornamewithlimits{\textsf{minimise}}_{\underline{x}}  \quad & \underline{1}^T\underline{t} \\
\textsf{subject to} \quad & \left[\begin{array}{}\left.\begin{array}{} 
A \\ -A \end{array}\right| & \left.\begin{array}{} -I_{19} \\ -I_{19} 
\end{array} \right| & I_{38} \end{array}\right]\left[\begin{array}{} 
\underline{x} \\ \underline{t} \\ \underline{s}
\end{array}\right] = \left[\begin{array}{} 
\underline{b} \\ -\underline{b}
\end{array}\right]\\
& \underline{t},\,\underline{s} \geq \underline{0}
\end{align*}$$


-----------------


**b)** We let:


$$X=\left[\begin{array}{}\begin{array}{} 
A \\ -A \end{array} & \bigg| & \begin{array}{} -I_{19} \\ -I_{19} 
\end{array} \end{array}\right]$$


$$\underline{c} = \left[\begin{array}{} 
\underline{b} \\ -\underline{b}\quad
\end{array}\right]$$


We define our objective function through the vector $\underline{f}$

$$\underline{f} = \left[\begin{array}{} 
0 \\ 0 \\ 1 \\ \vdots \\ 1
\end{array}\right]$$


Such that $\underline{f}$ is $21\times1$. 


As $\underline{1}^T \underline{t} = \underline{0}^T \underline{x} + \underline{1}^T\underline{t}$


The function $\textsf{linprog}(f,\,X,\,c)$ will minimise $f$, such that $X\underline{x}^* \leq \underline{c}$, where $\underline{x}^*$ is the vector of all our constraints. 

This is equivalent to:


$$\left[\begin{array}{}\begin{array}{} 
A \\ -A \end{array} & \bigg| & \begin{array}{} -I_{19} \\ -I_{19} 
\end{array} \end{array}\right]\left[\begin{array}{} 
\underline{x} \\ \underline{t}
\end{array}\right]\leq \left[\begin{array}{} 
\underline{b} \\ -\underline{b}
\end{array}\right]$$


Hence, using $\textsf{linprog}(f,\,X,\,c)$ will solve our L.P.

The solution to this LP using the $\textsf{linprog}$ function. which satisfies the non-negativity constraints of $\underline{t}$ and $\underline{s}$ is:

$$\left[\begin{array}{} 
-0.0083 \\ 26.4085 \\ 0 \\ 0.1331 \\ 0.2662 \\ 0.0008 \\ 0.1677 \\ 0.1846 \\ 0.0385 \\ 0.0085 \\ 0.2146 \\ 0.0123 \\ 0.0492 \\ 0.0238 \\ 0.0631 \\ 0 \\ 0.0131 \\ 0.1138 \\ 0.1408 \\ 0.0723 \\ 0.0954
\end{array}\right]$$


-------


**c)** Given our solution to the LP we find:

$$x_1=-0.0083\quad x_2=26.4085$$

This gives  a linear approximation of 

$$T=26.4085 - 0.0083x$$


![description](https://raw.githubusercontent.com/Reezy-git/303-Linear-Algebra/3671b084a4555314b92d25f4618ce1aad2987d39/plot_2.svg "Plot")

```Figure 1. - A linear prediction of winning 100m olympic sprint times.```

<br>

----------

**d)** Using the linear approximation for the line of best fit

$$T=26.4085 - 0.0083x$$

and inputting the year $x=2024$ yields

$$T=26.4085 - 0.0083\times2024=9.6093$$

Our linear approximation predicts a winning 100m spring time of $9.6093$ seconds at the 2024 Olympic Games.
