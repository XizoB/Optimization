### <font size=5 face=Times New Roman> *A Semismooth Newton Method for $\mathrm{L}^1$ Data Fitting with Automatic Choice of Regularization Parameters and Noise Calibration*</font>

$Analyze$  $By$  $Xizhou Bu$

#### 1.Abstract

文章考虑具有$\mathrm{L}^1$数据拟合项问题（缺乏可微性）的数值解。利用凸对偶性，将该问题重新表述为<font color=red>具有点态约束</font>的光滑函数最小化问题，该问题可以用半光滑牛顿法（Semismooth Newton Method）有效地求解。（把原问题转换成凸对偶问题，用半光滑牛顿法求解）

为了实现超线性收敛，对偶问题需要额外的正则化，并且对原问题和对偶问题的正则化参数的选择都是关键。文章中提出了选择这些参数的自适应策略。其中，原问题公式中的正则化参数是根据<font color=red>模型函数法的平衡原则</font>选择的，而对偶问题公式中的正则化参数是根据<font color=red>最优条件结构的路径跟踪策略</font>确定的。文章通过数值实验验证了该方法和自适应策略的有效性和鲁棒性。

#### 2.Introduction

问题描述：
$$
Kx=y^\delta
$$
这个问题在Hadamard的意义上是不适定的（ill-posed，就是说K的条件数大，即为病态的），特别是，解决方案常常不能持续地依赖于数据。现在的标准方法是<font color=red>Tikhonov正则化</font>：
$$
F(\alpha)=||Kx_\alpha-y^\delta||_\mathrm{L^2}^2+\alpha R(x)
$$
其中正则化项$R(x)$是根据应用选择的，本文$R(x)$选择Euclid范数（L2范数），适用于光滑的解决方案，可以得到一个（凸）二次优化问题：

<font color=red> $Primal \ Problem:$ </font>
$$
\underset {x \in L^2}{min} \ \lbrace \ {\cal J}_\alpha(x)\equiv||Kx_\alpha-y^\delta||_\mathrm{L^1}+\frac{\alpha}{2}||x_\alpha||_\mathrm{L^2}^2 \rbrace
$$
关于值函数中$\alpha$参数的属性，我们转换为以下，
$$
F(\alpha)=||Kx_\alpha-y^\delta||_\mathrm{L^1}^2+\frac{\alpha}{2}||x_\alpha||_\mathrm{L^2}^2
$$

由下式可知，$F(\alpha)$连续单调递增，其中由定理2.2可得，$||Kx_\alpha-y^\delta||_\mathrm{L^1}^2$关于$\alpha$参数连续单调递增，而$||x_\alpha||_\mathrm{L^2}^2$连续单调递减，
$$
F^{'}(\alpha)=\frac{1}{2}||x_\alpha||_\mathrm{L^2}^2
$$

**Theorem 2.2:**

函数$||Kx_{\alpha}-y^{\delta}||$和$||x_{\alpha}||^2_{\mathrm{L^2}}$是连续的，并且对于$\alpha$来说分别是单调递增和单调递减的
$$
\begin{aligned}
&(\alpha_1-\alpha_2)(||Kx_{\alpha_1}-y^\delta||_\mathrm{L^1}^2-||Kx_{\alpha_2}-y^\delta||_\mathrm{L^1}^2)\geq0
\\
&(\alpha_1-\alpha_2)(||x_{\alpha_1}||_\mathrm{L^2}^2-||x_{\alpha_2}||_\mathrm{L^2}^2)\leq0
\end{aligned}
$$


<font color=red> $Dual \ Problem:$ </font>
$$
\begin{cases}
        \underset {p\in \mathrm{L^2}}{min}\quad \frac{1}{2\alpha}||K^*p||^2_{\mathrm{L^2}}-\langle p,y^\delta \rangle_{\mathrm{L^2}}  \\
        s.t. \quad||p||_{\mathrm{L^\infin}}\leq1
        \end{cases}
$$


**Theorem 2.5:**

​		对偶问题至少有一个解，其中原问题解$x_\alpha\in\mathrm{L^2}$与对偶问题解$p_\alpha\in\mathrm{L^2}$的关系如下：
$$
\begin{cases}
\begin{aligned}
K^*p_\alpha &=\alpha x_\alpha
\\
0 &\leq \langle Kx_\alpha-y^\delta,p-p_\alpha \rangle_\mathrm{L^2}
\end{aligned}
\end{cases}
$$

对于$\forall p\in \mathrm{L^2}$有$||p||_\mathrm{L^{\infin}}\leq 1$。



#### 3.Sulution by Semismooth Newton method（半光滑牛顿法解决对偶问题）

##### 3.1 Regularization

如果$K$的转置是不适定的（病态），尽管在$p$上有点态边界，那么对偶问题（$P^*$）也是不适定的。为了对抗这种情况并确保求解约束优化问题的半光滑牛顿法超线性收敛，如下，文章引入了正则化问题
$$
\begin{cases}
        \underset {p\in \mathrm{L^2}}{min}\quad \frac{1}{2\alpha}||K^*p||^2_{\mathrm{L^2}}- \frac{\beta}{2}||\nabla p||_\mathrm{L^2}^{2} - \langle p,y^\delta \rangle_{\mathrm{L^2}}  
        \\
        s.t. \quad||p||_{\mathrm{L^\infin}}\leq1
\end{cases}
$$
对于$\beta>0$，$p$的点态界和半范数正则化之间的相互作用将使正则化参数$\beta$易于选择（文章中4.2节有详细解释），我们假设$ker\ K^*\cap ker\  \nabla=\{0\}$，也就是说常数函数不属于$K^*$的内核。在这个假设下内积$\frac{1}{\alpha} \langle K^*.,K^*. \rangle+\beta \langle \nabla.,\nabla.\rangle$在$\mathrm{L^1}$上推导出一个等价范数，并且问题（$P^*_\beta$）有一个唯一的解。如果半范数正则化被完整的$\mathrm{H^1}$规范所取代，这个假设可以被消除。

为了在数值上解决对偶问题（$P^*_\beta$），我们引入框约束的$Moreau-Yosida$正则化：
$$
\begin{aligned}
\underset {p\in \mathrm{H^1}}{min}\quad\frac{1}{2\alpha}||K^*p||^2_{\mathrm{L^2}}- &\frac{\beta}{2}||\nabla p||_\mathrm{L^2}^{2}- \langle p,y^\delta \rangle_{\mathrm{L^2}}
\\
&+\frac{1}{2c}||max(0,c(p-1))||^2_\mathrm{L^2}+\frac{1}{2c}||min(0,c(p+1))||^2_\mathrm{L^2}
\end{aligned}
$$
对于$c>0$，这里的最大值和最小值是逐点取得。对于固定的$\beta和c$，在上述假设下$\frac{1}{2\alpha}||K^*p||^2_{\mathrm{L^2}}- \frac{\beta}{2}||\nabla p||_\mathrm{L^2}^{2}$是严格凸的，因此对偶问题（$P^*_{\beta,c}$）有一个最小化解$p_c$，问题（$P^*_{\beta,c}$）的最优系统如下：
$$
\begin{cases}
\frac{1}{\alpha}KK^*p_c-\beta\Delta p_c-y^{\delta}+\lambda_c=0\\
\lambda_c=max(0,c(p_c-1))+min(0,c(p_c+1))

\label{eq3.1}
\end{cases}
$$


##### 3.2 Semismooth Newton method

正则化最优性系统可以使用超线性收敛的半光滑牛顿法进行有效求解，为此我们认为是一个非线性方程$F(p)=0$，且$F:\mathrm{H^1}\rightarrow(\mathrm{H^1})^*$，
$$
F(p):=\frac{1}{\alpha}KK^*p-\beta\Delta p+max(0,c(p-1))+min(0,c(p+1))-y^{\delta}
$$
已知投影算子
$$
P(p):=max(0,c(p-1))+min(0,c(p+1))
$$
它从$L_q$到$L_p$是半光滑的，当且仅当$q>p$的情况下，它是牛顿导数
$$
D_NP(p)h=h{\chi_{\{|p|>1\}}}:=\begin{cases}
h(x)&if|p(x)|>1
\\
0&if|p(x)|\leq1
\end{cases}
$$
由于$Frechet-differentiable$函数与半光滑函数的和是半光滑的（具有标准牛顿导数），所以$F$是半光滑的，它的牛顿导数是
$$
D_NP(p)h=\frac{1}{\alpha}KK^*h-\beta\Delta h+ch{\chi_{\{|p|\}}}
$$
半光滑牛顿法的迭代步骤在于求解$p^{k+1}\in\mathrm{H^1}$方程
$$
D_NF(p^k)(p^{k+1}-p^{k})=-F(p^k)
$$

##### 3.3 Algorithm 1 Semismooth Newton method for公式$\eqref{eq3.1}$

$$
\frac{1}{\alpha}KK^*p^{k+1}-\beta\Delta p^{k+1}+c\chi_{A_k}p^{k+1}=y^{\delta}+c(\chi_{A^+_k}-\chi^-_{A_k})
$$

<img src="C:\Users\HP\AppData\Roaming\Typora\typora-user-images\image-20211026120303291.png" alt="image-20211026120303291" style="zoom:50%;" />

求解$p^{k+1}$

#### 4.自适应正则化参数寻优策略

##### 4.1平衡原理

正则化参数选取的思想是平衡数据拟合项$\varphi(\alpha)= ||Kx_\alpha-y^\delta||_\mathrm{L^1}$和惩罚项$\alpha F^{'}(\alpha)=\frac{\alpha}{2}||x_\alpha||_\mathrm{L^2}^2$，因此平衡原理包括选择正则化参数$\alpha^*$作为方程的解：
$$
(\sigma-1)\varphi(\alpha^*)=\alpha^*F^{'}(\alpha^*)
\label{4.1}
$$
其中，$\sigma>1$是控制这两项之间的相对权重。类似的平衡思想存在于许多启发式参数选择规则之下。

**Theorem 4.1:**

对于$\sigma$足够接近与1并且$y^{\delta}\neq0$，平衡方程$\eqref{4.1}$至少存在一个正解

**Lemma 4.2:**

以下极限是成立的，
$$
\underset {\alpha\rightarrow 0^+}{lim}\  \frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}=\underset {\alpha\rightarrow +\infin}{lim}\  \frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}=0
$$
接下来，我们将反复利用到残差的概念，由定理2.2可知，函数$\varphi(\alpha)$和$F^{'}(\alpha)$是连续的，因此残差函数$r(\alpha)$也是连续的
$$
r(\alpha)=\alpha F^{'}(\alpha)-(\sigma-1)\varphi(\alpha)
\label{4.2}
$$
由引理4.2可知，$\sigma>1$时，有以下极限：
$$
\begin{aligned}
&\underset {\alpha\rightarrow 0^+}{lim}\ r(\alpha)=-(\sigma-1)\underset {\alpha\rightarrow 0^+}{lim}\ ||Kx_{\alpha}-y^\delta||_\mathrm{L^1}\leq0
\\
&\underset {\alpha\rightarrow +\infin}{lim}\ r(\alpha)=-(\sigma-1)\underset {\alpha\rightarrow 0^+}{lim}\ ||y^{\delta}||_\mathrm{L^1}<0
\end{aligned}
$$
且有$||Kx_{\alpha}-y^\delta||_\mathrm{L^1}\leq ||y^{\delta}||_\mathrm{L^1}$，以及$sup_{\alpha\in(0,+\infin)}\frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}>0$，故可以得到：
$$
r(\alpha)=\frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}-(\sigma-1)||Kx_{\alpha}-y^\delta||_\mathrm{L^1}\geq \frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}-(\sigma-1)||y^{\delta}||_{\mathrm{L^1}}
$$
所以，对于所有的$\sigma\in(1,\sigma_0)$，存在一个$\sigma_0>1使得$$sup_{\alpha\in(0,+\infin)}\frac{\alpha}{2}||x_{\alpha}||^2_{\mathrm{L^2}}>0$存在，即残差函数$r(\alpha)$存在一个正解，$r(\alpha)$函数的基本形式如下图所示：

<img src="D:\OneDriveX\OneDrive\桌面\Course_PostGraduate\凸优化\组号19：卜习州 张伟杰 熊婧伊\附件\r(a).jpg" alt="r(a)" style="zoom:25%;" />



##### 4.2模型方程与平衡点迭代

为了求解平衡方程的解，我们把公式$\eqref{4.1}$改写成
$$
F(\alpha^*)=\sigma(F(\alpha^*)-\alpha^*F^{'}(\alpha^*))
$$
其中转换步骤为把$F(\alpha^*)=\varphi(\alpha)+\alpha^*F^{'}(\alpha^*)$带入可得，
$$
\begin{aligned}
\varphi(\alpha)+\alpha^*F^{'}(\alpha^*)&=\sigma(\varphi(\alpha)+\alpha^*F^{'}(\alpha^*)-\alpha^*F^{'}(\alpha^*))
\\
(\sigma-1)\varphi(\alpha^*)&=\alpha^*F^{'}(\alpha^*)
\end{aligned}
$$
如下图所示，其中当$\sigma\rightarrow1$的时候，函数$F(\alpha)$在$\alpha^*$处的切线与纵轴之间的截距$(F(\alpha^*)-\alpha^*F^{'}(\alpha^*))$与其函数值$F(\alpha^*)$相等，这也就意味着$F(\alpha^*)$斜率为零，即$\alpha^*$为所需要找的平衡点。

<img src="D:\OneDriveX\OneDrive\桌面\Course_PostGraduate\凸优化\组号19：卜习州 张伟杰 熊婧伊\附件\F(a).jpg" alt="F(a)" style="zoom:9%;" />

考虑平衡点迭代，其中选取$\alpha^{k+1}$作为解
$$
F(\alpha^{k+1})=\sigma(F(\alpha^{k})-\alpha^kF^{'}(\alpha^{k}))
\label{4.3}
$$
为了计算这个解，我们使用[20]中提出的确定正则化参数的模型函数方法，该方法通过有理多项式局部逼近值函数$F(\alpha)$，在本文中考虑了这样的一个模型函数的形式
$$
m(\alpha)=b+\frac{s}{t+\alpha}
$$
注意对于$\alpha\rightarrow \infin$的时候，$x_{\alpha}\rightarrow0$，我们固定$b=||y^{\delta}||_{\mathrm{L^1}}$来使$m(\alpha)$匹配$F(\alpha)$的渐近特性（尽管$b$值较大，也适合我们的目的），而参数$s和t$由插值条件来决定
$$
m(\alpha)=F(\alpha),\quad m^{'}(\alpha)=F^{'}(\alpha)
$$
由$m(\alpha)$的定义可以得到
$$
b+\frac{s}{t+\alpha}=F(\alpha),\quad -\frac{s}{(t+\alpha)^2}=F^{'}(\alpha)
$$
参数$s和t$可以用显格式推导出，回想一下，根据定理2.2，我们有$F^{'}=\frac{1}{2}||x_{\alpha}||^2_{\mathrm{L^2}}$，并且该值不需要占用额外的计算工作。如果我们用模型函数$m_k(\alpha)$的值$m_k(\alpha^{k+1})$替换公式$\eqref{4.3}$的左侧，其中模型函数由$\alpha^k$处的插值条件推导出的。因此我们可以显式地计算一个新迭代。结果迭代由算法2给出。
$$
\begin{aligned}
\frac{s_k}{t_k+\alpha_k}&=F(\alpha_k)-b,\quad -\frac{s_k}{(t_k+\alpha_k)^2}=F^{'}(\alpha_k)
\\
s_k&=-\frac{(b-F(\alpha_k))^2}{F^{'}(a_k)}
\\
t_k&=\frac{b-F(\alpha_k)}{F^{'}(a_k)}-\alpha_k
\end{aligned}
$$
其显著的特征在于不需要了解噪声水平，事实上，由于$F(0)$表示噪声水平的下限，如果模型函数$m(\alpha)$在$\alpha=0$的邻域内合理地逼近价值函数$F(\alpha)$，则$m(0)$量可以被认为是噪声水平的有效估计。

找到算法二的平衡点$\bar{\alpha}$，我们发现---使用事实$\hat{m}=m(\hat{\alpha})$满足$\bar{\alpha}$处的插值条件---即$\bar{\alpha}$是平衡方程$\eqref{4.1}$的解。我们可以证明，在一些普遍的假设下，算法二局部收敛到这样的一个固定的点。为此，我们称公式$\eqref{4.1}$的解为正则吸引算子，如果存在一个$\epsilon>0$，使得对于$\alpha\in(\alpha^*-\epsilon,\alpha^*)使得r(\alpha)<0，以及对于\alpha\in(\alpha^*,\alpha^*+\epsilon)使得 r(\alpha)<0$。

**Theorem 4.3**

假设$\alpha_0$满足$r(\alpha_0)<0$并且足够接近正则吸引子$\alpha ^*$，由算法2生成的$\{\alpha_k\}$序列收敛于$\alpha ^*$

<img src="C:\Users\HP\AppData\Roaming\Typora\typora-user-images\image-20211026113712199.png" alt="image-20211026113712199" style="zoom:50%;" />

##### 4.3平衡点迭代的收敛性

为了分析平衡点迭代算法的收敛性，我们观察到
$$
\begin{aligned}
\alpha_{k+1}&=\frac{s_k}{\sigma \hat{m}-b}-t_k
\\
&=\frac{(F(\alpha_k)-b)^2-(b-F(\alpha_k)-\alpha_kF^{'}(\alpha_k))(b-\sigma\varphi(\alpha_k))}{F^{'}(\alpha_k)(b-\sigma\varphi(\alpha_k))}
\end{aligned}
$$
其中$\varphi(\alpha)=||Kx_\alpha-y^{\delta}||$表示残差的范数，分子可化简为：
$$
\begin{aligned}
(F(\alpha_k)-b)^2-(b-F(\alpha_k)-&\alpha_kF^{'}(\alpha_k))(b-\sigma\varphi(\alpha_k))
\\
&=(\alpha_kF^{'}(\alpha_k))^2+(\sigma-1)\varphi(\alpha_k)[b-F(\alpha_k)-\alpha_kF^{'}(\alpha_k)]
\end{aligned}
$$
因此平衡点迭代为
$$
\alpha_{k+1}=\frac{(\alpha_kF^{'}(\alpha_k))^2+(\sigma-1)\varphi(\alpha_k)[b-F(\alpha_k)-\alpha_kF^{'}(\alpha_k)]}{F^{'}(\alpha_k)(b-\sigma\varphi(\alpha_k))}=:\alpha_k\frac{N_k}{D_k}

\label{4.4}
$$
在$b>\sigma||y^{\delta}||_{\mathrm{L^1}}$的假设下，分母$D_k$是正的，这样迭代就好定义了：
$$
N_k-D_k=[(\sigma-1)\varphi(\alpha_k)-\alpha_kF^{'}(\alpha_k)][b-F(\alpha_k)]
$$
因此从公式$\eqref{4.4}$可得，如果$r(\alpha_k)=\alpha_kF^{'}(\alpha_k)-(\sigma-1)\varphi(\alpha_k)>0$，$N_k<D_k$，则$\alpha_{k+1}<\alpha_{k}$；否则$\alpha_{k+1}>\alpha_{k}$。公式$\eqref{4.4}$可以转换为：
$$
\begin{aligned}
\alpha_{k+1}-\alpha_{k}&=\alpha_{k}\frac{N_k}{D_k}-\alpha_{k}=\frac{N_k-D_k}{F^{'}(\alpha_k)(b-\sigma\varphi(\alpha_k))}
\\
&=\frac{[(\sigma-1)\varphi(\alpha_k)-\alpha_kF^{'}(\alpha_k)][b-F(\alpha_k)]}{F^{'}(\alpha_k)(b-\sigma\varphi(\alpha_k))}
\\
&=[(\sigma-1)\frac{\varphi(\alpha_k)}{F^{'}(\alpha_k)}-\alpha_k] \frac{b-F(\alpha_k)}{b-\sigma\varphi(\alpha_k)}
\\
&=[T(\alpha_k)-\alpha_k]\frac{b-F(\alpha_k)}{b-\sigma\varphi(\alpha_k)}
\end{aligned}
$$


定义辅助算子$T(\alpha)=(\sigma-1)\frac{\varphi(\alpha)}{F^{'}(\alpha)}，\omega_k:=\frac{b-F(\alpha_k)}{b-\sigma\varphi(\alpha_k)}$，其中辅助算子$T$可以被视为当标量$b\rightarrow+\infin$的时候，算子$\alpha_k\frac{N_k}{D_k}$的渐近解，而对于$b>\sigma||y^{\delta}||_{\mathrm{L^1}}$，$\omega_k>0$

公式$\eqref{4.4}$可以转换为
$$
\alpha_{k+1}=\omega_kT(\alpha_k)+(1-\omega_k)\alpha_k
\label{34}
$$
**Lemma 4.4**

如果$0<\alpha_0<\alpha_1$，算子$T$是单调的，有$T(\alpha_0)<T(\alpha_1)$

**Lemma 4.5**

对于任意初始设想的$\alpha_0$，$\{T^k(\alpha_0)\}$序列是单调的，此外当$r(\alpha_0)>0$的时候，它是单调递减的，当$r(\alpha_0)<0$的时候，它是单调递增的。

**Lemma 4.8**

对于任意$\alpha_0$，使$\{T^k(\alpha_0)\}$收敛于$\alpha^*$，函数$r(\alpha)$不会在开区间$(min(\alpha_0,\alpha^*),max(\alpha_0,\alpha^*))$上消失

**Theorem 4.9**

假设$\alpha_0$满足$r(\alpha_0)>0$，由平衡点迭代公式$\eqref{4.4}$生成的$\{\alpha_k\}$序列单调递减并且收敛于公式$\eqref{4.1}$

**Theorem 4.10**

假设初始猜想$\alpha_0$满足$r(\alpha_0)<0$，由平衡点迭代公式$\eqref{4.4}$生成的$\{\alpha_k\}$序列单调递增或者存在某些$k_0\in N$使得所有$k\geq k_0时，r(\alpha_k\geq 0)$

如下图所示迭代机理：

<img src="D:\OneDriveX\OneDrive\桌面\Course_PostGraduate\凸优化\组号19：卜习州 张伟杰 熊婧伊\附件\a.jpg" alt="a" style="zoom:30%;" />

第一种情况当$r(\alpha_0)>0$的时候，由定理4.9可以得出，$0<\omega_0<1$，迭代公式$\eqref{34}$生成的序列$\{\alpha_k \}$单调递减，由残差函数$r(\alpha)$的基本形状可知，$\alpha^*<\alpha_0$，由引理4.8可知$\alpha^*<\{T^k(\alpha_0)\}<\alpha_0$，此时辅助算子$T$相当于“方向引导”把$\alpha_0$“向左拉”到$\alpha_1$去靠近$\alpha^*$，而$\omega$相当于“拉的程度”，故有$\alpha^*<\{T^k(\alpha_0)\}<\alpha_1<\alpha_0$

第二种情况当$r(\alpha_0)<0$的时候，由定理4.10可以得出，$\omega_0>1$，迭代公式$\eqref{34}$生成的序列$\{\alpha_k \}$单调递增，由残差函数$r(\alpha)$的基本形状可知，$\alpha^*>\alpha_0$，由引理4.8可知$\alpha_0<\{T^k(\alpha_0)\}<\alpha^*$，此时辅助算子相当于“方向引导”把$\alpha_0$“向右拉”到$\alpha_1$去靠近$\alpha^*$。这时$\alpha_1$与$\alpha^*$的关系可能有两种，一种是$\alpha_1>\alpha^*，这时利用定理4.9，如第一种情况最终收敛于\alpha^*，第二种是$$\alpha_1<\alpha^*$，这样继续迭代直至$\alpha_1>\alpha^*$，再利用定理4.9，直至收敛于$\alpha^*$。



##### 4.3路径跟踪法选取参数$\beta$

<img src="C:\Users\HP\AppData\Roaming\Typora\typora-user-images\image-20211026120342878.png" alt="image-20211026120342878" style="zoom:50%;" />

