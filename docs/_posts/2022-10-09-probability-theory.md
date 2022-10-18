---
layout: post
title: 《概率论与数理统计》知识简单梳理
tags: [machine learning, probability theory]
---

- 随机试验与随机事件

试验，  
随机试验，相同条件下可重复、结果不止一个、无法预测  
事件，试验结果  
基本事件，  
复合事件，基本事件复合构成  
样本空间，所有基本事件的集合，$\Omega$  
样本点，样本空间的元素  
必然事件，  
随机事件，

- 事件间的关系

1、并集
$$A+B\supset A$$
$$A+A=A$$
$$A+\Phi=A$$
$$A+\Omega=\Omega$$
2、交集
$$AB\subset A$$
$$AA=A$$
$$A\Phi=\Phi$$
$$A\Omega=A$$
3、差集
$$A-B=A-AB$$
4、互斥事件，$AB=\Phi$  
5、对立事件，$AB=\Phi\ \&\ A+B=\Omega$
$$B=\overline A$$
6、完备事件组，$A_1,A_2,\dots,A_n$两两互不相容，且$\cup_{i=1}^nA_i=\Omega$

运算律：

1、结合
$$(A\cup B)\cup C=A\cup (B\cup C)$$
$$(A\cap B)\cap C=A\cap (B\cap C)$$
2、分配
$$(A\cup B)\cap C=(A\cap C)\cup(B\cap C)$$
$$(A\cap B)\cup C=(A\cup C)\cap(B\cup C)$$
3、对偶
$$\overline{A\cup B}=\overline A\cap\overline B$$
$$\overline{A\cap B}=\overline A\cup\overline B$$

- 古典概率模型

条件：1）有限个样本点；2）等可能性
$$P(A)=\frac{A中包含的基本事件数}{基本事件总数}=\frac{A的有利样本点}{\Omega的样本点总数}$$

排列，加法原理、乘法原理  
不重复排列，$P_n^m=n(n-1)(n-2)\dots(n-m+1)=\frac{n!}{(n-m)!}$  
全排列，$P_n^n=n!,\ 0!=1$

组合，  
$C_n^m=\frac{P_n^m}{m!}=\frac{n!}{m!(n-m)!}$  
$C_n^m=C_n^{n-m}$  

- 几何概率模型

核心在将问题转为几何描述。
$$P(A)=\frac{\mu(G)}{\mu(\Omega)}$$

- 公理

1）非负性，$0\le P(A)\le1$  
2）规范性，$P(\Omega)=1$  
3）完全可加，$P(A_1+A_2+\dots)=P(A_1)+P(A_2)+\dots$，其中$A_1,A_2,\dots$互不相容

- 基本公式

1、加法公式
$$P(A+B)=P(A)+P(B)-P(AB)$$
$$P(A+B+C)=P(A)+P(B)+P(C)-P(AB)-P(AC)-P(BC)+P(ABC)$$
2、条件概率
$$P(A|B)=\frac{P(AB)}{P(B)}$$
$$P(\sum_{i=1}^\infty A_i|B)=\sum_{i=1}^\infty P(A_i|B)$$
3、乘法公式
$$P(AB)=P(B)P(A|B)=P(A)P(B|A),\ P(A)>0\ \&\ P(B)>0$$
$$P(A_1A_2\dots A_n)=P(A_1)P(A_2|A_1)P(A_3|A_1A_2)\dots P(A_n|A_1A_2\dots A_{n-1})$$
4、全概率公式
$$P(B)=\sum_{i=1}^nP(A_i)P(B|A_i)$$
其中$A_1,A_2,\dots,A_n$是E的完备事件组。  
5、贝叶斯公式
$$P(A_k|B)=\frac{P(A_k)P(B|A_k)}{\sum_{i=1}^nP(A_i)P(B|A_i)}=\frac{P(A_k)P(B)}{P(B)}$$
其中$A_1,A_2,\dots,A_n$是E的完备事件组，B是任意随机事件。$P(A_i)$称为先验概率，$P(A_i|B)$称为后验概率。  
6、独立性
$$P(AB)=P(A)P(B),\ P(A)>0,P(B)>0$$
注意**独立**$P(AB)>0$和**互不相容**$P(AB)=0$的区别

- 伯努利模型
A的概率是p，重复n次，A发生k次，则
$$P_n(k)=C_n^kp_kq^{n-k},\ q=1-p$$

- 离散型随机变量及其概率分布

对于随机变量X的所有可列$x_k(k=1,2,\dots)$事件,
$$P\{X=x_k\}=p_k$$
使用概率分布表表示。

- 连续性随机变量及其概率密度

$$P\{a<X\le b\}=\int_a^bf(x)dx$$
X是连续随机变量，$f(x)$概率密度函数。

- 分布函数

离散变量和连续变量均适用，离散变量是分段函数；连续变量是连续函数，且为概率密度的积分

$$F(x)=P(X\le x)$$
X取值不超过x的概率，$x\in(-\infty,\infty),\ F(x)\in[0,1]$，性质：
$$\lim_{x\rightarrow+\infty}F(x)=F(+\infty)=1$$
$$\lim_{x\rightarrow-\infty}F(x)=F(-\infty)=1$$

- 常见随机变量的分布

1、0-1分布
$$P\{X=k\}=p^k(1-p)^{1-k},\ k=\{0,1\}$$
$EX=p,\ DX=pq$

2、几何分布  
第k次发生，
$$P\{X=k\}=(1-p)^{k-1}p,\ X\sim G(p)$$
$EX=1/p,\ DX=(1-p)/p^2$

3、二项分布  
发生k次，
$$P\{X=k\}=C_n^kp^k(1-p)^{n-k},\ k=0,1,2,\dots,n.\ X\sim B(n,p)$$
$EX=np,\ DX=npq$

4、泊松分布
$$P\{X=k\}=\frac{\lambda^k}{k!}e^{-\lambda},\ k=0,1,2,3,\dots,\  \lambda>0,\ X\sim P(\lambda)$$
$EX=\lambda,\ DX=\lambda$

5、超几何分布
$$P\{X=k\}=\frac{C_{N_1}^kC_{N_2}^{n-k}}{C_N^n},\ k=0,1,2,\dots,min(n,N_1)$$
可用于描述不放回试验

6、均匀分布  
$$f(x)=\begin{cases} \frac1{b-a} &a\le x\le b,\\0&else\end{cases},\ X\sim\cup[a,b]$$
$EX=(a+b)/2,\ DX=(b-a)^2/12$

7、指数分布
$$f(x)=\begin{cases} \lambda e^{-\lambda x} &x>0,\\0&x\le0\end{cases},\ \lambda>0,\ X\sim Exp(\lambda)$$
$EX=1/\lambda,\ DX=1/\lambda^2$

8、正态分布
$$\phi(x)=\frac1{\sqrt{2\pi}\sigma}e^{-\frac{(x-\mu)^2}{2\sigma^2}},-\infty<x<\infty,X\sim N(\mu,\sigma^2)$$
$EX=\mu,\ DX=\sigma^2$  
标准正态分布，
$$\phi_0(x)=\frac1{\sqrt{2\pi}}e^{-\frac{x^2}2},-\infty<x<\infty,X\sim N(0,1)$$

- 随机变量函数的分布

X的密度函数为$f_X(x)$，$Y=kx+b(k\ne0)$，则$f_Y(x)=\frac1{|k|}f_X(\frac{x-b}k)$
$$F_Y(x)\rightarrow F_X(x)$$
$$f_Y(x)\leftarrow f_X(x)$$
利用以上转换求得

- 多维随机变量及其分布

二维随机变量的分布函数  
联合分布，
$$F(x,y)=P\{X\le x,Y\le y\}$$
边缘分布，
$$F_X(x)=P\{X\le x,Y<+\infty\}$$
$$F_Y(y)=P\{X<+\infty,Y\le y\}$$

二维正态分布的两个边缘分布也是正正态分布，反之不一定。

- 条件分布

$$F(x|A)=P\{X\le x|A\}$$

- 数学期望

1、离散型
$$P(X=x_k)=P_k,\ 则EX=\sum_{k=1}^\infty x_kP_k$$
满足绝对收敛

2、连续型
连续变量X的概率密度函数$f(x)$，则
$$EX=\int_{-\infty}^{+\infty}xf(x)dx$$

3、性质
$$EC=C$$
$$E(X+C)=EX+C$$
$$E(CX)=CEX$$
$$E(kX+C)=kEX+C$$
$$E(X\pm Y)=EX\pm EY$$
$$E(XY)=EX\cdot EY,X和Y相互独立$$

- 方差

$$DX=E(X-EX)^2$$
可变换为$DX=EX^2-(EX)^2$，$\sqrt {DX}$为标准差

- 协方差

$$Cov(X,Y)=[(X-EX)(Y-EY)]=E(XY)-EXEY$$
如果X,Y相互独立，则$Cov(X,Y)=0$。

- 相关系数

$$\rho=\frac{Cov(X,Y)}{\sqrt{DX}\sqrt{DY}}$$
$|\rho|\le1$，$\rho$反应的是线性关系。  
1）$\rho=1$，X和Y完全正相关；  
2）$\rho=-1$，X和Y完全负相关；  
3）$|\rho|$接近0，X和Y线性关系趋弱；  
4）$\rho=0$，X和Y不存在线性关系

- 大数定律

大量重复试验的平均结果的稳定性。  
切比雪夫不等式：  
对于X，若$EX和DX$存在，则$\forall\epsilon>0$
$$P(|X-EX|\ge\epsilon)\le\frac{DX}{\epsilon^2}$$
切比雪夫大数定律：  
$X_1,X_2,\dots,X_n$两两不相关，$EX_i和DX_i$都存在，且$DX_i\le M$，则$\forall\epsilon>0$，有
$$\lim_{n\rightarrow\infty}P\left\{|\frac1n\sum_{i=1}^nX_i-\frac1n\sum_{i=1}^nEX_i|<\epsilon\right\}=1$$
辛钦大数定律：  
$X_1,X_2,\dots,X_n$独立同分布，$EX_i=\mu$，则$\forall\epsilon>0$，有
$$\lim_{n\rightarrow\infty}P\left\{|\frac1n\sum_{i=1}^nX_i-\mu|<\epsilon\right\}=1$$

- 中心极限定理

若一个随机现象由大量相互独立的因素影响，大量独立同分布的变量和的极限分布是正态分布。
$X_1,X_2,\dots,X_n$独立同分布（不论是和何种分布），$EX_i=\mu,\ DX_i=\sigma^2$，则
$$\lim_{n\rightarrow\infty}P\left(\frac{\sum_{i-1}^nX_i-n\mu}{\sqrt n\sigma}\le x\right)=\Phi_0(x),\ 0<\sigma^2<+\infty$$

- 常用统计量

样本$(X_1,X_2,\dots,X_n)$来自总体$X$，$EX=\mu,\ DX=\sigma^2$，则有：  
样本均值：
$$\overline X=\frac1n\sum_{i=1}^nX_i$$
未修正的样本方差：
$$S_0^2=\frac1n\sum_{i=1}^n(X_i-\overline X)^2$$
样本方差：
$$S^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X)^2$$
样本均值和样本方差的性质：
$$E\overline X=\mu$$
$$D\overline X=\frac1n\sigma^2$$
$$ES^2=\sigma^2$$

- 抽样分布

即统计量的分布。

1、卡方分布  
若$X_1,X_2,\dots,X_n$相互独立，且$X_i\sim N(0,1)$，则
$$\chi^2=\sum_{i=1}^nX_i^2\sim\chi^2(n)$$
$\chi^2$的$EX=n,\ DX=2n$；由中心极限定理，n充分大时，卡方分布近似为标准正态分布。  
若$X\sim\chi^2(n),\ Y\sim\chi^2(m)$，X,Y相互独立，则$X+Y\sim\chi^2(m+n)$。

2、t分布  
若$X\sim N(0,1),\ Y\sim\chi^2(n)$，X,Y相互独立，则$\frac X{\sqrt{Y/n}}\sim t(n)$。当$n>30$时，t分布近似标准正态分布。

3、F分布  
若$X\sim\chi^2(n_1),\ Y\sim\chi^2(n_2)$，X,Y相互独立，则$\frac {X/n_1}{Y/n_2}\sim F(n,n2)$。

- 参数估计

1、点估计

2、极大似然估计  
1）总体的概率（离散）或者密度函数（连续）；  
2）构造似然函数$L(\lambda)，\lambda$为参数；  
3）对数似然$\ln L(\lambda)$；  
4）对$\lambda$求导，令导数为0（取得最大值），多个参数时求偏导。

3、区间估计  
区间长度与概率的权衡。  
$P(\hat\theta_1\le\theta\le\hat\theta_2)=1-\alpha，[\hat\theta_1,\hat\theta_2]称为区间，1-\alpha称为置信度$

- 假设检验

若$X\sim N(\mu,\sigma^2)，X_1,X_2,\dots,X_n$是取自$X$的样本，检验水平为$\alpha$  

关于$\mu$的假设检验有，  
双边检验假设：$H_0,\mu=\mu_0;\ H_1,\mu\ne\mu_0$  
单边检验假设：$H_0,\mu\ge\mu_0;\ H_1,\mu\lt\mu_0$或者$H_0,\mu\le\mu_0;\ H_1,\mu\gt\mu_0$

1、U检验  
$\sigma^2已知$
$$U=\frac{\overline X-\mu_0}{\sigma_0/\sqrt n}\sim N(0,1),\ \overline X为样本均值$$
2、T检验  
$\sigma^2未知$
$$T=\frac{\overline X-\mu_0}{S/\sqrt n}\sim t(n-1),\ \overline X为样本均值，S为样本标准差$$

关于$\sigma^2$的假设检验有，  
双边检验假设：$H_0,\sigma^2=\sigma_0^2;\ H_1,\sigma^2\ne\sigma_0^2$  
单边检验假设：$H_0,\sigma^2\ge\sigma_0^2;\ H_1,\sigma^2\lt\sigma_0^2$或者$H_0,\sigma^2\le\sigma_0^2;\ H_1,\sigma^2\gt\sigma_0^2$

3、$\Chi$检验  

$\mu已知$
$$\frac{\sum_{i=1}^n(X_i-\mu_0)^2}{\sigma_0^2}\sim\chi^2(n)$$

$\mu未知$
$$\frac{\sum_{i=1}^n(X_i-\overline X)^2}{\sigma_0^2}\sim\chi^2(n-1)$$

参考：[概率论与数理统计-宋浩讲解](https://www.bilibili.com/video/BV1ot411y7mU?p=1&vd_source=78c8191328a642a8039d336c1a740ac1)
