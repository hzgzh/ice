# 内燃机工作过程数值模拟

## 公式说明
1. 缸内传热系数计算
    缸内传热采用Woschni公式

$H=820D^{-0.2}p^{0.8}T^{-0.53}[C_1c_m+C_2\frac{V_sT_a}{p_aV_a}(p-p_0)]^{0.8} W/m^2K$

符号说明  
$D:$气缸直径 $m$  
$p:$缸内压力 $MPa$  
$T:$气缸温度 $K$  
$p_a,V_a,T_a:$为进气阀关闭时的初始压力、温度、体积  
$p:$为气缸压力,$p_0:$为倒拖时气缸压力，可通过以下公式求取  
$p_0=(\frac{V_{ca0}}{V_{ca}})^kP_a$  
$V_{ca0}:$压缩开始时的体积，即进气阀关闭点  
$V_{ca}:$缸内体积

2. 计算传热量对曲轴角的微分  
$\frac{dQ}{d\phi}=\frac{1}{6n}h_tA_v\sum\limits_{i=1}^3(T-T_i)$  
$n,A_v,h_t,T_i,T$分别为转速 r/min,面积 $m^2$,换热系数 W/m2/K，缸内壁温，缸内气体温度。  

3. 




