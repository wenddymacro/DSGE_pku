* 北大经济学院《数量经济学》Guest Lecture：DSGE与编程导论
* 经典货币经济RBC模型
* stata代码 
* @许文立，2022-3-31


* 加载数据
use "/Users/xuwenli/Desktop/北京大学-DSGE与编程/example_code/usmacro2.dta",clear

* 校准参数
constraint 1 _b[alpha]=0.33        // 资本收入份额
constraint 2 _b[delta]=0.025      
constraint 3 _b[beta]=0.8
constraint 4 _b[Pibar]=1
constraint 5 _b[phi_pi]=1.5
constraint 6 _b[rho]=0.9



* 声明模型
dsgenl (1/c={beta}/F.c*({alpha}*F.A*(F.k)^({alpha}-1)+(1-{delta})))  ///
       (1/c={beta}/F.c*(r/F.p))    ///
	   (A*k^{alpha}=c+F.k-(1-{delta}*k))     ///
	   (y=A*k^{alpha})     ///
	   (r*{beta}=(p/{Pibar})^{phi_pi})      ///
	   (ln(F.A)={rho}*ln(A))      ///
	   ,observed(p) unobserved(c y r) exostate(A) endostate(k)  ///
	   constraint(1/6)

* 政策函数
estat policy


* 转移函数
estat transition

* 脉冲响应
irf set rbc_cme.irf

irf create model1

irf graph irf, impulse(A) response(y c k A) byopts(yrescale)





