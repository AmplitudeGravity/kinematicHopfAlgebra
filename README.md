# kinematicHopfAlgebra
This program generates the BCJ numerator in HEFT,YM, ([arxiv:2111.15649](https://arxiv.org/abs/2111.15649),  [arxiv:2208.05519](https://arxiv.org/abs/2208.05519)).  YMS+ $\phi^3$ ([arxiv:2208.05886](https://arxiv.org/abs/2208.05886)), QCD with fermions, YM+ $F^3+F^4$[arxiv:2208.05519](https://arxiv.org/abs/2310.11943), $DF^2$+YM [2403.04614](https://arxiv.org/abs/2403.04614). It can be used in constructing EYM and GR+HEFT amplitude via double copy. 
Examples are included. Local BCJ numerator for YM is also included. 
## Install KiHA
1. Copy "KiHA.wl" to your computer with "mathematica >10" installed
2. Load the package use `Needs["KiHA`","your local directory/KiHA.wl"]`
3. enjoy!


## Heavy-mass effective theory (HEFT) and Yang-Mills (YM)
Related paper for KiHA in HEFT are  [arxiv:2111.15649](https://arxiv.org/abs/2111.15649),  [arxiv:2208.05519](https://arxiv.org/abs/2208.05519) .
In HEFT, kinematic algebra is taken as current algebra. The building blocks are 
  - vector current $T^{(i)}_{(i)}$: correspond to the vector currents and map to a vector current(just the velocity $v$) product with a polarisation vector $\varepsilon_i$ .
  - tensor current $T^{(\alpha)}_{(\tau_1),(\tau_2),\cdots, (\tau_r)}$: Tensor currents and map to all multiplicity  universal tensor product with multi polarisation vectors
  - fusion product $\star$: The fusion rules from lower order rank tensor currents to higher order tensor currents. This fusion product is bilinear and associative.
  - convolution map $\langle \bullet \rangle$: This is a linear map from abstract algebra general to physics kinematic expression, which is, in general, non-local and manifestly gauge invariant. 
      
The major output of this program is 

$$\widehat {\mathcal{N}} (12\ldots n{-}2)={T_{(1)}^{(1)}}\star T_{(2)}^{(2)}\star \ldots \star T_{(n{-}2)}^{(n{-}2)}$$

which is known as an algebraic pre-numerator. Another major function is the convolution map, which maps the abstract generator to the physical function of kinematic information. 
In HEFT, the algebraic pre-numerator is from 

$$\widehat{\mathcal N}(123)= T_{(1)}^{(1)} \star  T_{(2)}^{(2)} \star  T_{(3)}^{(3)}.$$

```
num = \[FivePointedStar] @@ (T /@ List /@ Range[1, 3]);
```
We get

$$T_{\text{(1)},\text{(2)},\text{(3)}}+T_{\text{(1)},\text{(3)},\text{(2)}}+T_{\text{(2)},\text{(1)},\text{(3)}}+T_{\text{(2)},\text{(3)},\text{(1)}}+T_{\text{(3)},\text{(1)},\text{(2)}}+T_{\text{(3)},\text{(2)},\text{(1)}}$$
$$-T_{\text{(1)},\text{(23)}}-T_{\text{(23)},\text{(1)}}-T_{\text{(12)},\text{(3)}}-T_{\text{(3)},\text{(12)}}-T_{\text{(13)},\text{(2)}}-T_{\text{(2)},\text{(13)}}+T_{\text{(123)}}$$

After taking the convolution map, we have 

$${\mathcal{N}}(123,v)=\langle\widehat {\mathcal{N}}(123)\rangle= \langle T_{(1)}^{(1)} \star  T_{(2)}^{(2)} \star  T_{(3)}^{(3)}\rangle.$$ 

```
preNumerator = num /. T -> Tp /. rmzero,
```
where the ```Tp``` is the convolution map. 
Then we get the output of the pre-numerator 

$${\mathcal{N}}(123,v)=-\frac{v\cdot F_1\cdot F_2\cdot v p_{1,2}\cdot F_3\cdot v}{3v\cdot p_1 v\cdot p_{1,2}}-\frac{v\cdot F_1\cdot F_3\cdot v p_1\cdot F_2\cdot v}{3v\cdot p_1 v\cdot p_{1,3}}+\frac{v\cdot F_1\cdot F_2\cdot F_3\cdot v}{3v\cdot p_1}$$

All other BCJ numerators are obtained directly from the BCJ numerator by crossing symmetry. For the $n$ point YM amplitude,  you only need to replace the velocity by the polarisation vector of the last line $\varepsilon_n$

## Yang-Mills-scalar theory
The related paper for KiHA in Yang-Mills-scalar [arxiv:2208.05886](https://arxiv.org/abs/2208.05886). The kinematic algebra is taken as field algebra in the full theory of Yang-Mills-scalar+ $\phi^3$. The building blocks are 
  * vector field ${\mathsf K_i}=T_{(i)}^{(i)}$
  * scalar field ${\mathsf K_j}=T^{(j)}$
  * tensor field $T^{(\alpha)}_{(\tau_1),(\tau_2),\cdots, (\tau_r)}$: fields for multi-particle states lie on the interline, which is all multiplicity university mapping to the gauge invariant functions.
  * fusion product $\star$: The fusion rules from a fewer-particle field to a more-particle field. This fusion product is bilinear and associative.
  * convolution map $\langle \bullet \rangle$: This is a linear map from abstract algebra general to physics kinematic expression. This is the inner product between multi-particle states and single outgoing particle states. For each algebraic generator, the mapping value is in general non-local and manifestly gauge invariant. 
  
The major output of this program is 

$$\widehat {\mathcal N}(1,2,\ldots, n{-}1)= {\mathsf K_1}\star  {\mathsf K_2} \star  \ldots \star {\mathsf K_{n-1}}$$
which is known as an algebraic pre-numerator. Another major function is the convolution map, which maps the abstract generator to the physical function of kinematic information. 
For the amplitude with two scalars, 
```
preNumerator = \[FivePointedStar][\[ScriptCapitalK][1, 
       1], \[ScriptCapitalK][2, 1], \[ScriptCapitalK][3, 0]] /. 
     ET[f__] :> ET2F2s[ET[f]] /. CenterDot[f__] :> tr[f, t^a[n]] /. 
   rmzero //. niceF
```
one can get 

$${\mathcal N}(1,2,\overline 3,\overline 4)=\langle {\mathsf K_1} \star  {\mathsf K_2} \star {\mathsf K_3} \rangle = \langle T_{(1)}^{(1)}\star T_{(2)}^{(2)}\star T^{(3)} \rangle=-\frac{p_3\cdot F_1\cdot F_2\cdot p_3 \text{tr}\left(t^{a_3},t^{a_4}\right)}{p_{3,1}\cdot p_{3,1}}$$

For the amplitude with more than three scalars, 
```
preNumerator = \[FivePointedStar][\[ScriptCapitalK][1, 
      0], \[ScriptCapitalK][3, 1], \[ScriptCapitalK][2, 0]] /. 
    ET[f__] :> ET2F[ET[f]] /. CenterDot[f__] :> tr[f, t^a[n]] //. 
  niceF
```
you get 

$${\mathcal N}(\overline 1,\overline 2,3,\overline 4)=\langle {\mathsf K_1} \star {\mathsf K_2} \star {\mathsf K_3} \rangle= \langle T^{(1)}\star T^{(2)}\star T_{(3)}^{(3)} \rangle=\frac{2 p_1\cdot F_3\cdot p_2 \text{tr}\left(t^{a_1},t^{a_2},t^{a_4}\right)}{p_{1,2}\cdot p_{1,2}}.$$

## BCJ numerator for YM+two massive scalar amplitude
For the Yang-Mills amplitude with two massive scalar, we can use the following more efficient code
```
n = 3;
Timing[numpreP = \[FivePointedStar] @@ (T /@ List /@ Range[1, n]);]
bpeff2 = BinaryProduct[Range[n]]
nod = numpreP /. T -> TScalar /. rmzero
```
See more details in the Notebook "YMSAlgebraNew.nb"

## higher-derivative gauge field theory
The related paper for KiHA in higher-derivative gauge field theory is  [2310.11943](https://arxiv.org/abs/2310.11943) .
We consider the gauge field theory with higher order contraction of the strengthen tensor.

$$ \int \mathrm{d}^D x \text{Tr}\{\frac{1}{4} F_{\mu \nu} F^{\mu \nu}+\frac{2 \alpha^{\prime}}{3} F_\mu^\nu F_\nu^\lambda F_\lambda^\mu+\frac{\alpha^{\prime 2}}{4}[F_{\mu \nu}, F_{\lambda \rho}][F^{\mu \nu}, F^{\lambda \rho}] \} $$

```
num = \[FivePointedStar][T[{1}], T[{2}], T[{3}]] /. rmzeroT /. 
   T[f__] :> T2FF3F4[T[f]];
num=num /. W -> WFun /. F[i__] :> Sequence @@ (F /@ {i});
```
using the nice function
```
ClearAll[nicesp, niceFT]
nicesp = 
  Join[{\[Epsilon][p[i_]] :> \[Epsilon][i], 
    p[f__] :> Subscript[p, StringJoin[ToString /@ {f}]], 
    F[f__] :> Subscript[F, StringJoin[ToString /@ {f}]], 
    CenterDot[Subscript[p, f_], Subscript[p, f_]] :> 
\!\(\*SubsuperscriptBox[\(p\), \(f\), \(2\)]\)}, nice, {v[i_] :> v}];
niceFT = 
  Join[ niceT, {List[f_List, i_Integer] :> 
     StringJoin["(", ToString /@ f, ")"]^ToString[i], 
    FT[f__] :> Subscript[J, f]}];
```

you can see readable form of the pre-numerator
```
num//.nicesp
```
## BCJ numerator in the DF^2+YM theory
This part is to generate the BCJ numerator in the DF2+YM theory. This theory contains massless gluon, massive gluon and tachyon, see [1803.05452](https://arxiv.org/abs/1803.05452)  for reference. 

generate the HEFT BCJ numerator
```
n=7;
num = \[FivePointedStar][T[{1}], T[{2}], T[{3}],T[{4}],T[{5}]] /. rmzeroT /. 
   T[f__] :> T2FF3F4[T[f]];
num = num /. W -> WFun0;
num = num /. W[od__] :> WFunDF2[od] /. repNormalOrder /. 
     W[f__] :> W2basis[W[f]] //. dotRules // Expand;
ng=n-2;
Do[
 Print[Length[wfun]];
 Monitor[
  num = Sum[
     num[[jjj]] /. W -> WFunDF2 /. repNormalOrder /. 
         W[f__] :> W2basis[W[f]] /. W0[f__] :> W02basis[W0[f]] //. 
       dotRules // Expand, {jjj, Length@num}];, jjj];
 , {id, ng - 2}]
```

generate the W prime function
```
ng=5;
wfun = (-1)^(ng - 1) (WFunDF2 @@ Range[ng]) /. repNormalOrder /. 
     W[f__] :> W2basis[W[f]] //. dotRules // Expand;
Do[
 Print[Length[wfun]];
 Monitor[
  wfun = 
    Sum[wfun[[jjj]] /. W -> WFunDF2 /. repNormalOrder /. 
         W[f__] :> W2basis[W[f]] /. W0[f__] :> W02basis[W0[f]] //. 
       dotRules // Expand, {jjj, Length@wfun}];, jjj];
 , {id, ng - 2}]
```
More details can be found in the file "DF2YMExample.nb".


## QCD BCJ numerator
The amplitude with two fermion line and multi gluon lines are also of color-kinematic duality. The kinematic algebra is also quasi-shuffle hopf algebra. 
The pre-numerator is obtained as 
```
n = 5
numJ = \[FivePointedStar] @@ Table[FT[{{i}, 1}], {i, 1, n - 2}] /. 
    rmzeroT /. FT[f__] :> FT2F[FT[f]];
numJ = % /. sp[f__] :> spAB[p[n], f, p[n - 1]] // Expand;
numJ = numJ /. F[i__] :> Sequence @@ (F /@ {i});
```

Another version of the bcj numerator can be obtained recursively.
```
numRecQCD[5] // contractSp
```


## Local BCJ numerator for Yang-Mills
For the local BCJ numerator in Yang-Mills theory,  the kinematic algebra is also quasi-shuffle hopf algebra. 
The BCJ numerator is obtained as 
```
ng = 8;
rmzeroDF2 = {T[{i1_, ils__}, g__] :> 0 /; i1 != 1};
Timing[num2 = (\[FivePointedStar] @@ (Table[
        T[{ii}], {ii, ng - 1}])) /. 
    T[{i1__}, g___] :> T[{i1, ng}, g] /. rmzeroDF2;
 num2 = num2 /. T[f__] :> T2YMLocal[T[f]];]
```

Another version of the bcj numerator can be obtained recursively. We first set the replace rules for the factors in the  W' function
```
repdotnum2 = {W0[1, f___, 
      ng] :> -1/\[Alpha] dot[\[Epsilon][1], 
       Sequence @@ F /@ {f}, \[Epsilon][ng]] /. 
    dot[i1_, F[], i2_] :> dot[i1, i2]};
repdotden2 = {(1 - \[Alpha] dot[p[f__], p[f__]]) :> 
    2 \[Alpha] Length[Complement[Range[ng], {f}]] x};
repdotFnum2 = {dot[f__, F[i_], 
     p[gs___, ng]] :> -dot[f, \[Epsilon][i]] x, 
   dot[f__, F[rs__, i_], 
     p[gs___, ng]] :> -dot[f, F[rs], \[Epsilon][i]] x};
```
Then we use the recursive rules that contribute to the leading order of $\alpha'$ and in local numerator limit, after applying the above replacement rules, we can 
get the local BCJ numerator as "wfunYM"
```
Timing[wfun = (-1)^(ng - 1) (WFun2YM @@ Range[ng]);
 wfun = Sum[
   wfun[[ii]]*(1 - \[Alpha] dot[p @@ Range[ng]]) // Simplify, {ii, 
    Length@wfun}];
 wfun = wfun //. W -> WFun2YM // Expand;
 wfunYM = wfun /. repdotnum2 /. repdotden2 /. repdotFnum2;
 ]
```

