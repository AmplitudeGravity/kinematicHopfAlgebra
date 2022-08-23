# kinematicHopfAlgebra
This is the program to generate the BCJ numerator in HEFT, YM, YMS+ $\phi^3$. It can be used in constructing EYM and GR+HEFT amplitude via double copy. 
Examples are included. 

## Heavy-mass effective theory (HEFT) and Yang-Mills (YM)
In HEFT, the kinematic algebra is taken as current algebra. The building blocks are 
  - vector current $T^{(i)}_{(i)}$: correspond to the vector currents and map to a vector current(just the velocity $v$) product with a polarisation vector $\varepsilon_i$ .
  - tensor current $T^{(\alpha)}_{(\tau_1),(\tau_2),\cdots, (\tau_r)}$: Tensor currents and map to a all multiplicity  universal tensor product with multi polarisation vectors
  - fusion product $\star$: The fusion rules from lower order rank tensor currents to higher order tensor currents. This fusion product is bilinear and associative.
  - convolution map $\langle \bullet \rangle$: This is a linear map from abstract algebra general to physics kinematic expression, which is, in general, non-local and manifestly gauge invariant. 
      
The major output of this program is 

$$\widehat {\mathcal{N}} (12\ldots n{-}2)={T_{(1)}^{(1)}}\star T_{(2)}^{(2)}\star \ldots \star T_{(n{-}2)}^{(n{-}2)}$$

which is known as algebraic pre-numerator. Another major function is the convolution map, which maps the abstract generator to the physical function of kinematic information. 
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

All other BCJ numerators are obtained directly from the BCJ numerator by the crossing symmetry. For the $n$ point YM amplitude,  you only need to replace the velocity by the polarisation vector of the last line $\varepsilon_n$

## Yang-Mills-scalar theory
In the full theory of Yang-Mills-scalar+ $\phi^3$, the kinematic algebra is taken as field algebra. The building blocks are 
  * vector field ${\mathsf K_i}=T_{(i)}^{(i)}$
  * scalar field ${\mathsf K_j}=T^{(j)}$
  * tensor field $T^{(\alpha)}_{(\tau_1),(\tau_2),\cdots, (\tau_r)}$: fields for multi-particle states lie on the interline, which is all multiplicity university mapping to the gauge invariant functions.
  * fusion product $\star$: The fusion rules from a fewer-particle field to a more-particle field. This fusion product is bilinear and associative.
  * convolution map $\langle \bullet \rangle$: This is a linear map from abstract algebra general to physics kinematic expression. This is the inner product between multi-particle states and single outgoing particle states. For each algebraic generator, the mapping value is in general non-local and manifestly gauge invariant. 
  
The major output of this program is 

$$\widehat {\mathcal N}(1,2,\ldots, n{-}1)= {\mathsf K_1}\star  {\mathsf K_2} \star  \ldots \star {\mathsf K_{n-1}}$$
which is known as algebraic pre-numerator. Another major function is the convolution map, which maps the abstract generator to the physical function of kinematic information. 
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



# Citation 
If you use **kinematicHopfAlgebra.wl**, please cite the following three papers [arxiv:2111.15649](https://arxiv.org/abs/2111.15649),  [arxiv:2208.05519](https://arxiv.org/abs/2208.05519) and [arxiv:2208.05886](https://arxiv.org/abs/2208.05886)as following

```
@article{Brandhuber:2021bsf,
    author = "Brandhuber, Andreas and Chen, Gang and Johansson, Henrik and Travaglini, Gabriele and Wen, Congkao",
    title = "{Kinematic Hopf Algebra for Bern-Carrasco-Johansson Numerators in Heavy-Mass Effective Field Theory and Yang-Mills Theory}",
    eprint = "2111.15649",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "NORDITA 2021-091, QMUL-PH-21-45, SAGEX-21-34, UUITP-60/21",
    doi = "10.1103/PhysRevLett.128.121601",
    journal = "Phys. Rev. Lett.",
    volume = "128",
    number = "12",
    pages = "121601",
    year = "2022"
}
```
```
@article{Brandhuber:2022enp,
    author = "Brandhuber, Andreas and Brown, Graham R. and Chen, Gang and Gowdy, Joshua and Travaglini, Gabriele and Wen, Congkao",
    title = "{Amplitudes, Hopf algebras and the colour-kinematics duality}",
    eprint = "2208.05886",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "QMUL-PH-22-18, SAGEX-22-27",
    month = "8",
    year = "2022"
}
```
```
@article{Chen:2022nei,
    author = "Chen, Gang and Lin, Guanda and Wen, Congkao",
    title = "{Kinematic Hopf algebra for amplitudes and form factors}",
    eprint = "2208.05519",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "QMUL-PH-22-22",
    month = "8",
    year = "2022"
}
```

