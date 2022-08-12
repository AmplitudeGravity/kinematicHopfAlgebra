# kinematicHopfAlgebra
This is the program to generate the BCJ numerator in HEFT, YM, YMS+ $\phi^3$. It can be used in construct EYM, GR+HEFT amplitdue. 
Examples is included.

The major output of this program is 
$$\widehat N(1,2,\ldots, n{-}1)= K_1 \star  K_1 \star  \ldots \star K_{n{-}1}$$
which is known as algebraic pre-numerator. Another major function is the convolution map which is to map the abstract generator to physical function of kinematic information. 

## Heavy-mass effective theory (HEFT) and Yang-Mills (YM)
In HEFT, 

```
num = \[FivePointedStar] @@ (T /@ List /@ Range[1, 3]);
```

The output is
$$T_{\text{(1)},\text{(2)},\text{(3)}}+T_{\text{(1)},\text{(3)},\text{(2)}}+T_{\text{(2)},\text{(1)},\text{(3)}}+T_{\text{(2)},\text{(3)},\text{(1)}}+T_{\text{(3)},\text{(1)},\text{(2)}}+T_{\text{(3)},\text{(2)},\text{(1)}}$$
$$-T_{\text{(1)},\text{(23)}}-T_{\text{(23)},\text{(1)}}-T_{\text{(12)},\text{(3)}}-T_{\text{(3)},\text{(12)}}-T_{\text{(13)},\text{(2)}}-T_{\text{(2)},\text{(13)}}+T_{\text{(123)}}$$
After taking the convolution map and remove the trivial terms 
```
preNumerator = num /. T -> Tp /. rmzero
```
we get the output of the pre-numerator 
$$-\frac{v\cdot F_1\cdot F_2\cdot v p_{1,2}\cdot F_3\cdot v}{v\cdot p_1 v\cdot p_{1,2}}-\frac{v\cdot F_1\cdot F_3\cdot v p_1\cdot F_2\cdot v}{v\cdot p_1 v\cdot p_{1,3}}+\frac{v\cdot F_1\cdot F_2\cdot F_3\cdot v}{v\cdot p_1}$$

All other BCJ numerator is obtained directly from the BCJ numerator. For the YM you only need to replace the velocity by the polarisation vector of the last line $\varepsilon_n$

## Yang-Mills-scalar theory
For the amplitude with two scalars, 
```
preNumerator = \[FivePointedStar][\[ScriptCapitalK][1, 
       1], \[ScriptCapitalK][2, 1], \[ScriptCapitalK][3, 0]] /. 
     ET[f__] :> ET2F2s[ET[f]] /. CenterDot[f__] :> tr[f, t^a[n]] /. 
   rmzero //. niceF
```
one can get 
$$-\frac{p_3\cdot F_1\cdot F_2\cdot p_3 \text{tr}\left(t^{a_3},t^{a_4}\right)}{p_{3,1}\cdot p_{3,1}}$$

For the amplitude with more than three scalars, 
```
preNumerator = \[FivePointedStar][\[ScriptCapitalK][1, 
      0], \[ScriptCapitalK][3, 1], \[ScriptCapitalK][2, 0]] /. 
    ET[f__] :> ET2F[ET[f]] /. CenterDot[f__] :> tr[f, t^a[n]] //. 
  niceF
```
you get 

$$\frac{2 p_1\cdot F_3\cdot p_2 \text{tr}\left(t^{a_1},t^{a_2},t^{a_4}\right)}{p_{1,2}\cdot p_{1,2}}.$$



# Citation 
If you use **kinematicHopfAlgebra.wl**, please cite the two papers [arxiv:2111.15649](https://arxiv.org/abs/2111.15649),  [arxiv:2208.05519](https://arxiv.org/abs/2208.05519) and [arxiv:2208.05886](https://arxiv.org/abs/2208.05886)as following

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

