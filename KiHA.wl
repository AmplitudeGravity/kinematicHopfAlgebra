(* ::Package:: *)

(*Copyright 2022,Gang Chen AmplitudeGravity/kinematicHofpAlgebra is licensed under the
GNU General Public License v3.0 .  Permissions of this strong copyleft license are 
conditioned on making available complete source code of licensed works and modifications,
which include larger works using a licensed work,under the same license.
Copyright and license notices must be preserved.Contributors provide an express grant 
of patent rights.*)
(*KiHA is dependent on some basis functions which are written by Gustav Mogull and Gregor Kaelin in Uppsala theoritical physics group 
*)(*  Updating version can be find in  
https://github.com/AmplitudeGravity/kinematicHofpAlgebra
  *)


BeginPackage["KiHA`"]
Unprotect@@Names["KiHA`*"];
(* Exported symbols added here with SymbolName::usage *) 


If[$Notebooks,
  CellPrint[Cell[#,"Print",CellFrame->0.5,FontColor->Blue]]&,
  Print][" KiHA(v5.6),Copyright 2022,author Gang Chen. It is licensed under the GNU General Public License v3.0. 
 KiHA is based on the work of Kinematic Hopf Algebra in CTP of Queen Mary University of London. It generates the duality all-n numerator for colour-kinematic duality and double copy in heavy mass effective 
 theory(HEFT), YM, YM-Scalar, YM-fermion, F3+F4 theory and DF2+YM theory. KiHA is built on some basic functions written by Gustav Mogull and Gregor Kaelin.
 use ?KiHA`* for help and a series of papers (2111.15649, 2208.05519, 2208.05886,2310.11943) for more reference." ]


(* ::Subsection:: *)
(*KiHA Output function and variables*)


nice::usage ="nice the output form, //.nice"
niceT::usage ="nice the output form, //.niceT"
rmzero::usage ="remove the trivial terms dot prodoct functions"
rmzeroT::usage ="remove the trivial terms in T generator"
BinaryProduct::usage ="the binary products in a given order"
GBinaryProduct::usage ="all the binary products"
X::usage="Empty bracket"
ET::usage="extended generator of GT with flavor ET[GT[{1,3,2},{{3}}],{a[1],a[2]}]"
FT::usage="extended generator of GT with flavor FT[{{1,2},3}]"
GT::usage="extended generator of T with colour group GT[{1,3,2},{{3}}], first arg is colour order, sencond arg is the kinematic partitions"
T::usage="algebra generator, e.g. T[{1},{2,3}],"
\[DoubleStruckA]::usage="generator index of the flavor group"
p::usage="momentum vector"
\[Epsilon]::usage="polarization vector"
F::usage="strengthen tensor"
m::usage="mass of the field,"
\[DoubleStruckT]::usage="generator of the flavor group"
v::usage="velocity of the heavy particle or black hole"
T2F::usage ="transform abelian generators to strengthen tensor form in HEFT amplitude"
GT2F::usage ="transform non-abelian generators to strengthen tensor form in HEFT amplitude"
ET2F::usage ="transform generators to strengthen tensor form in amplitude with three more scalars"
ET2F2s ::usage ="transform generators to strengthen tensor form in amplitude with two scalars"
\[FivePointedStar] ::usage ="funsion product of the generators"
\[ScriptCapitalK] ::usage ="generators for the external states \[ScriptCapitalK][stateIndex_Integer, Spin_Integer]"
\[CapitalOmega] ::usage ="left nested commutators"
\[DoubleStruckCapitalS] ::usage ="antipode"
\[CapitalDelta] ::usage ="coproduct"
\[DoubleStruckCapitalI] ::usage ="Identity"
sp::usage="spin chain product"
Diamond::usage="a distributive general product"


KSetPartitions::usage="KSetPartitions[list,number of subsets]"
KSetPartitionsEmpty::usage="including empty sets KSetPartitionsEmpty[list,number of subsets]"
SetPartitionsOrdered::usage="SetPartitionsOrdered[number of subsets]"


ddmBasis::usage ="generate the ddm basis, ddmBasis[n]"
KLTR::usage ="generate the one elemment of the KLT matrix KLTR[f_List]"
KLTRM::usage ="generate the first row of the KLT matrix, e.g. KLTRM[5] generate the KLT matrix element for the colour order{1,2,3,4,5} "


ExpandNCM::usage ="Expand the non-commutative multiply"
NC::usage="commutator of the non-commutative multiply"


(* ::Subsection::Closed:: *)
(*QCD current*)


numRecQCD::usage="BCJ numerator for QCD amplitude numRecQCD[n] or numRecQCD[{1,2,3...}] "


NumRec::usage="BCJ numerator for QCD amplitude NumRec[{1,2,3...}] "


FT2F::usage="transform fermion algebra to F form "
FRank2::usage="G2 function "


JQCD::usage="Generate the qcd current from current algebra JQCD[{1,2,3,4}]"
X2Prop::usage="The propagators from the binary product. The last line is taken as on-shell"


(* ::Subsection:: *)
(*alpha' higher orders*)


T2FF3F4::usage="transform F^3+F^4 kinematic algebra to F form"
T2YMLocal::usage="transform YM local kinematic algebra to F form "
WFun::usage="W prime function in F^3+F^4  evaluation map"
W::usage="Abstract W prime function"
\[Alpha]::usage="string tension constant"
W0::usage="Abstract W0 function"
WFunDF2::usage="W prime function in the DF^2+YM evaluation map,WFunDF2[1,2,3,4,5],WFunDF2[1,2,3,4,5,{{{1,2},{3,4,5}}}]"
MultiTrace::usage="generate terms in the BCJ numerator corresponding to the multi trace partitions, e.g. MultiTrace[5]
MultiTrace[1,2,3,4,5],MultiTrace[5,{{{1,2},{3,4,5}}}],MultiTrace[1,2,3,4,5,{{{1,2},{3,4,5}}}]"
WFunDF2SingleTr::usage="W prime single trace function in the DF^2+YM evaluation map,WFunDF2SingleTr[1,2,3,4,5]"
SingleTrace::usage="generate single trace terms in the BCJ numerator corresponding to the multi trace partitions, e.g. SingleTrace[5]
SingleTrace[1,2,3,4,5]"
W0Relations::usage="generate all the relations of W0 or W prime function, e.g. W0Relations[5],W0Relations[1,2,3,4,5]"
repW0::usage="replace all the W0 function to their basis"
repWp::usage="replace all the W prime function to their basis"
trFLA::usage="value of W0 function, e.g. trFLA[1,2,3,4,5]"
WFun0::usage="value of W function, e.g. WFun0[1,2,3,4,5]"
W2basis::usage="efficient replace all the W prime function to their basis,W2basis[W[2,1,3,5,4]]"
W02basis::usage="efficient replace all the W0 function to their basis"
WFun2YM::usage="Generate the W function for pure Yang-Mills BCJ numerator, e.g WFun2YM[1,2,3,4,5]"


(* ::Subsection:: *)
(*Scalars*)


setGlobalDim::usage = "setGlobalDim[d] sets d as the global dimension."
setGlobalDim::ssle = "Protected objects cannot be set as the global dimension."
setGlobalDim::author = "Gustav Mogull"
declareAntisymmetric::usage="set function is antisymmetric under exchange the order of the variables"
declareSymmetric::usage="set function do not dependent order of the variables"


(* ::Subsection::Closed:: *)
(*Vectors*)


tensorQ::usage = "tensorQ[p] yields True if p is a tensor, and yields False otherwise."
tensorQ::author = "Gustav Mogull"


vectorQ::usage = "vectorQ[p] yields True if p is a vector (rank-1 tensor), and yields False otherwise."
vectorQ::author = "Gustav Mogull"
tensorQR2::usage = "vectorQ[p] yields True if p is a vector (rank-1 tensor), and yields False otherwise."
tensorQR2::author = "Gustav Mogull"


antiTensorQ::usage="antiTensorQ[p[1]] yields True if p is tensor head"


tensorDim::usage = "tensorDim[p] returns the dimensionality of the declared vector p."
tensorDim::argx = "tensorDim should only be called on tensors."
tensorDim::author = "Gustav Mogull"


tensorRank::usage = "tensorRank[p] yields the rank of the declared tensor p."
tensorRank::argx = "tensorRank should only be called on tensors."
tensorRank::author = "Gustav Mogull"


declareVector::usage = "declareVector[p] declares p as a vector (rank-1 tensor).
declareVector[{p1,p2,...}] declares p1, p2, ... as vectors.
declareVector has options 'dim' and 'verbose'.  By default, 'dim' is the global dimension and 'verbose' is true."
declareVector::safe = "Only unprotected symbols can be declared as vectors."
declareVector::author = "Gustav Mogull"


declareVectorHead::usage = "declareVectorHead[p] declares that all objects with header p should be treated as vectors (rank-1 tensors).
declareVectorHead[{p1,p2,...}] declares p1, p2, ... as vector headers.
declareVectorHead has options 'dim' and 'verbose'.  By default, 'dim' is the global dimension and 'verbose' is true."
declareVectorHead::safe = "Only unprotected symbols can be declared as vector headers."
declareVectorHead::author = "Gustav Mogull"


declareTensor::usage = "declareTensor[F] declares F as a tensor.
declareTensor[{F1,F2,...}] declares F1, F2, ... as tensors.
declareTensor has options 'dim', 'rank' and 'verbose'.  By default, 'dim' is the global dimension, 'rank' is 1 and 'verbose' is true."
declareTensor::safe = "Only unprotected symbols can be declared as tensors."
declareTensor::author = "Gustav Mogull"


declareTensorHead::usage = "declareTensorHead[p] declares that all objects with header p should be treated as tensors.
declareTensorHead[{F1,F2,...}] declares F1, F2, ... as tensor headers.
declareTensorHead has options 'dim', 'rank' and 'verbose'.  By default, 'dim' is the global dimension, 'rank' is 1 and 'verbose' is true."
declareTensorHead::safe = "Only unprotected symbols can be declared as tensor headers."
declareTensorHead::author = "Gustav Mogull"


declareAntiTensorHead::usage="declareAntiTensorHead[p], declareAntiTensorHead[{F1, F2}] declare the antiSymmetric tensor "


undeclareTensor::usage = "undeclareTensor[p] undeclares p as a vector.
undeclareTensor has the option 'verbose', by default taken as True."
undeclareTensor::author = "Gustav Mogull"


undeclareTensorHead::usage = "undeclareTensorHead[p] undeclares p as a tensor header.
undeclareTensorHead has the option 'verbose', by default taken as True."
undeclareTensorHead::author = "Gustav Mogull"


dot::usage = "dot[p1,p2] represents the symmetric Lorentz inner product between the two vectors p1 and p2.
dot[p1,F1,F2,...,p2] represents an inner product where p1, p2 are vectors and F1, F2, ... are rank-2 tensors.
dot[p] evaluates as dot[p,p]."
dot::badrank = "The first and last arguments of dot should be rank-1 tensors (vectors); any inbetween should be rank-2 tensors."
dot::author = "Gustav Mogull"


tr::usage = "tr[F1,F2,...] represents the trace of rank-2 tensors F1, F2, ...."
tr::author = "Gustav Mogull"


eps::usage = "eps[p1,p2,...,pd] represents the integer d-dimensional Levi-Civita symbol acting on vectors p1, p2, ..., pd."
eps::author = "Gustav Mogull"


outer::usage = "outer[p1,p2] is a rank-2 tensor representing the outer product between two vectors p1, p2."
outer::author = "Gustav Mogull"


eta::usage = "eta is a rank-2 tensor representing the flat (mostly-minus) metric."
eta::author = "Gustav Mogull"


mu::usage = "mu[l1,l2] represents the symmetric Lorentz inner product between the extra-dimensional parts of the two vectors l1 and l2.
mu[l] evaluates as mu[l,l]."
mu::argx = "mu called with `1` arguments; 1 or 2 arguments are expected."
mu::badrank = "The arguments of mu should be rank-1 tensors (vectors)."
mu::author = "Gustav Mogull"


aMu::usage = "aMu[l1,l2] represents the antisymmetric product of the extra-dimensional part of two vectors l1 and l2"
aMu::author = "Gregor Kaelin"


(* ::Subsection::Closed:: *)
(*Spinors*)


spQ::usage = "spQ[p] yields True if p is a spinor, and yields False otherwise."
spQ::author = "Gustav Mogull"

spA::usage = "spA[p] represents the angle spinor of a vector p."
spB::usage = "spB[p] represents the square spinor of a vector p."
spA::author = "Gustav Mogull"
spB::author = "Gustav Mogull"

spAA::usage = "spAA[p1,p2,...,pn] yields the spinor bracket sp[spA[p1],p2,...,spA[pn]]."
spAB::usage = "spAB[p1,p2,...,pn] yields the spinor bracket sp[spA[p1],p2,...,spB[pn]]."
spBA::usage = "spBA[p1,p2,...,pn] yields the spinor bracket sp[spB[p1],p2,...,spA[pn]]."
spBB::usage = "spBB[p1,p2,...,pn] yields the spinor bracket sp[spB[p1],p2,...,spB[pn]]."
spAA::author = "Gustav Mogull"
spBB::author = "Gustav Mogull"
spAB::author = "Gustav Mogull"
spBA::author = "Gustav Mogull"


spp::usage = "spp[sp1,sp2] represents an antisymmetric inner product between the two spinors sp1 and sp2 (angle or square bracket).
spp[sp1,p1,p2,...,sp2] represents a spinor product with vectors p1, p2, ... inserted."
spp::author = "Gustav Mogull"


spOuter::usage = "spOuter[sp1,sp2] is a vector formed from the outer product of an angle and a square spinor.
spOuter[p] returns spOuter[spA[p],spB[p]]."
spOuter::notvector = "spOuter with one argument should only be applied to vectors."
spOuter::badspinors = "spOuter should act on one angle and one square spinor."
spOuter::argx = "spOuter called with `1` arguments; 1 or 2 arguments are expected."
spOuter::author = "Gustav Mogull"


toSpinors::usage = "toSpinors[expr,vecs] converts all instances of the four-dimensional vectors vecs in expr to spinors."
toSpinors::ssle = "Only four-dimensional vectors can be converted to spinors."
toSpinors::author = "Gustav Mogull"


(* ::Subsection::Closed:: *)
(*Free Indices*)


lI::usage = "lI[i] represents the i'th Lorentz (spacetime) index."
lI::author = "Gustav Mogull"


ui::usage = "ui[j] represents the i'th general metric  up index."
ui::author = "Gang Chen"


di::usage = "di[j] represents the i'th general metric down index."
di::author = "Gang Chen"


(* spI::usage = "spI[i] represents the i'th spinor index."
spI::author = "Gustav Mogull" *)


contract::usage = "contract[expr] contracts all pairs of raised and lowered Lorentz indices in expr."
contract::author = "Gustav Mogull"


contractAnti::usage = "contractAnti[expr] contract in expr with the anti-symmetric tensor or vectors"
contractAnti::author = "Gang Chen"


contractspp::usage = "contractspp[expr] contract in expr with index in spp"
contractspp::author = "Gang Chen"


contractSp::usage = "contractSp[expr] nice form of the contract in expr with index in spp"
contractSp::author = "Gang Chen"


expandTensor::usage="expandTensor[dot[f]] expands all the tensor into the conponent"
expandTensor::author="Gang Chen"


dressIndex::usage = "dressIndex[expr,i] dresses the expression expr with the free Lorentz index i."
dressIndex::author = "Gustav Mogull"


freeIndices::usage = "freeIndices[expr] finds the overall free index structure in expr."
freeIndices::author = "Gustav Mogull"


exposeIndex::author = "Gustav Mogull"


indexCoefficient::author = "Gustav Mogull"


(* ::Subsection::Closed:: *)
(*Other*)


parkeTaylor::usage = "parkeTaylor[p1,p2,...,pn] gives the Parke-Taylor factor for n four-dimensional vectors p1, p2, ..., pn."
parkeTaylor::ssle = "parkeTaylor should only be called on four-dimensional vectors."
parkeTaylor::author = "Gustav Mogull"


basisExpand::author = "Gustav Mogull"


lgScaling::usage = "lgScaling[expr,vec] counts the little-group scaling of expr with respect to the on-shell vector vec."
lgScaling::author = "Gustav Mogull"


massDim::usage = "massDim[expr] counts the mass dimension of expr."
massDim::author = "Gustav Mogull"


protectedQ::usage = "protectedQ[expr] checks whether expr is an unprotected symbol."
protectedQ::author = "Gustav Mogull"
(* safeVariableQ::usage = "safeVariableQ[expr] checks whether expr is an unprotected symbol."
safeVariableQ::author = "Gustav Mogull" *)
distributive::usage = "distributive is an option used for distributive functions to determine their behavior."
distributive::author = "Gustav Mogull"
verbose::usage = "verbose is an option used to determine whether certain functions give printed output statements."
verbose::author = "Gustav Mogull"
declareDistributive::usage = "declareDistributive[s,test] declares that the symbol s should linearly distribute its arguments over objects satisfying test.
declareDistributive[s,test, n1, n2] specifies that the first n1 and last n2 arguments of s should be excluded.
declareDistributive[s,test, n1, n2, default] specifies that the default value of distributive should be default."
declareDistributive::author = "Gustav Mogull, Gregor Kaelin"
declareDistributiveMiddle::author = "Gustav Mogull, Gregor Kaelin"
declareDistributiveFirst::author = "Gustav Mogull, Gregor Kaelin"
declareDistributiveLast::author = "Gustav Mogull, Gregor Kaelin"
distributiveRules::author = "Gustav Mogull, Gregor Kaelin"
patternFreeQ::usage = "Checks if an expression is free of any any pattern (i.e. Pattern, Blank, BlankSequence, BlankNullSequence)"
patternFreeQ::author = "Gregor Kaelin"


(* ::Section:: *)
(*Private*)


Begin["`Private`"]


(*declareVectorHead[{p}];
declareTensorHead[{F},{"rank"-> 2}];*)


(* ::Subsection::Closed:: *)
(*Simplifications*)


nice={p[i__]:> Subscript[p, i],F[i_]:> Subscript[F, i],\[DoubleStruckA][i_]:> Subscript[\[DoubleStruckA], i],
dot[f___]:> CenterDot[f],CenterDot[p[i_],\[Epsilon][j_]]:>CenterDot[\[Epsilon][j],p[i]],
a_[i_]?vectorQ:> Subscript[a, i],a_[i_]?tensorQ:> Subscript[a, i],spp-> Diamond,J[f1_,f2_,f3_]:> Subscript[J, f2],J[f1_,f2_]:> Subscript[J, f2]};
(*niceF={dot[f_]:> CenterDot[f,f],dot[f__]:> CenterDot[f],p[i__]:> Subscript[p, i],F[i_]:> Subscript[F, i],a[i_]:> Subscript[a, i]};*)


(*niceT={T[f___]:>(Subscript[T, f]/. List->L),L[f1___]:>"("<>ToString/@{f1}<>")"}*)
niceT={T[f___]:>(Subscript[T, f]/. List->L),\[DoubleStruckA][i_]:>\[DoubleStruckT]^Subscript[\[DoubleStruckA], i],GT[{f1___},{f2___}]:>(\!\(\*SuperscriptBox[
SubscriptBox[\(T\), \(f2\)], \("\<(\>" <> ToString /@ {f1} <> "\<)\>"\)]\)/. List->L),L[f1___]:>"("<>ToString/@{f1}<>")",ET[f1_GT,{f2__}]:>f1 tr[CenterDot[f2]]}


rmzero={dot[a_,F[i_],a_]:> 0,dot[p[],f__,v]:> 0,dot[p[],g___]:> 0,dot[g___,p[]]:> 0};
rmzeroT={T[{i_},g___]:> 0,T[f__,{1,h___},g___]:>0,FT[f1_,g__]:>0/;(Length[f1[[1]]]==1||f1[[1,1]]=!=1)};


SetAttributes[Diamond,{NHoldAll}];
declareDistributive[Diamond,tensorQ];
declareDistributive[Diamond,vectorQ];


(*Generating the binary product*)


(* ::Subsection::Closed:: *)
(*BinaryProduct and DDM*)


shuffle[s1_List,s2_List]:=shuffle[s1,s2]=Module[{p,tp,ord},p=Permutations@Join[1&/@s1,0&/@s2]\[Transpose];
tp=BitXor[p,1];
ord=Accumulate[p] p+(Accumulate[tp]+Length[s1]) tp;
Outer[Part,{Join[s1,s2]},ord,1][[1]]\[Transpose]]


BinaryProduct[n_]:=Block[{v,bpX,bpvars2,h,h0,tt},
v=Range[n];
If[n>=2,bpX={X[v[[1]],v[[2]]]};bpvars2={{v[[2]]}},Return[v[[1]]]];
Do[h={};
	h0={};(*Print[i];*)
 tt=addOne[v[[i]],bpX,bpvars2];bpX=tt[[1]];bpvars2=tt[[2]],{i,3,Length[v]}];
Return[bpX]
]
addOne[vi_,bpX_,bpvars2_]:=Block[{sexpr,vars2,lv2,newvars,h,h0},h={};h0={};Do[(*Print[ii];*)sexpr=bpX[[ii]];vars2=bpvars2[[ii]];lv2=Length[vars2];
Do[h=Join[h,{sexpr/.vars2[[jj]]-> X[vars2[[jj]],vi]}];
newvars=Table[vars2[[k]]/.vars2[[jj]]-> X[vars2[[jj]],vi],{k,jj}];
newvars=Join[newvars,{vi}];
h0=Join[h0,{newvars}],{jj,lv2}];
h=Join[h,{X[sexpr,vi]}];
h0=Join[h0,{{vi}}],
{ii,Length[bpX]}];
(*Print[h,h0];*)
Return[{h,h0}]
]


BinaryProduct[f_List]:=BinaryProduct[Length[f]]/.i_/;IntegerQ[i]:> f[[i]]


GBinaryProduct[h_List/;Length[h]>2]:=Module[{pts,res,r1,r2},pts=partition[h]//Cases[f_/;Length[f]==2];(*pts=Drop[pts,{-1}];*)Table[r1=(GBinaryProduct[pts[[ii,1]]]);r2=(GBinaryProduct[pts[[ii,2]]]);res=Table[X[r1[[ii1]],r2[[ii2]]],{ii1,Length@r1},{ii2,Length@r2}]//Flatten,{ii,Length@pts}]//Flatten]
GBinaryProduct[h_List/;Length[h]==2]:={X[h[[1]],h[[2]] ]}
GBinaryProduct[h_List/;Length[h]==1]:={h[[1]]}


diagNumber[n_]:=1/(n-2+1)*Binomial[2*(n-2),n-2]


(*Programms on the closed form and convolution map*)


ddmBasis[momlist_List] /; If[Length[momlist]>=3,True,Message[ddmBasis::notenoughpoints];False]:=Block[{},
If[!vectorQ[#],declareTensor[#]]&/@momlist;
Prepend[#,momlist[[1]]]&/@(Append[#,momlist[[-1]]]&/@Permutations[momlist[[2;;-2]]])
]


replace[otherDDMOrder_,n_]:=Table[(p/@Range[n-1])[[i]]-> otherDDMOrder[[i]],{i,n-1}];


KLTR[f_List]:=KLTR[f]=Which[Length[f]>2,Module[{imax=Max[f],idmax,fr,fl},idmax=Position[f,imax]//Flatten;fr=Drop[f,idmax];
fl=f[[1;;(idmax[[1]]-1)]];2dot[p[imax],(p/@fl)//Total]KLTR[fr]],Length[f]==2,2dot@@(p/@f),Length[f]<2,Print["Length of f should be large than 2"]]
KLTRM[n_]:=KLTR/@(Drop[#,-1]&/@ddmBasis[Range[n]])


(* ::Subsection:: *)
(*SetPartitions*)


KSetPartitions[0,0]:={{}}
KSetPartitions[{},0]:={{}}
KSetPartitions[s_List,0]:={}
KSetPartitions[s_List,k_Integer]:={}/;k>Length[s]
KSetPartitions[s_List,k_Integer]:={({#1}&)/@s}/;k===Length[s]
KSetPartitions[s_List,k_Integer]:=Block[{$RecursionLimit=Infinity,j},Join[(Prepend[#1,{First[s]}]&)/@KSetPartitions[Rest[s],k-1],Flatten[(Table[Prepend[Delete[#1,j],Prepend[#1[[j]],s[[1]]]],{j,Length[#1]}]&)/@KSetPartitions[Rest[s],k],1]]]/;k>0&&k<Length[s]
KSetPartitions[0,(k_Integer)?Positive]:={}
KSetPartitions[(n_Integer)?Positive,0]:={}
KSetPartitions[(n_Integer)?Positive,(k_Integer)?Positive]:=KSetPartitions[Range[n],k]
NSetPartitions[s_List,n_Integer]:=Select[SetPartitions[s],Union[Length/@#]=={2}&];
KSetPartitionsEmpty[s_List,k_Integer]:=Module[{app,KMod},app=KSetPartitions[s,k];
Do[KMod=KSetPartitions[s,q];
Do[Do[AppendTo[KMod[[i]],{}],{j,1,k-q}];
AppendTo[app,KMod[[i]]],{i,Length[KSetPartitions[s,q]]}],{q,1,k-1}];Return[app]];
KSetPartitionsEmptyPerm[s_List,k_Integer]:=Module[{app,KMod,app2={}},app=KSetPartitions[s,k];
Do[KMod=KSetPartitions[s,q];
Do[Do[AppendTo[KMod[[i]],{}],{j,1,k-q}];
AppendTo[app,KMod[[i]]],{i,Length[KSetPartitions[s,q]]}],{q,1,k-1}];
Do[AppendTo[app2,Permutations[app[[ii]]]],{ii,Length[app]}];
app2=Flatten[app2,1]//DeleteDuplicates;Return[app2]];
SetPartitionsOrdered[1]={{{1}}};
SetPartitionsOrdered[n_Integer]:=Module[{app={},previous=SetPartitionsOrdered[n-1]},Do[app=Join[app,{Append[previous[[ii,1;;-2]],Append[previous[[ii,-1]],n]],AppendTo[previous[[ii]],{n}]}],{ii,Length[previous]}];
Return[app]];


(* ::Subsection::Closed:: *)
(*Algebraic generators to dot product of F*)


(*Tt[od_]:=Module[{phat=Range[Length@od],leftv,num,den},phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],#<od[[i,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],v]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=(-1)^Length@od dot[v,p[1]]Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]*)
T2F[{i_}]:=dot[\[Epsilon][i],v]
T2F[f__]:=Module[{od={f},phat,leftv,num,den},phat=Range[Length@od];
(*If[od===T[{i_}]];*)
phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],#<od[[i,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],v]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=dot[v,p[1]]Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]

GT2F[{g__},{{i_}}]:=dot[\[Epsilon][i],v]
GT2F[f1_List,f2_List]:=Module[{cod=f1,od=f2,phat,leftv,num,den},phat=Range[Length@od];phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],Position[cod,#][[1,1]]<Position[cod,od[[i,1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],v]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=dot[v,p[MinimalBy[od//Flatten,Position[cod,#1]&]//First]]Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]
GTRM0[f1_List,f2_List]:=Module[{cod=f1,od=f2,phat},phat=Range[Length@od];phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],Position[cod,#][[1,1]]<Position[cod,od[[i,1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];
(*Print[phat];*)
If[MemberQ[phat,p[]],0,GT[f1,f2]]
]

TScalar[{i_}]:=dot[\[Epsilon][i],p[0]]
TScalar[f__]:=Module[{od={f},phat,leftv,rightv,num,den},phat=Range[Length@od];
(*If[od===T[{i_}]];*)
phat[[1]]=p[0];Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],#<od[[i,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
rightv=Range[Length@od];
rightv[[1]]=p[0];
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)rightv[[i]]=p[0]+p@@Select[Flatten[od[[1;;(i-1)]]],#>od[[i,-1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],rightv[[i]]]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=1/2 dot[p[0]+p[1]]Product[1/2 dot[p[0]+p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]
ET2F[f_GT,flavors_List]:=Module[{flavor=flavors,od=f[[2]],cod=f[[1]],phat,leftv,rightv,num,den,scalars},phat=Range[Length@od];
scalars=Complement[Flatten[cod],Flatten[od]];
(*If[od===T[{i_}]];*)
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[Join[scalars,od[[1;;(i-1)]]]],Position[cod,#][[1,1]]<Position[cod,od[[i,1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,1,Length@od}];leftv=phat;
rightv=Range[Length@od];
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)rightv[[i]]=p@@Select[Flatten[Join[scalars,od[[1;;(i-1)]]]],Position[cod,#][[1,1]]>Position[cod,od[[i,-1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,1,Length@od}];
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],rightv[[i]]]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=1/2 dot[p@@scalars]Product[1/2 dot[(p@@scalars)+p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
den=den/.dot[gg__]:>dot[gg]-m^2;
flavor=CenterDot@@(flavor/.\[DoubleStruckA][i_]:>\[DoubleStruckT]^\[DoubleStruckA][i]);
 flavor num/den
]
ET2F2s[f_GT,flavors_List]:=Module[{flavor=flavors,od=f[[2]],cod=f[[1]],phat,leftv,rightv,num,den,refs,sc,refg},leftv=Range[Length@od];
refs=Join[Complement[Flatten[cod],Flatten[od]],{Min[Flatten[od]]}];
sc=First@Complement[Flatten[cod],Flatten[od]];
refg=Min[Flatten[od]];
If[Not[(cod[[1]]==sc&&cod[[-1]]==refg)||(cod[[1]]==refg&&cod[[-1]]==sc)],num=0;Return[num]];
(*If[od===T[{i_}]];*)
leftv[[1]]=p@@Complement[Flatten[cod],Flatten[od]];
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)leftv[[i]]=p@@Select[DeleteDuplicates[Flatten[Join[refs,od[[1;;(i-1)]]]]],Position[cod,#][[1,1]]<Position[cod,od[[i,1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];
rightv=Range[Length@od];
rightv[[1]]=p@@Complement[Flatten[cod],Flatten[od]];
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)rightv[[i]]=p@@Select[DeleteDuplicates[Flatten[Join[refs,od[[1;;(i-1)]]]]],Position[cod,#][[1,1]]>Position[cod,od[[i,-1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],rightv[[i]]]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,1,Length@od}];
den=dot[p@@refs]Product[1/2 dot[(p@@DeleteDuplicates[Flatten[Join[refs,od[[1;;(i-1)]]]]])],{i,2,Length@od}];
den=den/.dot[gg__]:>dot[gg]-m^2;
flavor=CenterDot@@(flavor/.\[DoubleStruckA][i_]:>\[DoubleStruckT]^\[DoubleStruckA][i]);
 flavor num/den
]


FT2F[f_FT]:=Module[{fls,od,odl,n,phat,leftv,rightv,num,num1,den,refs,sc,refg},
fls=List@@f;
od=fls[[All,1]];
n=2+Length@(od//Flatten);
leftv=Range[Length@od];
     leftv[[1]]=0;
rightv=Range[Length@od];
     rightv[[1]]=0;
Table[
odl=Flatten[od[[1;;(i-1)]]];
leftv[[i]]=p@@Select[odl,#<od[[i,1]]&];
rightv[[i]]=p@@Join[{n},Select[odl,#>od[[i,-1]]&]],{i,2,Length@od}];
num=dot[p[n],F[Sequence@@od[[1]]],lI[1]]sp[lI[1]]+1/4 dot[p[1],F[Sequence@@Rest[od[[1]]]],lI[1]] sp[F[1], lI[1]];
Do[If[fls[[ii,2]]==2,num=num/.sp[gg__]:>sp[gg,F[leftv[[ii]],od[[ii]]]],num=num/.sp[gg__]:>2dot[leftv[[ii]],F[Sequence@@od[[ii]]],rightv[[ii]]]sp[gg]],{ii,2,Length@fls}];
den=dot[p[n,1]]Product[ dot[(p@@Flatten[Join[{n},od[[1;;(i-1)]]]])],{i,2,Length@od}];
den=den;
  (num/den)//.sp[f1___,F[p[g__],{i_,h___}],f2___]:>FRank2[p[g],{i,h}]sp[f1,lI[i],lI[i,0],f2]
]
(*FRank2[pv_,{i_}]:=dot[pv,F[i],lI[i]]p[i][lI[i,0]]
FRank2[pv_,od_List]:=FRank2[pv,od]=Module[{fewF,newF1,newF2},fewF=FRank2[pv,od[[1;;-2]]]//Expand;newF1=fewF/.dot[f__,lI[i_,0]]:>dot[f,F[od[[-1]]],lI[i,0]]/.p[ii__][lI[i_,0]]:>dot[p[ii],F[od[[-1]]],lI[i,0]];
newF2=fewF/.dot[h__,lI[i_]]dot[p[f__],g__,lI[i_,0]]:>0(*/.dot[h__,lI[i_]]dot[p[f__],g__,lI[i_,0]]:>dot[h,Sequence@@{First[{g}]},lI[i]]dot[p[f],Sequence@@Rest[{g}],F[od[[-1]]],lI[i,0]]*)/.dot[h__,lI[i_]]p[ii__][lI[i_,0]]:>dot[h,F[od[[-1]]],lI[i]]dot[p[ii,od[[-1]]],lI[i,0]];
newF1+newF2
]*)
FRank2[pv_,od_List]:=Module[{pts,lv,res},pts=KSetPartitionsEmpty[od,2];res=Sum[lv=p@@Select[pts[[ii,1]],#<pts[[ii,2,1]]&];dot[pv,F@@pts[[ii,1]],lI[od[[1]]]]dot[lv,F@@pts[[ii,2]],lI[od[[1]],0]],{ii,-1+Length@pts}];
res=res+dot[pv,F@@od,lI[od[[1]]]]dot[p@@od,lI[od[[1]],0]]
]


partition[{x_}]:={{{x}}}
partition[{r__,x_}]:=partition[{r,x}]=Join@@(ReplaceList[#,{{b___,{S__},a___}:>{b,{S,x},a},{S__}:>{S,{x}}}]&/@partition[{r}])
mypartition[ls_,j_]:=Module[{allpt,tb1,tb2},allpt=partition[ls];
allpt=Permutations/@allpt//Flatten[#,1]&;tb1=Table[Join[{Join[{j},allpt[[ii,1]]]},allpt[[ii,2;;Length@allpt[[ii]]]]],{ii,Length@allpt}];
tb2=Join[{{j}},#]&/@allpt;
Join[tb1,tb2]]
mypartition[ls_]:=Module[{allpt,tb1,tb2},allpt=partition[ls];
allpt=Permutations/@allpt//Flatten[#,1]&;tb1=Table[Join[{Join[{1},allpt[[ii,1]]]},allpt[[ii,2;;Length@allpt[[ii]]]]],{ii,Length@allpt}]
(*tb2=Join[{{1}},#]&/@allpt;*)
(*Join[tb1,tb2]*)]
partitionPreNum[ls_]:=mypartition[ls[[2;;-1]],ls[[1]]]
partitionOrdered[ls_]:=Module[{allpt},allpt=partition[ls];
(*allpt=Permutations/@allpt//Flatten[#,1]&;*)Table[Join[{Join[{1},allpt[[ii,1]]]},allpt[[ii,2;;Length@allpt[[ii]]]]],{ii,Length@allpt}]]


partitionOd[fl_]:= Module[{allpt},allpt=partition[fl];Cases[allpt,x_/;OrderedQ[Flatten[x]]]]
MyPartitionOd[f_]:=Cases[partitionOd[Sort[f]]/.MapThread[Rule,{Sort[f],f}],x_/;And@@OrderedQ/@x]


inLeft[i_, ls_, top_]:=MemberQ[ls[[1;;Flatten[Position[ls,i]][[1]] ]],(top//Cases[Rule[ii_,i]:> ii])[[1]]]
inLeft[ls_, top_]:=And@@Table[inLeft[ls[[ii]],ls,top],{ii,2,Length@ls}]
ConsistentOrder[nm1_,top_]:=Module[{allorder},allorder=Join[{2},#]&/@Permutations[Range[3,nm1]];allorder=allorder//Cases[#,f_/;inLeft[f,top]]&]


list2top[ls_List]:= If[Length[ls]==1,{},Table[Rule[ls[[ii-1]],ls[[ii]]],{ii,2,Length@ls}]]
order2Mfun[od_,top_]:=Module[{ls,allowedPartitions},allowedPartitions=Cases[MyPartitionOd[od],f_/;And@@(MemberQ[top,#]&/@(list2top/@f//Flatten))]; ls=(M@@@allowedPartitions)/.M[f1_,f___]:> M[L@@Join[{1},f1],f];
Table[ls[[iii]]/.List[f1_,g___]:> R[(top//Cases[Rule[ii_,f1]:> ii])[[1]],f1,g],{iii,Length@ls}]]


Mt[ls_M]:=Module[{od=List@@ls,phat,leftv,num,den},phat=od/.L[f___]:>{}/.R[f1_,f___]:> f1//Flatten; phat=Join[{v},p/@phat];leftv=phat;
od=od/.L[f___]:>{f}/.R[f1_,f___]:> {f};
(*Print[leftv,od];*)
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],v]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=(*(-1)^Length@od*)dot[v,p[1]]Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den]
Mt[a_ ls_M]:=a Mt[ls]


NumPre[n_]:=Module[{alltop,lenTop,allods,res,res2},alltop=List[(KLTRM[n]//.nice/.CenterDot[Subscript[p, i_],Subscript[p, j_]]:> s[i,j]/.s[1,i_]:> 0/;i>2//Union//Expand//Total)/.a_ b_s:> b/;NumericQ[a]]/.Plus-> List//Flatten//Union;
alltop=(List@@@alltop)/.s-> Rule;
lenTop=Length/@((Keys@Select[Counts[#1],#==1&])&/@(Flatten/@(alltop/.Rule-> List)));
allods=Table[R[ConsistentOrder[n-1,alltop[[ii]]],alltop[[ii]]],{ii,Length@alltop}];
res=Table[allods[[ii]]/.R[f1_,f2_]:> (order2Mfun[#,f2]&/@f1)//Flatten,{ii,Length@alltop}];
(*Print[res,lenTop];*)
res=(Total/@res) . (Power[#,-1]&/@lenTop)/.M[f___]:> (-1)^Length[{f}] M[f];
res=res/.M[f___]:>  Mt[M[f]]
]


(*ClearAll[NumRec]*)
NumRec[{1}]:=sp[\[Epsilon][1]]
NumRec[{1,i2_}]:=-(1/(4dot[p[nv,1]]))dot[p[1],F[i2],lI[1]]sp[F[1],lI[1]]-1/dot[p[nv,1]] dot[p[nv],F[1],F[i2],lI[1]]sp[lI[1]]
NumRec[od_List/;Length[od]>2]:=NumRec[od]=Module[{allsets,leftall,rightall,leftv,res,remains,partitions},
partitions=KSetPartitions[od,2];
(*Print[partitions];*)
allsets=Drop[partitions,1];
remains=partitions[[1]];
res=Sum[
leftall=allsets[[ii,1]];
rightall=allsets[[ii,2]];
leftv=p@@Join[{nv},Select[od,#<rightall[[1]]&]];
NumRec[leftall]/.sp[f__]:>(-1)^(Length[rightall]+1) (dot[leftv[[2;;-1]],F[rightall],lI[rightall[[1]]]]/dot[p@@Join[{nv},leftall]])sp[f,p@@Join[{nv},leftall],lI[rightall[[1]]]],{ii,Length@allsets}];
res=res+(-1)^Length[remains[[2]]]/(4dot[p[nv,1]]) dot[p[1],F[remains[[2]]],lI[1]]sp[F[1],lI[1]]+(-1)^Length[remains[[2]]]/dot[p[nv,1]] dot[p[nv],F[1],F[remains[[2]]],lI[1]]sp[lI[1]]
]
numRecQCD[n_]:=NumRec[Range[n-2]]/.sp[f__]:>spAB[p[n],f,p[n-1]]/.nv->n/.F[f_List]:>T@@(F/@f)/.T[f__]:>f
numRecQCD[od_List]:=NumRec[Range[Length@od]]/.sp[f__]:>spAB[p[Length@od],f,p[-1+Length@od]]/.nv->Length@od/.F[f_List]:>T@@(F/@f)/.T[f__]:>f


(*Programs on the Hopf algebra*)


(* ::Subsection::Closed:: *)
(*Hopf algebra*)


mC[f_]:=If[Head[f]=!=Plus&&Length[f]>=2,True,False]
nMC[f_]:= Not[mC[f]]&& Head[f]=!=Plus


ExpandNCM[(h:NonCommutativeMultiply)[a___,b_Plus,c___]]:=Distribute[h[a,b,c],Plus,h,Plus,ExpandNCM[h[##]]&]
ExpandNCM[(h:NonCommutativeMultiply)[a___,b_Times,c___]]:=Most[b]ExpandNCM[h[a,Last[b],c]]
ExpandNCM[a_Plus]:=ExpandNCM[#]&/@a
ExpandNCM[c_Times]:=c[[1]]*ExpandNCM[c[[2]]]
ExpandNCM[a_]:=ExpandAll[a]
ExpandT[f_]:=((NonCommutativeMultiply@@f//ExpandNCM)/.T[hh__]:> hh/.NonCommutativeMultiply-> T)
ExpandCT[f_]:=((NonCommutativeMultiply@@f//ExpandNCM)/.NonCommutativeMultiply-> CircleTimes)
GetTId[f_]:=(f/.T[hh__]:> hh)
ExpandCircleTimes[f_]:=(f/.CircleTimes->NonCommutativeMultiply//ExpandNCM)/.NonCommutativeMultiply->CircleTimes


NC[i_,j_]:=NonCommutativeMultiply[i,j]-NonCommutativeMultiply[j,i]


\[FivePointedStar][f_T?mC,g_T?mC]:=(T[f[[1]],\[FivePointedStar][f[[2;;-1]],g]]//ExpandT)+(T[g[[1]],\[FivePointedStar][f,g[[2;;-1]]]]//ExpandT)-(T[Join[f[[1]],g[[1]]],\[FivePointedStar][f[[2;;-1]],g[[2;;-1]]]]//ExpandT)/.T[gg___]:>T@@(Sort/@{gg})
\[FivePointedStar][f_T?nMC,g_T?mC]:=T[f[[1]],(g[[1;;-1]]//GetTId)]+(T[g[[1]],\[FivePointedStar][f,g[[2;;-1]]]]//ExpandT)-T[Join[f[[1]],g[[1]]],(g[[2;;-1]]//GetTId)]/.T[gg___]:>T@@(Sort/@{gg})
\[FivePointedStar][f_T?mC,g_T?nMC]:=(T[f[[1]],\[FivePointedStar][f[[2;;-1]],g]]//ExpandT)+T[g[[1]],(f[[1;;-1]]//GetTId)]-T[Join[f[[1]],g[[1]]],(f[[2;;-1]]/.T[hh__]:> hh)]/.T[gg___]:>T@@(Sort/@{gg})
\[FivePointedStar][f_T?nMC,g_T?nMC]:=T[f[[1]],g[[1]]]+T[g[[1]],f[[1]]]-T[Join[f[[1]],g[[1]]]]/.T[gg___]:>T@@(Sort/@{gg})
\[FivePointedStar][f_Plus,g_]:=(NonCommutativeMultiply[f,g]//ExpandNCM)/.NonCommutativeMultiply-> \[FivePointedStar]
\[FivePointedStar][f_,g_Plus]:=(NonCommutativeMultiply[f,g]//ExpandNCM)/.NonCommutativeMultiply-> \[FivePointedStar]
\[FivePointedStar][f_Times,g_]:=(NonCommutativeMultiply[f,g]//ExpandNCM)/.NonCommutativeMultiply-> \[FivePointedStar]
\[FivePointedStar][f_,g_Times]:=(NonCommutativeMultiply[f,g]//ExpandNCM)/.NonCommutativeMultiply-> \[FivePointedStar]
\[FivePointedStar][\[DoubleStruckCapitalI],g_]:=g
\[FivePointedStar][g_,\[DoubleStruckCapitalI]]:=g
\[FivePointedStar][f1_,f2_,g__]:=\[FivePointedStar][\[FivePointedStar][f1,f2],g]
\[FivePointedStar][g_T]:=g
\[FivePointedStar][f_GT,g_GT]:=If[f[[2]]=={}||g[[2]]=={},T@@Join[f[[2]],g[[2]]]/.T[h___]:> GT[Join[f[[1]],g[[1]]],SortBy[#1,Position[Join[f[[1]],g[[1]]],#]&]&/@{h}],\[FivePointedStar][T@@f[[2]],T@@g[[2]]]/.T[h__]:> GT[Join[f[[1]],g[[1]]],SortBy[#1,Position[Join[f[[1]],g[[1]]],#]&]&/@{h}]
]
\[FivePointedStar][f_ET ,g_ET]:=Module[{groupGens},groupGens=Join[f[[2]],g[[2]]] ;\[FivePointedStar][f[[1]] ,g[[1]]]/.GT[gg__]:>ET[GT[gg],groupGens]]
\[FivePointedStar][J1_FT,J2_FT]:=Module[{tausJ1,tausJ2,ranks0a,ranks0,ranks2,tlist,res0,res2,res},
tausJ1=T@@J1[[All,1]];
tausJ2=T@@J2[[All,1]];
ranks0a=List@@J1[[All,2]];
ranks0=Join[List@@J1[[All,2]],{0}];
ranks2=Join[List@@J1[[All,2]],{2}];
tlist=\[FivePointedStar][tausJ1,tausJ2];
res0=tlist/.T[f__]:>FT@@MapThread[List,{{f},ranks0}]/;Length[{f}]==Length[ranks0]/.T[f__]:>FT@@MapThread[List,{{f},ranks0a}]/;Length[{f}]==(Length[ranks0]-1);
res2=tlist/.T[f__]:>FT@@MapThread[List,{{f},ranks2}]/;Length[{f}]==Length[ranks2]/.T[f__]:>0/;Length[{f}]==(Length[ranks0]-1);
res0+res2
]





\[ScriptCapitalK][i_,0]:=ET[GT[{i},{}],{\[DoubleStruckA][i]}]
\[ScriptCapitalK][i_,1]:=ET[GT[{i},{{i}}],{}]


\[CapitalOmega][f_,g_]:=\[FivePointedStar][f,g]-\[FivePointedStar][g,f]
\[CapitalOmega][f1_,f2_,g__]:=\[CapitalOmega][\[CapitalOmega][f1,f2],g]


\[DoubleStruckCapitalS][f_T]:=If[Length[f]>1,-f-Sum[\[FivePointedStar][\[DoubleStruckCapitalS][f[[1;;i]]],f[[i+1;;Length[f]]]],{i,1,Length[f]-1}],-f]/.T[gg___]:>T@@(Sort/@{gg})
\[DoubleStruckCapitalS][f_Plus]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> \[DoubleStruckCapitalS]
\[DoubleStruckCapitalS][f_Times]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> \[DoubleStruckCapitalS]
\[DoubleStruckCapitalS][\[DoubleStruckCapitalI]]:=\[DoubleStruckCapitalI]


\[CapitalDelta][f_T]:=If[Length[f]>1,CircleTimes[\[DoubleStruckCapitalI],f]+CircleTimes[f,\[DoubleStruckCapitalI]]+Sum[CircleTimes[f[[1;;i]],f[[i+1;;Length[f]]]]//ExpandCT,{i,1,Length[f]-1}],CircleTimes[\[DoubleStruckCapitalI],f]+CircleTimes[f,\[DoubleStruckCapitalI]]]
\[CapitalDelta][f_Plus]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> \[CapitalDelta]
\[CapitalDelta][f_Times]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> \[CapitalDelta]
\[CapitalDelta][f_GT]:=\[CapitalDelta][T@@f[[2]]]/.T[h__]:>GT[f[[1]],{h}]
\[CapitalDelta][\[DoubleStruckCapitalI]]:=CircleTimes[\[DoubleStruckCapitalI],\[DoubleStruckCapitalI]]


picCTR[ls_]:=CircleTimes[f__,g_]:> 0/;g===\[DoubleStruckCapitalI]||(Sort[Flatten[List@@g]]=!=ls)
picCTL[ls_]:=CircleTimes[g_,f__]:> 0/;g===\[DoubleStruckCapitalI]||(Sort[Flatten[List@@g]]=!=ls)
picTR[h_T]:=CircleTimes[f__,g_]:> 0/;g=!=h
picTL[h_T]:=CircleTimes[g_,f__]:> 0/;g=!=h


CtCut[f__]:=Module[{od={f},phat,leftv,num,den},phat=Range[Length@od];
If[MemberQ[od,\[DoubleStruckCapitalI]],Return[0]];
phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],#<First[Min[od[[i]]]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],p[od[[i,1]]]],{i,2,Length@od}];
den=Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]
CTCut[f__]:= (CtCut@@(({f}/.T[g__]:>Flatten[{g}])))Times[f]


GT2Cut[fc_,fk_]:=Module[{od=fk,odc=fc,phat,leftv,num,den},phat=Range[Length@od];
phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],Position[odc,#][[1,1]]<Position[odc,MinimalBy[od[[i]],Position[odc,#1]&]//First][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],p[MinimalBy[od[[i]],Position[odc,#1]&]//First]],{i,2,Length@od}];
den=Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]
CT2GT[f__]:= If[MemberQ[{f},\[DoubleStruckCapitalI]],Return[0],(GT2Cut@@{({f}/.GT[c_,g_]:>c)//Union//Flatten,({f}/.GT[c_,g_]:>Flatten[g])})Times[f]]


(* ::Subsection::Closed:: *)
(*basic functions*)


$globalDim = 4;


setGlobalDim[dim_] /; If[!NumericQ[dim] && protectedQ[dim],Message[setGlobalDim::ssle]; False, True] := Set[$globalDim,dim];


protectedQ[expr_Symbol] := MemberQ[Attributes[expr],Protected]
(* protectedQ[expr_] := TrueQ[MemberQ[Quiet[Attributes[expr]], Protected] || !FreeQ[expr, _?(MemberQ[Quiet[Attributes[#]], Protected] &), Infinity]] *)
protectedQ[_] := False


patternFreeQ[expr_] := FreeQ[expr,Pattern]&&FreeQ[expr,Blank]&&FreeQ[expr,BlankSequence]&&FreeQ[expr,BlankNullSequence]


declareDistributive[function_Symbol,test_,default_:True] := (
	declareDistributiveMiddle[function,test,default];
	declareDistributiveFirst[function,test,default];
	declareDistributiveLast[function,test,default];
)


declareDistributiveMiddle[function_Symbol,test_,default_:True] := (
	If[!MemberQ[Options[function][[All,1]],distributive],AppendTo[Options[function],distributive->default]];
	function[__,0,__]  := 0;
	function[p1__,n_ p2_?test,p3__,OptionsPattern[]] := n*function[p1, p2, p3]  /; OptionValue[distributive];
	(*function[p1__,p2_Plus,p3__,OptionsPattern[]] := Total[function[p1, #, p3] & /@ List @@ p2]  /; OptionValue[distributive];*)
	(*function[p1__,Times[a___, p2_Plus, b___],p3__,OptionsPattern[]] := function[p1,Expand[a*p2*b],p3]  /; OptionValue[distributive];*)
	function[p1__, a_. p2_Plus ,p3__,OptionsPattern[]] := function[p1, #, p3]& /@ Distribute[a p2]  /; OptionValue[distributive];
)


declareDistributiveFirst[function_Symbol,test_,default_:True] := (
	If[!MemberQ[Options[function][[All,1]],distributive],AppendTo[Options[function],distributive->default]];
	function[0,___] := 0;
	function[n_ p1_?test,p2___,OptionsPattern[]] := n*function[p1, p2]  /; OptionValue[distributive];
	(*function[p1_Plus,p2___,OptionsPattern[]] := Total[function[#,p2] & /@ List @@ p1]  /; OptionValue[distributive];*)
	(*function[Times[a___, p1_Plus, b___],p2___,OptionsPattern[]] := function[Expand[a*p1*b],p2]  /; OptionValue[distributive];*)
	function[a_. p1_Plus ,p2___,OptionsPattern[]] := function[#, p2]& /@ Distribute[a p1]  /; OptionValue[distributive];
)


declareDistributiveLast[function_Symbol,test_,default_:True] := (
	If[!MemberQ[Options[function][[All,1]],distributive],AppendTo[Options[function],distributive->default]];
	function[___,0] := 0;
	function[p1___,n_ p2_?test,OptionsPattern[]] := n*function[p1, p2]  /; OptionValue[distributive];
	(*function[p1___,p2_Plus,OptionsPattern[]] := Total[function[p1,#] & /@ List @@ p2]  /; OptionValue[distributive];*)
	(*function[p1___,Times[a___, p2_Plus, b___],OptionsPattern[]] := function[p1,Expand[a*p2*b]]  /; OptionValue[distributive];*)
	function[p1___, a_. p2_Plus, OptionsPattern[]] := function[p1, #]& /@ Distribute[a p2]  /; OptionValue[distributive];
)


distributiveRules[function_Symbol,test_,excfront_:0,excback_:0] := With[
	{uniqueFront = Sequence @@ (Unique[]& /@ Range[excfront]), uniqueBack = Sequence @@ (Unique[]& /@ Range[excback])},
	{
		HoldPattern[function][Sequence @@ (Pattern[#, _]& /@ {uniqueFront}), ___, 0, ___, Sequence @@ (Pattern[#, _]& /@ {uniqueBack})] :> 0,
		HoldPattern[function][Sequence @@ (Pattern[#, _]& /@ {uniqueFront}), p1___, n_ p2_?test, p3___, Sequence @@ (Pattern[#, _]& /@ {uniqueBack})] :> n*function[uniqueFront, p1, p2, p3, uniqueBack],
		HoldPattern[function][Sequence @@ (Pattern[#, _]& /@ {uniqueFront}), p1___, a_. p2_Plus, p3___, Sequence @@ (Pattern[#, _]& /@ {uniqueBack})] :> (function[uniqueFront, p1, #, p3, uniqueBack]& /@ Distribute[a p2])
	}
]


declareSymmetric[function_Symbol] := SetAttributes[function,Orderless];


(* Todo: give an option to not automatically antisymmetrise *)
declareAntisymmetric[function_Symbol] := (
	(*orderes arguments and puts correct sign if it's not used as a pattern!*)
	function[x___] := Signature[{x}] (function @@ Sort@{x}) /; (! OrderedQ[{x}] && FreeQ[{x}, Pattern]);
	(*two equal arguments make the function vanish*)
	function[___,x_,x_,___] := 0;
)


(* ::Subsection::Closed:: *)
(*Vectors*)


(* Todo: fix dimensionality of dot and eps *)
(* $declaredTensors = {{dot,$globalDim,2},{eps,4,4}}; *)
$declaredTensors = {{eta,4,2}};
$declaredTensorHeads = {{spOuter,4,1},{outer,4,2}};
$declaredAntiTensorHeads={};
(* Todo: there should be an (optional) function that declares these...? *)
$extMomLabel = "p";
$loopLabel = "l";


antiTensorQ[tens_]:=tensorQ[tens]&&MemberQ[$declaredAntiTensorHeads,Head[tens]]


tensorQ[tens_] := If[patternFreeQ[tens],MemberQ[$declaredTensors[[All,1]],tens] || MemberQ[$declaredTensorHeads[[All,1]],Head[tens]],False]


vectorQ[vec_] := If[patternFreeQ[vec],tensorQ[vec] && (tensorRank[vec]===1),False]


tensorQR2[ten_] := If[patternFreeQ[ten],tensorQ[ten] && (tensorRank[ten]===2),False]


tensorDim[tens_] /; If[tensorQ[tens],True,Message[tensorDim::argx]; False] :=
	FirstCase[$declaredTensors,_?(#[[1]]===tens&),FirstCase[$declaredTensorHeads,_?(#[[1]]===Head[tens]&)]][[2]]


tensorRank[tens_] /; If[tensorQ[tens],True,Message[tensorRank::argx]; False] :=
	FirstCase[$declaredTensors,_?(#[[1]]===tens&),FirstCase[$declaredTensorHeads,_?(#[[1]]===Head[tens]&)]][[3]]


expandTensor[f_dot]:=Module[{lt=Length[f],ls,tlsid,lsT},tlsid=Range[lt];ls=List@@f;lsT=dot@@MapThread[T,{ls,tlsid}];(lsT//.dot[ff__,T[g_?tensorQR2,i_],h__]:> dot[ff,lI[i,1]]dressIndex[g,{lI[i,1],lI[i,2]}]dot[lI[i,2],h])/.T[gg_,i_]:> gg]


(* Todo: should I really distinguish between declareVector and declareTensor? *)
Options[declareVector] := {"dim"->$globalDim,"verbose"->True}
declareVector[expr_List,options___] := Scan[declareVector[#,options]&,expr];
declareVector[vec_Symbol,OptionsPattern[]] /; If[!protectedQ[vec],True,Message[declareVector::safe]; False] := (
	If[tensorQ[vec],undeclareTensor[vec];];
	vec /: MakeBoxes[vec[KiHA`lI[KiHA`Private`i_]],TraditionalForm] :=
		FormBox[SuperscriptBox[MakeBoxes[vec,TraditionalForm],SubscriptBox["\[Nu]",ToString[KiHA`Private`i]]],TraditionalForm];
	AppendTo[$declaredTensors,{vec,OptionValue["dim"],1}];
	If[OptionValue["verbose"],Print[ToString[vec]<>" declared as a "<>ToString[OptionValue["dim"]]<>"-dimensional vector."]];
)
declareVector[__] := $Failed


Options[declareVectorHead] := {"dim"->$globalDim,"verbose"->True}
declareVectorHead[expr_List,options___] := Scan[declareVectorHead[#,options]&,expr];
declareVectorHead[vec_Symbol,OptionsPattern[]] /; If[!protectedQ[vec],True,Message[declareVectorHead::safe]; False] := (
	If[tensorQ[vec[1]],undeclareTensorHead[vec];];
	vec /: MakeBoxes[vec[i_],TraditionalForm] := FormBox[SubscriptBox[ToString[vec],MakeBoxes[i,TraditionalForm]],TraditionalForm];
	MakeBoxes[vec[KiHA`Private`i_][basicfun`lI[KiHA`Private`j_]],TraditionalForm] :=
		FormBox[SuperscriptBox[MakeBoxes[vec[KiHA`Private`i],TraditionalForm],SubscriptBox["\[Nu]",ToString[KiHA`Private`j]]],TraditionalForm];
	AppendTo[$declaredTensorHeads,{vec,OptionValue["dim"],1}];
	If[OptionValue["verbose"],Print[ToString[vec]<>" declared as a "<>ToString[OptionValue["dim"]]<>"-dimensional vector header."]];
)
declareVectorHead[__] := $Failed


Options[declareTensor] := {"dim"->$globalDim,"rank"->1,"verbose"->True}
declareTensor[expr_List,options___] := Scan[declareTensor[#,options]&,expr];
declareTensor[tens_Symbol,OptionsPattern[]] /; If[!protectedQ[tens],True,Message[declareTensor::safe]; False] := (
	If[tensorQ[tens],undeclareTensor[tens];];
	(* Todo: formatting for free indices on tensors *)
	AppendTo[$declaredTensors,{tens,OptionValue["dim"],OptionValue["rank"]}];
	If[OptionValue["verbose"],Print[ToString[tens]<>" declared as a "<>ToString[OptionValue["dim"]]<>"-dimensional rank-"<>ToString[OptionValue["rank"]]<>" tensor."]];
)
declareTensor[__] := $Failed


Options[declareTensorHead] := {"dim"->$globalDim,"rank"->1,"verbose"->True}
declareTensorHead[expr_List,options___] := Scan[declareTensorHead[#,options]&,expr];
declareTensorHead[tens_Symbol,OptionsPattern[]] /; If[!protectedQ[tens],True,Message[declareTensorHead::safe]; False] := (
	If[tensorQ[tens],undeclareTensorHead[tens];];
	tens /: MakeBoxes[tens[i_],TraditionalForm] := FormBox[SubscriptBox[ToString[tens],MakeBoxes[i,TraditionalForm]],TraditionalForm];
	AppendTo[$declaredTensorHeads,{tens,OptionValue["dim"],OptionValue["rank"]}];
	If[OptionValue["verbose"],Print[ToString[tens]<>" declared as a "<>ToString[OptionValue["dim"]]<>"-dimensional rank-"<>ToString[OptionValue["rank"]]<>" tensor header."]];
)
declareTensorHead[__] := $Failed


declareAntiTensorHead[expr_List] := Scan[declareAntiTensorHead[#]&,expr];
declareAntiTensorHead[tens_Symbol]  := (
	AppendTo[$declaredAntiTensorHeads,tens];
)


Options[undeclareTensor] = {"verbose"->True};
undeclareTensor[args_List,options___] := Scan[undeclareTensor[#,options]&,args];
undeclareTensor[vec_?tensorQ,OptionsPattern[]] /; !MemberQ[$declaredTensorHeads,Head[vec]] := (
	$declaredTensors = DeleteCases[$declaredTensors,_?(#[[1]]==vec &)];
(* 	If[Head[vec]===Symbol,
		FormatValues[vec] = {}; DownValues[vec] = {}; UpValues[vec] = {}; SubValues[vec] = {};,
		FormatValues[Evaluate[Head[vec]]] = DeleteCases[FormatValues[Evaluate[Head[vec]]],_?(MemberQ[#,vec,Infinity]&)];
		UpValues[Evaluate[Head[vec]]] = DeleteCases[UpValues[Evaluate[Head[vec]]],_?(MemberQ[#,vec,Infinity]&)];
		DownValues[Evaluate[Head[vec]]] = DeleteCases[DownValues[Evaluate[Head[vec]]],_?(MemberQ[#,vec,Infinity]&)];
		SubValues[Evaluate[Head[vec]]] = DeleteCases[SubValues[Evaluate[Head[vec]]],_?(MemberQ[#,vec,Infinity]&)];
	]; *)
	If[OptionValue["verbose"],Print[ToString[vec]<>" undeclared as a tensor."]];
)
undeclareTensor[__] := Null


Options[undeclareTensorHead] = {"verbose"->True};
undeclareTensorHead[args_List,options___] := Scan[undeclareTensorHead[#,options]&,args];
undeclareTensorHead[vec_Symbol,OptionsPattern[]] := (
	$declaredTensorHeads = DeleteCases[$declaredTensorHeads,_?(#[[1]]==vec &)];
	If[OptionValue["verbose"],Print[ToString[vec]<>" undeclared as a tensor header."]];
)
undeclareTensorHead[__] := Null


SetAttributes[dot,{NHoldAll}];
dot[expr_] := dot[expr,expr] /; patternFreeQ[expr]
dot[p_?tensorQ,__] /; If[vectorQ[p],False,Message[dot::badrank]; True] := $Failed
dot[__,p_?tensorQ] /; If[vectorQ[p],False,Message[dot::badrank]; True] := $Failed
dot[__,p_?tensorQ,__] /; If[tensorRank[p]===2,False,Message[dot::badrank]; True] := $Failed
declareDistributive[dot,tensorQ];
dot[p1_,p2_] := dot[p2,p1] /; (!OrderedQ[{p1,p2}] && patternFreeQ[{p1,p2}])

dot[vec_?vectorQ,ix_lI] := vec[ix]
dot[ix_lI,vec_?vectorQ] := vec[ix]
dot[ix1_lI,ix2_lI] := eta[ix1,ix2]

dot[vec_?vectorQ,ix_ui] := vec[ix]
dot[ix_ui,vec_?vectorQ] := vec[ix]

dot[vec_?vectorQ,ix_di] := vec[ix]
dot[ix_di,vec_?vectorQ] := vec[ix]

dot /: MakeBoxes[dot[p1_Plus,p1_Plus],TraditionalForm] := FormBox[SuperscriptBox[RowBox[{"(",MakeBoxes[p1,TraditionalForm],")"}],"2"],TraditionalForm]
dot /: MakeBoxes[dot[p1_,p1_],TraditionalForm] := FormBox[SuperscriptBox[MakeBoxes[p1,TraditionalForm],"2"],TraditionalForm]
dot /: MakeBoxes[dot[p1_Plus,p2_Plus],TraditionalForm] := FormBox[RowBox[{"(",MakeBoxes[p1,TraditionalForm],")","\[CenterDot]","(",MakeBoxes[p2,TraditionalForm],")"}],TraditionalForm]
dot /: MakeBoxes[dot[p1_Plus,p2_],TraditionalForm] := FormBox[RowBox[{"(",MakeBoxes[p1,TraditionalForm],")","\[CenterDot]",MakeBoxes[p2,TraditionalForm]}],TraditionalForm]
dot /: MakeBoxes[dot[p1_,p2_Plus],TraditionalForm] := FormBox[RowBox[{MakeBoxes[p1,TraditionalForm],"\[CenterDot]","(",MakeBoxes[p2,TraditionalForm],")"}],TraditionalForm]
dot /: MakeBoxes[dot[p1_,p2_],TraditionalForm] := FormBox[RowBox[{MakeBoxes[p1,TraditionalForm],"\[CenterDot]",MakeBoxes[p2,TraditionalForm]}],TraditionalForm]

(* Todo: fix for arguments as sums *)
dot /: MakeBoxes[dot[args__],TraditionalForm] := FormBox[RowBox[Join[Riffle[MakeBoxes[#,TraditionalForm]&/@{args},"\[CenterDot]"]]], TraditionalForm]


SetAttributes[tr, {NHoldAll}];
declareDistributive[tr, tensorQ];
tr[] := $globalDim
tr[args__] /; patternFreeQ[{args}] && {args}[[1]]=!=Sort[{args}][[1]] := tr@@RotateLeft@{args}


declareAntisymmetric[eps];
SetAttributes[eps,NHoldAll];
declareDistributive[eps,tensorQ];

eps /: MakeBoxes[eps[expr__],TraditionalForm] := FormBox[RowBox[Join[{"\[Epsilon]("}, Riffle[MakeBoxes[#,TraditionalForm]&/@{expr}, ","],{")"}]], TraditionalForm]


declareDistributive[outer,vectorQ,tensorQR2];

outer[args___] /; Length[{args}]=!=2 := (Message[outer::argx,Length[{args}]]; $Failed)
outer[v1_,v2_] /; (!vectorQ[v1] || !vectorQ[v2]) && patternFreeQ[{v1,v2}] := (Message[outer::notvectors]; $Failed)

outer[v1_?vectorQ,v2_?vectorQ][ix1_lI,ix2_lI] := v1[ix1]*v2[ix2]
outer /: dot[p1__,outer[v1_,v2_],p2__] := dot[p1,v1]*dot[v2,p2]
outer /: tr[p1___,outer[v1_,v2_],p2___] := dot[v2,p2,p1,v1]

outer /: MakeBoxes[outer[p1_,p2_],TraditionalForm] := FormBox[RowBox[{MakeBoxes[p1, TraditionalForm], "\[TensorProduct]", MakeBoxes[p2, TraditionalForm]}], TraditionalForm]


SetAttributes[eta,{Orderless}];
eta /: dot[p1__,eta,p2__] := dot[p1,p2]
eta /: tr[p1___,eta,p2___] := tr[p1,p2]

eta /: MakeBoxes[eta,TraditionalForm] := "\[Eta]"


SetAttributes[mu,{Orderless,NHoldAll}];
mu[expr_] := mu[expr,expr] /; patternFreeQ[{expr}]
mu[vec_?vectorQ,__] /; tensorDim[vec]==4 := 0
mu[expr__] := 0 /; (patternFreeQ[{expr}] && $globalDim==4)
mu[expr__,OptionsPattern[]] /; (patternFreeQ[{expr}] && Length[{expr}]>2) := (Message[mu::argx,Length[{expr}]]; $Failed)
mu[OptionsPattern[]] := (Message[mu::argx,0]; $Failed)
declareDistributive[mu,vectorQ];

mu /: MakeBoxes[mu[p1_,p2_],TraditionalForm] :=
	If[AllTrue[{p1,p2},ToString[Head[#]]==$loopLabel &],
		FormBox[SubscriptBox["\[Mu]",ToString[p1[[1]]]<>ToString[p2[[1]]]], TraditionalForm],
		FormBox[RowBox[{"\[Mu](",MakeBoxes[p1,TraditionalForm],",",MakeBoxes[p2,TraditionalForm],")"}],TraditionalForm]
	]
	
	
declareAntisymmetric[aMu];
SetAttributes[aMu,NHoldAll];
declareDistributive[aMu,tensorQ];


(* ::Subsection::Closed:: *)
(*Spinors*)


spQ[expr_] := MemberQ[{spA,spB},Head[expr]]


SetAttributes[{spA,spB},NHoldAll];
(* spA[p_?tensorQ] /; If[tensorDim[p]===4 && onShellMass[p]===0,False,Message[spA::offshell]; True] := $Failed *)

(* Analytic continuation in 4D *)
(* Todo: remove this as a requirement... *)
(* TODO: which continuation is better? there's a problem with Zvi's *)
(*spA[-p_?vectorQ] := I*spA[p]
DownValues[spB] = DownValues[spA] /. {HoldPattern[spA] :> spB};*)
spA[-p_?vectorQ] := - spA[p]
spB[-p_?vectorQ] := spB[p]

spA /: MakeBoxes[spA[vec_],TraditionalForm] := If[ToString[Head[vec]]===$extMomLabel,
		FormBox[SubscriptBox["\[Lambda]",ToString[vec[[1]]]],TraditionalForm],
		FormBox[SubscriptBox["\[Lambda]",MakeBoxes[vec,TraditionalForm]],TraditionalForm]
]
spB /: MakeBoxes[spB[vec_],TraditionalForm] := If[ToString[Head[vec]]===$extMomLabel,
		FormBox[SubscriptBox[MakeBoxes[OverTilde["\[Lambda]"],TraditionalForm],ToString[vec[[1]]]],TraditionalForm],
		FormBox[SubscriptBox[MakeBoxes[OverTilde["\[Lambda]"],TraditionalForm],MakeBoxes[vec,TraditionalForm]],TraditionalForm]
]


spAA[p1_,p2___,p3_] := spp[spA[p1],p2,spA[p3]]
spBB[p1_,p2___,p3_] := spp[spB[p1],p2,spB[p3]]
spAB[p1_,p2___,p3_] := spp[spA[p1],p2,spB[p3]]
spBA[p1_,p2___,p3_] := spp[spB[p1],p2,spA[p3]]


SetAttributes[spp,NHoldAll];

spp[sp_?spQ,sp_?spQ] := 0 /; patternFreeQ[{sp}]
spp[sp1_spA,sp2_spA] := -spp[sp2,sp1] /; (!OrderedQ[{sp1,sp2}] && patternFreeQ[{sp1,sp2}])
spp[sp1_spB,sp2_spB] := -spp[sp2,sp1] /; (!OrderedQ[{sp1,sp2}] && patternFreeQ[{sp1,sp2}])
spp[spA[p1_],args___,spA[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && OddQ[GLength[{args}]] := 0
spp[spB[p1_],args___,spB[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && OddQ[GLength[{args}]] := 0
spp[spA[p1_],args___,spB[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && EvenQ[GLength[{args}]] := 0
spp[spB[p1_],args___,spA[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && EvenQ[GLength[{args}]] := 0

(*spp[___,p_?tensorQ,___] /; If[tensorRank[p]===1,False,Message[spp::badrank]; True] := $Failed*)

spp[spA[p_],p_,___] := 0
spp[spB[p_],p_,___] := 0
spp[___,p_,spA[p_]] := 0
spp[___,p_,spB[p_]] := 0
spp[p1___,p2_,p2_,p3___] /; patternFreeQ[p1,p2,p3] := dot[p2]*spp[p1,p3] 

spp[] := 1

declareDistributiveFirst[spp,(tensorQ[#] || spQ[#] &)];
declareDistributiveLast[spp,(tensorQ[#] || spQ[#] &)];
declareDistributiveMiddle[spp,tensorQ];

spp /: MakeBoxes[spp[spA[p1_],spA[p2_]],TraditionalForm] :=
	If[AllTrue[{p1,p2},ToString[Head[#]]==$extMomLabel &],
		"\[LeftAngleBracket]"<>ToString[p1[[1]]]<>ToString[p2[[1]]]<>"\[RightAngleBracket]",
		FormBox[RowBox[{"\[LeftAngleBracket]",MakeBoxes[p1,TraditionalForm],MakeBoxes[p2,TraditionalForm],"\[RightAngleBracket]"}], TraditionalForm]
	]
spp /: MakeBoxes[spp[spB[p1_],spB[p2_]],TraditionalForm] :=
	If[AllTrue[{p1,p2},ToString[Head[#]]==$extMomLabel &],
		"["<>ToString[p1[[1]]]<>ToString[p2[[1]]]<>"]",
		FormBox[RowBox[{"[",MakeBoxes[p1,TraditionalForm],MakeBoxes[p2,TraditionalForm],"]"}], TraditionalForm]
	]
spp /: MakeBoxes[spp[spA[p1_],p2__,spA[p3_]],TraditionalForm] :=
	If[AllTrue[{p1,p3},ToString[Head[#]]==$extMomLabel &],
		FormBox[RowBox[Join[{"\[LeftAngleBracket]"<>ToString[p1[[1]]]<>"|"},Riffle[MakeBoxes[#,TraditionalForm]&/@{p2},"|"],{"|"<>ToString[p3[[1]]]<>"\[RightAngleBracket]"}]],TraditionalForm],
		FormBox[RowBox[Join[{"\[LeftAngleBracket]"}, Riffle[MakeBoxes[#,TraditionalForm]&/@{p1,p2,p3},"|"],{"\[RightAngleBracket]"}]], TraditionalForm]
	]
spp /: MakeBoxes[spp[spB[p1_],p2__,spB[p3_]],TraditionalForm] :=
	If[AllTrue[{p1,p3},ToString[Head[#]]==$extMomLabel &],
		FormBox[RowBox[Join[{"["<>ToString[p1[[1]]]<>"|"},Riffle[MakeBoxes[#,TraditionalForm]&/@{p2},"|"],{"|"<>ToString[p3[[1]]]<>"]"}]],TraditionalForm],
		FormBox[RowBox[Join[{"["}, Riffle[MakeBoxes[#,TraditionalForm]&/@{p1,p2,p3},"|"],{"]"}]], TraditionalForm]
	]
spp /: MakeBoxes[spp[spA[p1_],p2__,spB[p3_]],TraditionalForm] :=
	If[AllTrue[{p1,p3},ToString[Head[#]]==$extMomLabel &],
		FormBox[RowBox[Join[{"\[LeftAngleBracket]"<>ToString[p1[[1]]]<>"|"},Riffle[MakeBoxes[#,TraditionalForm]&/@{p2},"|"],{"|"<>ToString[p3[[1]]]<>"]"}]],TraditionalForm],
		FormBox[RowBox[Join[{"\[LeftAngleBracket]"}, Riffle[MakeBoxes[#,TraditionalForm]&/@{p1,p2,p3},"|"],{"]"}]], TraditionalForm]
	]
spp /: MakeBoxes[spp[spB[p1_],p2__,spA[p3_]],TraditionalForm] :=
	If[AllTrue[{p1,p3},ToString[Head[#]]==$extMomLabel &],
		FormBox[RowBox[Join[{"["<>ToString[p1[[1]]]<>"|"},Riffle[MakeBoxes[#,TraditionalForm]&/@{p2},"|"],{"|"<>ToString[p3[[1]]]<>"\[RightAngleBracket]"}]],TraditionalForm],
		FormBox[RowBox[Join[{"["}, Riffle[MakeBoxes[#,TraditionalForm]&/@{p1,p2,p3},"|"],{"\[RightAngleBracket]"}]], TraditionalForm]
	] 


GLength[f_List]:=Sum[If[Head[f[[ii]]]===lI,1,tensorRank[(f[[ii]]//Cases[#,x_/;tensorQ[x],{0,\[Infinity]}]&)//First]],{ii,Length@f}]


declareDistributive[spOuter,spQ[#] || vectorQ[#] &]

spOuter[vec_?vectorQ] := spOuter[spA[vec],spB[vec]]
spOuter[expr_] /; patternFreeQ[{expr}] := (Message[spOuter::notvector]; $Failed)
spOuter[args___] /; Length[{args}]>2 && patternFreeQ[{args}] := (Message[spOuter::argx,Length[{args}]]; $Failed)
spOuter[] := (Message[spOuter::argx,0]; $Failed)

spOuter[sp1_,sp2_] /; (!spQ[sp1] || !spQ[sp2]) && patternFreeQ[{sp1,sp2}] := (Message[spOuter::badspinors]; $Failed)
spOuter[sp1_spB,sp2_spA] := spOuter[sp2,sp1]
spOuter[sp1_spA,sp2_spA] := (Message[spOuter::badspinors]; $Failed)
spOuter[sp1_spB,sp2_spB] := (Message[spOuter::badspinors]; $Failed)

spOuter[sp1_?spQ,sp2_?spQ][ix_lI] := spp[sp1,ix,sp2]

spOuter /: dot[spOuter[sp1_spA,sp2_spB],spOuter[sp3_spA,sp4_spB]] := 2*spp[sp1,sp3]*spp[sp4,sp2]

spOuter /: dot[vec_,spOuter[sp1_?spQ,sp2_?spQ]] /; vectorQ[vec] || Head[vec]===lI := spp[sp1,vec,sp2]
spOuter /: dot[spOuter[sp1_?spQ,sp2_?spQ],vec_] /; vectorQ[vec] || Head[vec]===lI := spp[sp1,vec,sp2]

(* spOuter /: sp[spOuter[sp1_,sp2_],sp3__] := sp1*sp[sp2,sp3]
spOuter /: sp[sp1__,spOuter[sp2_,sp3_]] := sp[sp1,sp2]*sp3 *)
spOuter /: spp[sp1__,spOuter[sp2_,sp3_],sp4__] := 2spp[sp1,sp2]*spp[sp3,sp4]+2spp[sp1,sp3]*spp[sp2,sp4]

spOuter[args___] /; Length[{args}]=!=2 := (Message[spOuter::argx,Length[{args}]]; $Failed)

spOuter /: MakeBoxes[spOuter[spA[p1_],spB[p2_]],TraditionalForm] :=
	If[AllTrue[{p1,p2},ToString[Head[#]]==$extMomLabel &],
		"\[LeftAngleBracket]"<>ToString[p1[[1]]]<>"|"<>"\[Gamma]"<>"|"<>ToString[p2[[1]]]<>"]",
		FormBox[RowBox[{"\[LeftAngleBracket]",MakeBoxes[p1,TraditionalForm],"|"<>"\[Gamma]"<>"|",MakeBoxes[p2,TraditionalForm],"]"}], TraditionalForm]
	]
(* spOuter /: MakeBoxes[spOuter[spB[p1_],spA[p2_]],TraditionalForm] :=
	If[AllTrue[{p1,p2},ToString[Head[#]]==$extMomLabel &],
		"|"<>ToString[p1[[1]]]<>"]"<>"\[LeftAngleBracket]"<>ToString[p2[[1]]]<>"|",
		FormBox[RowBox[{"|",MakeBoxes[p1,TraditionalForm],"]"<>"\[LeftAngleBracket]",MakeBoxes[p2,TraditionalForm],"|"}], TraditionalForm]
	] *)


toSpinors[expr_,vec_?tensorQ] := toSpinors[expr,{vec}]
toSpinors[expr_,vecs_] /; If[TrueQ[Quiet[AllTrue[vecs,tensorDim[#]===4 &]]],True, Message[toSpinors::ssle]; False] := expr //. {
		dot[p1_,p2_] /; (MemberQ[vecs,p1] && MemberQ[vecs,p2]) :> spp[spA[p1],spA[p2]]*spp[spB[p2],spB[p1]]/2,
		eps[p1___,p2_,p3___] /; (MemberQ[vecs,p2] && Length[{p1,p2,p3}]==4) :> Power[-1,Length[{p1}]]/(4I)*(spp[spB[p2],p1,p3,spA[p2]]-spp[spA[p2],p1,p3,spB[p2]]),
		spp[p1__,p2_,p3__] /; MemberQ[vecs,p2] :> spp[p1,spA[p2]]*spp[spB[p2],p3]+spp[p1,spB[p2]]*spp[spA[p2],p3],
		spp[p1_,p2_,p3__] /; MemberQ[vecs,p1] :> spA[p1]spp[spB[p1],p2,p3]+spB[p1]spp[spA[p1],p2,p3],
		spp[p1__,p2_,p3_] /; MemberQ[vecs,p3] :> spp[p1,p2,spA[p3]]spB[p3]+spp[p1,p2,spB[p3]]spA[p3]
}


(* ::Subsection::Closed:: *)
(*Free Indices*)


SetAttributes[lI,NHoldAll]
SetAttributes[ui,NHoldAll]
SetAttributes[di,NHoldAll]
(* SetAttributes[spI,NHoldAll] *)


(* raisedIndexQ[index_lI] := !MemberQ[index,-1]
raisedIndexQ[index_spI] := !MemberQ[index,-1]
raisedIndexQ[__] := (Message[raisedIndexQ::notindex]; $Failed) *)


contract[expr_] := FixedPoint[Expand[#] //. {
	Power[p1_?vectorQ[ix_lI],2] :> dot[p1],
	p1_?vectorQ[lI[i_]]*p2_?vectorQ[lI[i_]] :> dot[p1,p2],
	vec_?vectorQ[ix_lI]*spp[p1__,ix_lI,p2__] :> spp[p1,vec,p2],
	vec_?vectorQ[ix_lI]*eps[p1___,ix_lI,p2___] :> eps[p1,vec,p2],
	vec_?vectorQ[ix_lI]*dot[ix_lI,p__] :> dot[vec,p],
	vec_?vectorQ[ix_lI]*dot[p__,ix_lI] :> dot[p,vec],

	f_?tensorQ[ix1_lI,ix2_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[p,f,ix2],
	dot[p_,ix2_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[ix1,f,p],

	dot[p1__,ix_lI]*dot[ix_lI,p2__] :> dot[p1,p2],

	eta[ix1_lI,ix2_lI]*spp[p1__,ix2_lI,p2__] :> spp[p1,ix1,p2],
	eta[ix1_lI,ix2_lI]*eps[p1___,ix2_lI,p2___] :> eps[p1,ix1,p2],
	eta[ix1_lI,ix2_lI]*eta[ix2_lI,ix3_lI] :> eta[ix1,ix3],
	eta[ix1_lI,ix2_lI]*f_?tensorQ[p1___,ix2_lI,p2___] :> f[p1,ix1,p2],

	eps[___,lI[i_],___,lI[i_],___] :> 0,
	eta[ix_lI,ix_lI] :> $globalDim,
	mu[ix_lI,ix_lI] :> $globalDim-4,
	Power[eta[ix1_lI,ix2_lI],2] :> $globalDim,
	Power[mu[ix1_lI,ix2_lI],2] :> $globalDim-4,

	spp[sp1_spA,lI[i_],sp2_spB]spp[sp3_spA,lI[i_],sp4_spB] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp1_spA,lI[i_],sp2_spB]spp[sp4_spB,lI[i_],sp3_spA] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp2_spB,lI[i_],sp1_spA]spp[sp3_spA,lI[i_],sp4_spB] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp2_spB,lI[i_],sp1_spA]spp[sp4_spB,lI[i_],sp3_spA] :> 2spp[sp1,sp3]spp[sp4,sp2],
 	eps[p1___,lI[i_],p2___]*spp[p3___,lI[i_],p4___] :> toSpinors[eps[p1,lI[i],p2],{p1,p2}]*spp[p3,lI[-i],p4],

 	spp[sp1__,lI[i_],lI[i_],sp2__] :> 4spp[sp1,sp2],

 	mu[x_,ix_lI]*eta[ix_lI,y_] :> -mu[x,y],
 	mu[x_,ix_lI]*mu[ix_lI,y_] :> mu[x,y],
 	mu[x_,ix_lI]*p_?vectorQ[ix_lI] :> mu[x,p],

 	f_?tensorQ[ix_lI, ix_lI] :> tr[f],

 	f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],

 	v1_?vectorQ[ix1_lI]*f_?tensorQ[ix1_lI, ix2_lI]*v2_?vectorQ[ix2_lI] :> dot[v1, f, v2],

 	dot[ix_lI,args__,ix_lI] :> tr[args],
     f1_?vectorQ[ix1_lI]*f2_?vectorQ[ix1_lI] :> dot[f1,f2],
 	(* Rules only for symmetric rank-2 tensors! *)
 	(*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix3_lI, ix2_lI] :> dot[ix1,f1,f2,ix3],
 	f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],
 	f_?tensorQ[ix2_lI,ix1_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix2_lI,ix1_lI] :> dot[p,f,ix2],*)
	(* Rules only for general rank-2 tensors! *)
     f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix3_lI, ix2_lI] :> dot[ix3,f2,f1,ix1],
 	f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],
 	f_?tensorQ[ix1_lI,ix2_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[p,f,ix2],

	dot[p1_,ix_lI]*dot[p2__,ix_lI] :> dot[p2,p1],
	dot[ix_lI,p1_]*dot[ix_lI,p2__] :> dot[p1,p2]


(*  	f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix1_lI, ix2_lI] /; AllTrue[{f1, f2}, tensorRank[#] === 2 &] :> tr[f1, f2],

 	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix2_lI, ix3_lI]*v2_?vectorQ[ix3_lI] /; AllTrue[{f1,f2}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,v2],
 	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix2_lI, ix3_lI]*v2_?vectorQ[ix3_lI] /; AllTrue[{f1, f2}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix3_lI, ix4_lI]*f3_?tensorQ[ix5_lI, ix6_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix3_lI, ix4_lI]*f3_?tensorQ[ix6_lI, ix5_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix4_lI, ix3_lI]*f3_?tensorQ[ix5_lI, ix6_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix4_lI, ix3_lI]*f3_?tensorQ[ix6_lI, ix5_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix3_lI, ix4_lI]*f3_?tensorQ[ix5_lI, ix6_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix3_lI, ix4_lI]*f3_?tensorQ[ix6_lI, ix5_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix4_lI, ix3_lI]*f3_?tensorQ[ix5_lI, ix6_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2],
	v1_?vectorQ[ix1_lI]*f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix4_lI, ix3_lI]*f3_?tensorQ[ix6_lI, ix5_lI]*v2_?vectorQ[ix7_lI] /; AllTrue[{f1,f2,f3}, tensorRank[#] === 2 &] :> dot[v1,f1,f2,f3,v2] *)
} &, expr]


contractAnti[expr_] := FixedPoint[Expand[#] //. {
	Power[p1_?vectorQ[ix_lI],2] :> dot[p1],
	p1_?vectorQ[lI[i_]]*p2_?vectorQ[lI[i_]] :> dot[p1,p2],
	vec_?vectorQ[ix_lI]*spp[p1__,ix_lI,p2__] :> spp[p1,vec,p2],
	vec_?vectorQ[ix_lI]*eps[p1___,ix_lI,p2___] :> eps[p1,vec,p2],
	vec_?vectorQ[ix_lI]*dot[ix_lI,p__] :> dot[vec,p],
	vec_?vectorQ[ix_lI]*dot[p__,ix_lI] :> dot[p,vec],

	f_?tensorQ[ix1_lI,ix2_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[p,f,ix2],
	dot[p_,ix2_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[ix1,f,p],

	dot[p1__,ix_lI]*dot[ix_lI,p2__] :> dot[p1,p2],

	eta[ix1_lI,ix2_lI]*spp[p1__,ix2_lI,p2__] :> spp[p1,ix1,p2],
	eta[ix1_lI,ix2_lI]*eps[p1___,ix2_lI,p2___] :> eps[p1,ix1,p2],
	eta[ix1_lI,ix2_lI]*eta[ix2_lI,ix3_lI] :> eta[ix1,ix3],
	eta[ix1_lI,ix2_lI]*f_?tensorQ[p1___,ix2_lI,p2___] :> f[p1,ix1,p2],

	eps[___,lI[i_],___,lI[i_],___] :> 0,
	eta[ix_lI,ix_lI] :> $globalDim,
	mu[ix_lI,ix_lI] :> $globalDim-4,
	Power[eta[ix1_lI,ix2_lI],2] :> $globalDim,
	Power[mu[ix1_lI,ix2_lI],2] :> $globalDim-4,

	spp[sp1_spA,lI[i_],sp2_spB]spp[sp3_spA,lI[i_],sp4_spB] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp1_spA,lI[i_],sp2_spB]spp[sp4_spB,lI[i_],sp3_spA] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp2_spB,lI[i_],sp1_spA]spp[sp3_spA,lI[i_],sp4_spB] :> 2spp[sp1,sp3]spp[sp4,sp2],
	spp[sp2_spB,lI[i_],sp1_spA]spp[sp4_spB,lI[i_],sp3_spA] :> 2spp[sp1,sp3]spp[sp4,sp2],
 	eps[p1___,lI[i_],p2___]*spp[p3___,lI[i_],p4___] :> toSpinors[eps[p1,lI[i],p2],{p1,p2}]*spp[p3,lI[-i],p4],

 	spp[sp1__,lI[i_],lI[i_],sp2__] :> 4spp[sp1,sp2],

 	mu[x_,ix_lI]*eta[ix_lI,y_] :> -mu[x,y],
 	mu[x_,ix_lI]*mu[ix_lI,y_] :> mu[x,y],
 	mu[x_,ix_lI]*p_?vectorQ[ix_lI] :> mu[x,p],

 	f_?tensorQ[ix_lI, ix_lI] :> tr[f],

 	f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],

 	v1_?vectorQ[ix1_lI]*f_?tensorQ[ix1_lI, ix2_lI]*v2_?vectorQ[ix2_lI] :> dot[v1, f, v2],

 	dot[ix_lI,args__,ix_lI] :> tr[args],
     f1_?vectorQ[ix1_lI]*f2_?vectorQ[ix1_lI] :> dot[f1,f2],
     eta[ix1_lI,ix2_lI]*dot[ix1_lI,p__]:>dot[ix2,p],
     eta[ix1_lI,ix2_lI]*dot[p__,ix2_lI]:>dot[p,ix1],
 	(* Rules only for symmetric rank-2 tensors! *)
 	(*f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix3_lI, ix2_lI] :> dot[ix1,f1,f2,ix3],
 	f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],
 	f_?tensorQ[ix2_lI,ix1_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix2_lI,ix1_lI] :> dot[p,f,ix2],*)
	(* Rules only for general rank-2 tensors! *)
     f1_?tensorQ[ix2_lI, ix1_lI]*f2_?tensorQ[ix3_lI, ix2_lI] :> dot[ix3,f2,f1,ix1],
 	f1_?tensorQ[ix1_lI, ix2_lI]*f2_?tensorQ[ix2_lI, ix3_lI] :> dot[ix1,f1,f2,ix3],
 	f_?tensorQ[ix1_lI,ix2_lI]*dot[ix2_lI,p__] :> dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?tensorQ[ix1_lI,ix2_lI] :> dot[p,f,ix2],
    dot[p1__,ix_lI]*dot[p2__,ix_lI] :>(-1)^(Length[Cases[{p2},_?antiTensorQ]]) ((dot[p1,##]&)@@Reverse@{p2}),
	dot[ix_lI,p1__]*dot[ix_lI,p2__] :> (-1)^(Length[Cases[{p1},_?antiTensorQ]]) ((dot[##,p2]&)@@Reverse@{p1}),
Power[dot[p1__,ix_lI],2]:>(-1)^(Length[Cases[{p1},_?antiTensorQ]]) ((dot[p1,##]&)@@Reverse@{p1}),
Power[dot[ix_lI,p1__],2]:>(-1)^(Length[Cases[{p1},_?antiTensorQ]]) ((dot[##,p1]&)@@Reverse@{p1}),
	(*dot[p1_,ix_lI]*dot[p2__,ix_lI] :> dot[p2,p1],
	dot[ix_lI,p1_]*dot[ix_lI,p2__] :> dot[p1,p2],*)
(* Rules only for antisymmetric rank-2 tensors! *)
         f1_?antiTensorQ[ix1_lI, ix2_lI]*f2_?antiTensorQ[ix3_lI, ix2_lI] :> - dot[ix1,f1,f2,ix3],
         f1_?antiTensorQ[ix2_lI, ix1_lI]*f2_?antiTensorQ[ix2_lI, ix3_lI] :> - dot[ix1,f1,f2,ix3],
         f_?antiTensorQ[ix2_lI,ix1_lI]*dot[ix2_lI,p__] :> - dot[ix1,f,p],
	dot[p__,ix1_lI]*f_?antiTensorQ[ix2_lI,ix1_lI] :> - dot[p,f,ix2]


} &, expr]



contractspp[f_]:=f//.vector_[lI[i__]] spp[a___,lI[i__],b___]:>spp[a,vector,b]
contractSp[f_]:=f//.dot[lI[i_],g__,lI[j_]]spp[a___,lI[i_],lI[j_],b___]:>sp[a,CenterDot[g],b]//.dot[g__,lI[i__]]spp[a___,lI[i__],b___]:>sp[a,CenterDot[g],b]//.dot[lI[i__],g__]spp[a___,lI[i__],b___]:>sp[a,CenterDot[g],b]//.dot[lI[i__],g__]sp[a___,lI[i__],b___]:>sp[a,CenterDot[g],b]//.dot[g__,lI[i__]]sp[a___,lI[i__],b___]:>sp[a,CenterDot[g],b]


(*contractv2[expr_] /; FreeQ[expr, lI, Infinity] := expr 
contractv2[expr_lI] := expr
contractv2[expr_Times] := Module[
	{temp, recursive},
	(* get top position of lI's inside expr *)
	temp = topPosLI[expr];
	(* top position of lI's that are only contained in a single subproduct *)
	recursive = DeleteDuplicates[{First[#]}&/@Cases[temp, {_}]];
	If[Length[recursive] =!= 0,
		(*then*)
		contractTwo[MapAt[contractv2, expr, recursive]], 
		(*else*)
		contractTwo[expr, DeleteCases[temp,{_}]]
	]
]
contractv2[Power[expr_, 2]] := contractv2 /@ contractAtom[expr, expr] (* TODO: find more efficient solution here? *)
eps /: contractv2[eps[___, lI[i_], ___, lI[i_], ___]] := 0
eta /: contractv2[eta[lI[i_],lI[i_]]] := $globalDim
mu /: contractv2[mu[lI[i_],lI[i_]]] := $globalDim - 4
contractv2[expr_] := contractv2 /@ expr*)


(* Todo: build in a check for badly-formed index structures? *)
(* Todo: implementation for higher-rank tensors *)
dressIndex[expr_?(MemberQ[{Times,Plus},Head[#]]&),index_] := dressIndex[#,index]&/@expr
dressIndex[vec_?vectorQ,index_lI] := vec[index]
dressIndex[expr_,_] := expr
dressIndex[ten_?tensorQ,index_List] := ten@@index


freeIndices[expr_Times] := Sort[Join @@ freeIndices /@ List @@ expr] //. {ix1___, ix2_, ix2_, ix3___} :> {ix1, ix3}
freeIndices[expr_Plus] := freeIndices[expr[[1]]]
freeIndices[vec_?tensorQ[ix__]] := Sort[Cases[{ix}, lI[_]]]
(* freeIndices[vec_?(MemberQ[{sp,eps,dot}, Head[#]] &)] := Sort[Join @@ freeIndices /@ List @@ Cases[vec, _lI]] //. {ix1___, ix2_, ix2_, ix3___} :> {ix1, ix3} *)
freeIndices[vec_?(MemberQ[{spp,eps,dot}, Head[#]] &)] := Sort[Cases[vec, _lI]] //. {ix1___, ix2_, ix2_, ix3___} :> {ix1, ix3}
freeIndices[_] := {}


(* Todo: implementation for higher-rank tensors *)
exposeIndex[expr_,vec_?vectorQ,lI[index_]] := expr //. {
	dot[vec,p_] :> vec[lI[index]]*p[lI[index]],
	spp[p1__,vec,p2__] :> vec[lI[index]]*spp[p1,lI[index],p2],
	eps[p1___,vec,p2___] :> vec[lI[index]]*eps[p1,lI[index],p2]
}


(* Todo: implementation for lower-rank tensors *)
(* Todo: for now, this assumes that only one instance of the relevant tensor appears in the expression. *)

indexCoefficient[expr_Plus, tens_] := indexCoefficient[#, tens] & /@ expr
indexCoefficient[expr_Times, tens_] := indexCoefficient[#, tens] & /@ expr
indexCoefficient[Power[expr_, n_], tens_] := Power[indexCoefficient[expr, tens], n]
indexCoefficient[tr[p1___,f_,p2___],f_[ix1_lI,ix2_lI]] := dot[ix2,p2,p1,ix1]
indexCoefficient[dot[p1__,f_,p2__],f_[ix1_lI,ix2_lI]] := dot[p1,ix1]*dot[ix2,p2]
indexCoefficient[expr_,_] := expr


(* ::Subsection::Closed:: *)
(*Other*)


parkeTaylor[mom_List] := parkeTaylor@@mom
parkeTaylor[mom__] (* /; If[TrueQ[AllTrue[{mom},tensorDim[#]==4 &]],True,Message[parkeTaylor::ssle]; False] *) :=
	1/(Times@@Append[(spp[spA[{mom}[[#]]],spA[{mom}[[#+1]]]]&)/@Range[Length[{mom}]-1],spp[spA[{mom}[[-1]]],spA[{mom}[[1]]]]])


basisExpand[expr_, basis_List] := Catch[Module[{index = freeIndices[expr], x = Table[Unique[], {Length[basis]}]},
	Do[If[freeIndices[i] =!= index, Message["basisExpand::badindices"]; Throw[$Failed];];, {i, basis}];
   	If[index === {},
    	Return[basis . Solve[(0 == dot[#, (-expr + x . basis)] &)  /@ basis, x][[1, All, 2]]];,
    	Return[basis . Solve[(0 == contract[#*(-expr + x . basis)] &)  /@ basis, x][[1, All, 2]]];
	];
]
]


lgScaling[expr_Plus, vec_] := lgScaling[expr[[1]], vec]
lgScaling[expr_Times, vec_] := Plus @@ (lgScaling[#, vec] & /@ (List @@ expr))
lgScaling[Power[a_, p_], vec_] := p*lgScaling[a, vec]
lgScaling[expr_sp, vec_] := Plus @@ (lgScaling[#, vec] & /@ {expr[[1]], expr[[-1]]})
lgScaling[spA[vec_], vec_] := 1
lgScaling[spB[vec_], vec_] := -1
lgScaling[__] := 0


massDim[expr_Plus] := massDim[expr[[1]]]
massDim[expr_Times] := Plus @@ (massDim /@ (List @@ expr))
massDim[Power[a_, p_]] := p*massDim[a]
massDim[f_[args__]] /; MemberQ[{dot, mu, eps, spp}, f] := Plus @@ massDim /@ {args}
massDim[spA[p_?vectorQ]] := 1/2
massDim[spB[p_?vectorQ]] := 1/2
massDim[p_?vectorQ] := 1
massDim[__] := 0


(* ::Subsection::Closed:: *)
(*QCD Current*)


V2spBracket={J[_,f_,m_]:> spBracket[spB[\[Beta]],f/.CircleTimes[g___]:> g,spA[\[Alpha]]],J[_,f_]:> spBracket[spB[\[Beta]],f/.CircleTimes[g___]:> g,spA[\[Alpha]]]};


expandBracket[f_]:=f/.{spBracket-> spp}/.{spp-> spBracket}
NumQQ[f_]:=Cases[f/.CircleTimes-> List,Q[___]]//Length
NumP[f_]:=Cases[{f}/.CircleTimes-> List//Flatten,k[___]]//Length


\[ScriptCapitalF][a_Plus,b_]:=\[ScriptCapitalF][a,b]=(\[ScriptCapitalF][#,b]&/@a)//Expand
\[ScriptCapitalF][a_,b_Plus]:=\[ScriptCapitalF][a,b]=(\[ScriptCapitalF][a,#]&/@b)//Expand
\[ScriptCapitalF][f_ a_J,b_J]:=\[ScriptCapitalF][f a,b]=f \[ScriptCapitalF][ a,b]//Expand
\[ScriptCapitalF][ a_J,f_ b_J]:=\[ScriptCapitalF][a,f b]=f \[ScriptCapitalF][ a,b]//Expand
\[ScriptCapitalF][ g_ a_J,f_ b_J]:=\[ScriptCapitalF][g a,f b]=g f \[ScriptCapitalF][ a,b]//Expand
\[ScriptCapitalF][0,0]:=0
\[ScriptCapitalF][0,f_]:=0
\[ScriptCapitalF][f_,0]:=0


\[ScriptCapitalF][V1_J,V2_J]:=\[ScriptCapitalF][V1,V2]=Which[
NumQQ[V1[[2]]]==0&&NumQQ[V2[[2]]]== 0,F0[V1,V2],
NumQQ[V1[[2]]]==1&&NumQQ[V2[[2]]]==0,FL[V1,V2],
NumQQ[V1[[2]]]==0&& NumQQ[V2[[2]]]==1,FR[V1,V2],
True,0Khigher[V1,V2]]
F0[A_J,B_J]:=Module[{Qp=Q[A[[1]]+B[[1]]],Pa=k[A[[1]]],Pb=k[B[[1]]],a=A[[2]],b=B[[2]]},(
(DOT[b,Pa] J[A[[1]]+B[[1]],a])-(DOT[a,Pb] J[A[[1]]+B[[1]],b])-1/2 DOT[a,b] J[A[[1]]+B[[1]],Pa]+1/2 DOT[a,b] J[A[[1]]+B[[1]],Pb]-1/2 J[A[[1]]+B[[1]],a\[CircleTimes]b\[CircleTimes]Qp]+1/2 J[A[[1]]+B[[1]],b\[CircleTimes]a\[CircleTimes]Qp])/.DOT-> dot//Expand
]
FL[A_J,B_J]:=Module[{Qp=Q[A[[1]]+B[[1]]],Pa=k[A[[1]]],Pb=k[B[[1]]],a=A[[2]],b=B[[2]]},(
-(1/4)(DOT[Pa,Pa] DOT[a[[1]],b]J[A[[1]]+B[[1]],a[[2]]])+1/4 (DOT[Pa,Pa] DOT[a[[2]],b]J[A[[1]]+B[[1]],a[[1]]]))/.DOT-> dot//Expand
]
FR[A_J,B_J]:=Module[{Qp=Q[A[[1]]+B[[1]]],Pa=k[A[[1]]],Pb=k[B[[1]]],a=A[[2]],b=B[[2]]},(
1/4 (DOT[Pb,Pb] DOT[a,b[[1]]]J[A[[1]]+B[[1]],b[[2]]])-1/4 (DOT[Pb,Pb] DOT[a,b[[2]]]J[A[[1]]+B[[1]],b[[1]]]))/.DOT-> dot//Expand
]


PG[f1_p,f2_p]:=dot[P[(f1+f2)/.P[f_]:> f],P[(f1+f2)/.P[f_]:> f]]P[(f1+f2)/.P[f_]:> f//Expand]
PG[f1_p,f2_P]:=dot[P[(f1+f2)/.P[f_]:> f],P[(f1+f2)/.P[f_]:> f]]P[(f1+f2)/.P[f_]:> f//Expand]
PG[f1_P,f2_p]:=dot[P[(f1+f2)/.P[f_]:> f],P[(f1+f2)/.P[f_]:> f]]P[(f1+f2)/.P[f_]:> f//Expand]
PG[f1_P,f2_P]:=dot[P[(f1+f2)/.P[f_]:> f],P[(f1+f2)/.P[f_]:> f]]P[(f1+f2)/.P[f_]:> f//Expand]
PG[Times[s1_,f1_p],f2_p]:=s1 PG[f1,f2]
PG[f2_p,Times[s1_,f1_p]]:=s1 PG[f2,f1]
PG[Times[s1_,f1_P],f2_p]:=s1 PG[f1,f2]
PG[f2_p,Times[s1_,f1_P]]:=s1 PG[f2,f1]
PG[Times[s1_,f1_p],f2_P]:=s1 PG[f1,f2]
PG[f2_P,Times[s1_,f1_p]]:=s1 PG[f2,f1]
PG[Times[s1_,f1_P],f2_P]:=s1 PG[f1,f2]
PG[f2_P,Times[s1_,f1_P]]:=s1 PG[f2,f1]
PG[Times[s2_,f2_P],Times[s1_,f1_P]]:=s1 s2 PG[f2,f1]


propG[f1_P,f2_P]:=dot[p@@({List@@f1,List@@f2}//Flatten)]P@@({List@@f1,List@@f2}//Flatten)
propG[Times[s1_,f1_P],f2_P]:=s1 propG[f1,f2]
propG[f2_P,Times[s1_,f1_P]]:=s1 propG[f2,f1]
propG[Times[s2_,f2_P],Times[s1_,f1_P]]:=s1 s2 propG[f2,f1]


X2Prop[f_X]:=(f[[1]]/.(i_:> P[i]/;IntegerQ[i])/.X->propG)(f[[2]]/.(i_:> P[i]/;IntegerQ[i])/.X->propG)/.P[ft__]:>1


JQCD[gg_List]:=If[Length[gg]==1,EJ[gg[[1]]],Module[{numerVector,PPVector,Amp,bp},bp=BinaryProduct[gg];
numerVector=Table[((bp[[ii]]/.i_:> J[p[i],\[Epsilon][p[i]]]/;IntegerQ[i]/.X-> \[ScriptCapitalF]/.V2spBracket/.k[f_]:> f/.P[f_]:> f(*//expandBracket*))/.spBracket[spB[\[Beta]],f___,Q[g_],spA[\[Alpha]]]:>0(*//expandBracket*))/.spBracket[spB[\[Beta]],f___,Q[g_],spA[\[Alpha]]]:>0//Expand,{ii,Length@bp}];
(*Print[numerVector];*)
PPVector=bp/.(i_:> p[i]/;IntegerQ[i])/.X-> PG/.dot[P[f_],f1_]:> dot[k[f],f1]/.dot[f1_,P[f_]]:> dot[f1,k[f]]/.P[f_]:> 1(*/.P[Sum[p[gg[[i]]],{i,Length[gg]}]]\[RuleDelayed] 1/.dot[1,1]\[RuleDelayed] 1*);
Amp=(Sum[numerVector[[i]]/PPVector[[i]],{i,Length@PPVector}]/.P[f_]:> f/.k[f_]:> f/.spBracket[spB[\[Beta]],f___,Q[g_],spA[\[Alpha]]]:>0(*//expandBracket*))/.spBracket[spB[\[Beta]],f_,spA[\[Alpha]]]:> dot[f,v]
]]


FP[f_List,vv_]:=1/(2dot[(p/@f)//Total,vv])
BP[f_List]:=1/dot[(p/@f)//Total,(p/@f)//Total]
EJ[i_/;IntegerQ[i]]:=dot[\[Epsilon][p[i]],v]


Jh[f_List,vv_]:=Module[{len=Length[f]},JQCD[Range[len]]]/.v-> vv/.\[Epsilon][p[i_]]:> \[Epsilon][i]/.{p[i_]:> f[[i,1]],\[Epsilon][i_]:> f[[i,2]]}


(* ::Subsection:: *)
(*alpha' higher orders*)


cyclePermuteSet[od__]:=Module[{odlist={od}},Table[Join[odlist[[i;;-1]],odlist[[1;;i-1]]],{i,Length@odlist}]]


trFLA[od__,i2_ ]:=Module[{lnest},lnest=(ExpandNCM/@(NonCommutativeMultiply/@(List@@(\[CapitalOmega][od]))))/.\[FivePointedStar]->Sequence/.NonCommutativeMultiply->F//Total;
lnest=lnest/.F[ii__]:>tr[F[ii,i2]]
]


WFun0[n_]:=Module[{effpt,monos,res},If[n>=4,effpt=SetPartitionsOrdered[n]//Cases[{{i1_,g___,i2_},f___,{i3_,h___,i4_}}];
monos=Table[W@@@effpt[[ii]],{ii,Length@effpt}]/.W[i_]:>Sequence[]//.{g___,W[h1__,i_],W[j_,h2__],f___}:>{g,W[h1,i],dot[p[i],F@@Table[ii,{ii,i+1,j-1}],p[j]],W[j,h2],f}/.F[]->Sequence[];
monos=Times@@@monos,monos={0}];
res=W@@Range[n]+Total[Join[monos,{0}]]
]
WFun0[od__,ode_]:=WFun0[Length[{od,ode}]]/.dot[f__]:>(dot[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W[f__]:>(W[f]/.i_?IntegerQ:>{od,ode}[[i]])


WFun[i1_]:=0
WFun[i1_,i2_]:=\[Alpha] tr[F[i1],F[i2]] +\[Alpha]^2 dot[p[i1,i2]]tr[F[i1],F[i2]] 
WFun[i1_,od__,i2_]:=Module[{lnest,cyc,cycs,res1,allnest,signs,res2,tau,res3,psMax},
lnest=trFLA[i1,od,i2];
res1=\[Alpha] lnest +\[Alpha]^2 dot[p[i1,od,i2]]lnest;
(*Print[lnest];*)
allnest=List@@lnest;
signs=allnest/.tr[_]:>1;
allnest=allnest/.Times[a_,tr[f_]]:>List@@f/.tr[f_]:>List@@f;
(*Print[allnest];*)
res2=Sum[
cycs=cyclePermuteSet@@allnest[[kk]];
(*Print[cycs];*)
\[Alpha]^2 signs[[kk]]Table[cyc=cycs[[jj]];
Sum[dot[p[cyc[[-1]]],F@@cyc[[1;;jj1]],p[cyc[[jj1+1]]]]trFLA@@cyc[[jj1+1;;-1]],{jj1,1,-2+Length@cyc}],{jj,Length@cycs}]/.p[]->0//Total,{kk,Length[allnest]}];
tau={i1,od,i2};
If[Length[tau]<4,res3=0,
res3=\[Alpha]^2 Sum[trFLA@@tau[[1;;ii]]dot[p[tau[[ii]]],F@@tau[[ii+1;;jj-1]],p[tau[[jj]]]]trFLA@@tau[[jj;;-1]],{ii,2,Length[tau]-2},{jj,ii+1,Length[tau]-1}]/.F[]->Sequence[]];
res1+res2+res3
 ]


T2FF3F4[f_T]:=Module[{fls,od,odi,odl,n,phat,leftv,rightv,num,num1,den,refs,sc,refg},
fls=List@@f;
od=fls;
n=2+Length@(od//Flatten);
leftv=Range[Length@od];
     leftv[[1]]=v[0];
Table[
odl=Flatten[od[[1;;(i-1)]]];
leftv[[i]]=p@@Select[odl,#<od[[i,1]]&],{i,2,Length@od}];num=Product[odi=od[[ii]];
(*Print[odi];*)
dot[leftv[[ii]],F[Sequence@@odi],v[0]]+Sum[dot[leftv[[ii]],F[Sequence@@odi[[1;;j1-1]]],p[odi[[j1]]]]W@@odi[[j1;;j2]]dot[p[odi[[j2]]],F[Sequence@@odi[[j2+1;;-1]]],v[0]],{j2,1,Length@odi},{j1,1,j2-1}],{ii,Length[od]}]/.dot[g1__,F[],g2__]:>dot[g1,g2];
den=dot[v[0],p[1]]Product[ dot[v[0],(p@@Flatten[Join[{},od[[1;;(i-1)]]]])],{i,2,Length@od}];
  (num/den)
]


MultiTrace[n_?IntegerQ]:=Module[{effpt,monos,lengthList},
If[n>=2,
effpt=Drop[SetPartitionsOrdered[n],1];
effpt=Drop[effpt,-1];
lengthList=Table[1/Length[effpt[[ii]]],{ii,Length@effpt}];
monos=Table[(W@@@effpt[[ii]]),{ii,Length@effpt}]/.W[i_]:>Sequence[]//.{g___,W[h1__,i_],W[j_,h2__],f___}:>{g,W[h1,i],dot[p[i],F@@Table[ii,{ii,i+1,j-1}],p[j]],W[j,h2],f}/.{W[i_,h1__],f___,W[h2__,j_]}:>{W[i,h1],f,W[h2,j],dot[p[j],F@@Join[Table[ii,{ii,j+1,n}],Table[ii,{ii,1,i-1}]],p[i]]}/.F[]->Sequence[]/.{W[i_,h___,j_]}:>{W[i,h,j],dot[p[j],F@@Join[Table[ii,{ii,j+1,n}],Table[ii,{ii,1,i-1}]],p[i]]};
monos=Times@@@monos;monos=Table[lengthList[[ii]]*monos[[ii]],{ii,Length@effpt}],monos={0}]
]
MultiTrace[od__?IntegerQ,ode_?IntegerQ]:=MultiTrace[Length[{od,ode}]]/.dot[f__]:>(dot[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W[f__]:>(W[f]/.i_?IntegerQ:>{od,ode}[[i]])
MultiTrace[n_?IntegerQ,pts_List]:=Module[{effpt,monos},
If[n>=2,
effpt=pts;
monos=Table[W@@@effpt[[ii]],{ii,Length@effpt}]/.W[i_]:>Sequence[]//.{g___,W[h1__,i_],W[j_,h2__],f___}:>{g,W[h1,i],dot[p[i],F@@Table[ii,{ii,i+1,j-1}],p[j]],W[j,h2],f}/.{W[i_,h1__],f___,W[h2__,j_]}:>{W[i,h1],f,W[h2,j],dot[p[j],F@@Join[Table[ii,{ii,j+1,n}],Table[ii,{ii,1,i-1}]],p[i]]}/.F[]->Sequence[]/.{W[i_,h___,j_]}:>{W[i,h,j],dot[p[j],F@@Join[Table[ii,{ii,j+1,n}],Table[ii,{ii,1,i-1}]],p[i]]};
monos=Times@@@monos,monos={0}]
]
MultiTrace[od__?IntegerQ,ode_?IntegerQ,pts_List]:=MultiTrace[Length[{od,ode}],pts]/.dot[f__]:>(dot[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W[f__]:>(W[f]/.i_?IntegerQ:>{od,ode}[[i]])


WFunDF2[i1_]:=0
WFunDF2[od__]:=WFunDF2[od]=Module[{res1,res2,res3,cycs,cyc,lnest,signs,allnest},lnest=trFLA[od]+zero;
(*Print[lnest];*)
allnest=DeleteCases[List@@lnest,zero];
signs=allnest/.tr[_]:>1;
allnest=allnest/.Times[a_,tr[f_]]:>List@@f/.tr[f_]:>List@@f;
(*Print[allnest];*)
res1=\[Alpha]/(1-\[Alpha] dot[p[od]]) W0[od];
res3=Sum[cycs=cyclePermuteSet@@allnest[[kk]];signs[[kk]]*\[Alpha]/(1-\[Alpha] dot[p[od]]) Total[Table[cyc=cycs[[jj]];MultiTrace@@cyc,{jj,Length@cycs}]//Flatten],{kk,Length[allnest]}];
res1+res3
]
WFunDF2[od__,pts_List]:=WFunDF2[od,pts]=Module[{res1,res2,res3,cycs,cyc,lnest,signs,allnest},lnest=trFLA[od]+zero;
(*Print[lnest];*)
allnest=DeleteCases[List@@lnest,zero];
signs=allnest/.tr[_]:>1;
allnest=allnest/.Times[a_,tr[f_]]:>List@@f/.tr[f_]:>List@@f;
(*Print[allnest];*)
res3=Sum[cycs=cyclePermuteSet@@allnest[[kk]];signs[[kk]]*\[Alpha]/(1-\[Alpha] dot[p[od]]) Total[Table[cyc=cycs[[jj]];MultiTrace[Sequence@@cyc,pts],{jj,Length@cycs}]//Flatten],{kk,Length[allnest]}];
res3
]


SingleTrace[n_?IntegerQ]:=Module[{effpt,monos,lengthList},
If[n>=2,
effpt=Drop[SetPartitionsOrdered[n],1];
effpt=Drop[effpt,-1]//Cases[{{f__},g__}/;Length[{f}]+Length[{g}]==n];
monos=Table[(W@@@effpt[[ii]]),{ii,Length@effpt}]/.W[i_]:>Sequence[]/.{W[i_,h___,j_]}:>{W[i,h,j],dot[p[j],F@@Join[Table[ii,{ii,j+1,n}],Table[ii,{ii,1,i-1}]],p[i]]};
monos=Times@@@monos;monos=Table[monos[[ii]],{ii,Length@effpt}],monos={0}]
]
SingleTrace[od__?IntegerQ,ode_?IntegerQ]:=SingleTrace[Length[{od,ode}]]/.dot[f__]:>(dot[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W[f__]:>(W[f]/.i_?IntegerQ:>{od,ode}[[i]])


WFunDF2SingleTr[i1_]:=0
WFunDF2SingleTr[od__]:=WFunDF2SingleTr[od]=Module[{res1,res2,res3,cycs,cyc,lnest,signs,allnest},lnest=trFLA[od]+zero;
(*Print[lnest];*)
allnest=DeleteCases[List@@lnest,zero];
signs=allnest/.tr[_]:>1;
allnest=allnest/.Times[a_,tr[f_]]:>List@@f/.tr[f_]:>List@@f;
(*Print[allnest];*)
res1=\[Alpha]/(1-\[Alpha] dot[p[od]]) W0[od];
res3=Sum[cycs=cyclePermuteSet@@allnest[[kk]];signs[[kk]]*\[Alpha]/(1-\[Alpha] dot[p[od]]) Total[Table[cyc=cycs[[jj]];SingleTrace@@cyc,{jj,Length@cycs}]//Flatten],{kk,Length[allnest]}];
res1+res3
]


replace[otherDDMOrder_,n_]:=Table[(p/@Range[n-1])[[i]]-> otherDDMOrder[[i]],{i,n-1}];
replace[otherDDMOrder_List]:=Module[{gg,ggod},gg=p/@otherDDMOrder;ggod=Sort[gg];Table[(ggod)[[i]]-> (gg)[[i]],{i,Length@gg}]]


W0Relations[gluons_]:=Module[{vponshell,onshell2,numW123,numWNC123,comp,solW,n},n=gluons+2;
vponshell={dot[v,p[1]]:> (-Sum[dot[v,p[j]],{j,2,n-2}])}/.v->v[0];
onshell2={dot[p[i_],p[i_]]:>0,
dot[p[i_],\[Epsilon][p[i_]]]:>0,dot[\[Epsilon][p[i_]],p[i_]]:>0};
numW123=(W@@Range[n-2])dot[p[p[n-2]],v[0]]/.W[f__]:>W@@(p/@{f})/.F[f__]:>F@@(p/@{f});
numWNC123=(ExpandNCM/@(BinaryProduct[n-2][[-1]]/.X-> NC)/.NonCommutativeMultiply-> T)/.T[f__]:>(numW123/.replace[{f}])/.p[i_]:>i;
comp=numWNC123-(n-2)(numW123/.p[i_]:>i);
comp=comp/.vponshell//Collect[#,dot[__]]&;
comp=Table[0==(comp[[ii]])/.dot[p[i_],v[0]]:>1/.W->W0,{ii,Length@comp}];
comp
]
W0Relations[gid1_,gids__]:=Module[{idlist},idlist={gid1,gids};W0Relations[Length[idlist]]/.W0[f__]:>(W0@@{f}/.i_?IntegerQ:>idlist[[i]])]


repW0Relations[gids__]:=Module[{idlist,allorder,eqns,solW},idlist={gids};allorder=Permutations[idlist];eqns=Table[W0Relations@@allorder[[ii]],{ii,Length@allorder}]//Flatten//Union;
solW=Solve[eqns][[1]]
]


repW0=Join[{W0[i1_,g__,i2_]:>(W0[i1,g,i2]/.repW0Relations[i1,g,i2])/;i1>Min[{g}]||i2<Max[{g}]},{W0[i1_,i2_]:>Sort[W0[i1,i2]]}];


repWp=Join[{W[i1_,g__,i2_]:>((W[i1,g,i2]/.(repW0Relations[i1,g,i2]/.W0->W))/;i1>Min[{g}]||i2<Max[{g}])},{W[i1_,i2_]:>Sort[W[i1,i2]]}];


WFun2YM[ng_?IntegerQ]:=WFun2YM[ng]=Module[{ls=Range[ng],lsb,pts,pts2,pt,pleft,pright,res1,res2},lsb=ls[[2;;-2]];pts=KSetPartitions[lsb,2];
pts2=pts/.{f1_List,f2_List}:>{f2,f1};
pts=Join[pts,pts2,{{lsb,{}}}];
pts=pts/.{f1_List,f2_List}:>{f1,Join[{1},f2,{ng}]};
If[ng>2,res1=Sum[pt=pts[[ii]];
pleft=p@@Select[pt[[2]],#<pt[[1,1]]&];
pright=p@@Select[pt[[2]],#>pt[[1,-1]]&];
2 (\[Alpha] dot[pleft,F@@pt[[1]],pright])/(1-\[Alpha] dot[p@@ls]) W@@pt[[2]],{ii,Length@pts}],res1=0];
res2=\[Alpha] /(1-\[Alpha] dot[p@@ls]) W0@@ls;
res1+res2
]
WFun2YM[od__?IntegerQ,ode_?IntegerQ]:=WFun2YM[Length[{od,ode}]]/.dot[f__]:>(dot[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W[f__]:>(W[f]/.i_?IntegerQ:>{od,ode}[[i]])/.W0[f__]:>(W0[f]/.i_?IntegerQ:>{od,ode}[[i]])


W2basis[f_W]:=Module[{lnf=Length[f],max=Max[List@@f],min=Min[List@@f]},f/.W[g1___,max,g2___,min,g3___]:>(-1)^lnf W[Sequence@@Reverse[{g1,max,g2,min,g3}]]/.W[g1___,min,g2___,max,g3___]:>(-1)^If[{g1}==={},0,1] (-1)(-1)^(Length[{g3}]-1) ExpandNCM[(NonCommutativeMultiply[min,\[CapitalOmega]@@{g1}/.\[CapitalOmega][i_?IntegerQ]:>i/.\[FivePointedStar]->NonCommutativeMultiply,g2,\[CapitalOmega]@@Reverse[{g3}]/.\[CapitalOmega][i_?IntegerQ]:>i/.\[FivePointedStar]->NonCommutativeMultiply,max]/.\[CapitalOmega][]->Sequence[])]/.NonCommutativeMultiply->W]


W02basis[f_W0]:=Module[{lnf=Length[f],max=Max[List@@f],min=Min[List@@f]},f/.W0[g1___,max,g2___,min,g3___]:>(-1)^lnf W[Sequence@@Reverse[{g1,max,g2,min,g3}]]/.W0[g1___,min,g2___,max,g3___]:>(-1)^If[{g1}==={},0,1] (-1)(-1)^(Length[{g3}]-1) ExpandNCM[(NonCommutativeMultiply[min,\[CapitalOmega]@@{g1}/.\[CapitalOmega][i_?IntegerQ]:>i/.\[FivePointedStar]->NonCommutativeMultiply,g2,\[CapitalOmega]@@Reverse[{g3}]/.\[CapitalOmega][i_?IntegerQ]:>i/.\[FivePointedStar]->NonCommutativeMultiply,max]/.\[CapitalOmega][]->Sequence[])]/.NonCommutativeMultiply->W0]


T2YM[f_T]:=Module[{fls,od,odi,odl,n,phat,leftv,rightv,num,num1,den,refs,sc,refg},
fls=List@@f;
od=fls;
n=1+Length@(od//Flatten);
leftv=Range[Length@od];
leftv[[1]]=v[0];
Table[odl=Flatten[od[[1;;(i-1)]]];leftv[[i]]=If[Select[odl,#<od[[i,1]]&]==={},p[Min[odl]],p@@Select[odl,#<od[[i,1]]&]],{i,2,Length@od}];rightv=Range[Length@od];
rightv[[1]]=v[0];
Table[odl=Flatten[od[[1;;(i-1)]]];rightv[[i]]=If[Select[odl,#>od[[i,-1]]&]==={},p[Max[odl]],p@@Select[odl,#>od[[i,-1]]&]],{i,2,Length@od}];
num=(W0@@od[[1]])Product[odi=od[[ii]];2dot[leftv[[ii]],F[Sequence@@odi],rightv[[ii]]],{ii,2,Length[od]}];
den=(*(-1/\[Alpha]+dot[p@@(od//Flatten)])*)-1/\[Alpha]*Product[ -1/\[Alpha]+dot[(p@@Flatten[Join[{},od[[1;;(i-1)]]]])],{i,2,Length@od}];
  (num/den)
]


T2YMLocal[f_T]:=Module[{fls,od,odi,odl,n,phat,leftv,rightv,num,num1,den,refs,sc,refg},
fls=List@@f;
od=fls;
n=Length@(od//Flatten);
leftv=Range[Length@od];
leftv[[1]]=v[0];
Table[odl=Flatten[od[[1;;(i-1)]]];leftv[[i]]=If[Select[odl,#<od[[i,1]]&]==={},p[Min[odl]],p@@Select[odl,#<od[[i,1]]&]],{i,2,Length@od}];rightv=Range[Length@od];
num=(dot[\[Epsilon][od[[1]][[1]]],F@@od[[1]][[2;;-2]],\[Epsilon][od[[1]][[-1]]]])Product[odi=od[[ii]];dot[leftv[[ii]],F[Sequence@@odi[[1;;-2]]],\[Epsilon][odi[[-1]]]],{ii,2,Length[od]}]/.F[]->Sequence[];
den=(*(-1/\[Alpha]+dot[p@@(od//Flatten)])*) Product[Length[Flatten[Join[{},od[[i;;-1]]]]],{i,2,Length@od}];
  (num/den)
]


(*trFLA[od__,i2_ ]:=Module[{lnest},lnest=(ExpandNCM/@(NonCommutativeMultiply/@(List@@(\[CapitalOmega][od]))))/.\[FivePointedStar]->Sequence/.NonCommutativeMultiply->F//Total;
lnest=lnest/.F[ii__]:>tr[F[ii,i2]]
]*)


(* ::Subsection:: *)
(*Ends*)


End[] (* End Private Context *)


Protect@@Names["KiHa`*"];


EndPackage[]
