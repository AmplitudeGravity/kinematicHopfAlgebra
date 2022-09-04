(* ::Package:: *)

(* Wolfram Language Package *)



(* ::Section:: *)
(*Public Initializations*)


BeginPackage["basicfun`"]
Unprotect@@Names["basicfun`*"];
(* Exported symbols added here with SymbolName::usage *) 


(* ::Subsection:: *)
(*Options*)


(* dim::usage = "dim is an option used to specify dimensions."
dim::author = "Gustav Mogull" *)


(* mass::usage = "mass is an option used to specify the on-shell mass of declared vectors."
mass::author = "Gustav Mogull" *)


(* rank::usage = "rank is an option used to specify the rank of declared tensors."
rank::author = "Gustav Mogull" *)


(* ::Subsection:: *)
(*Scalars*)


setGlobalDim::usage = "setGlobalDim[d] sets d as the global dimension."
setGlobalDim::ssle = "Protected objects cannot be set as the global dimension."
setGlobalDim::author = "Gustav Mogull"


(* ::Subsection:: *)
(*Vectors*)


tensorQ::usage = "tensorQ[p] yields True if p is a tensor, and yields False otherwise."
tensorQ::author = "Gustav Mogull"


vectorQ::usage = "vectorQ[p] yields True if p is a vector (rank-1 tensor), and yields False otherwise."
vectorQ::author = "Gustav Mogull"
tensorQR2::usage = "vectorQ[p] yields True if p is a vector (rank-1 tensor), and yields False otherwise."
tensorQR2::author = "Gustav Mogull"


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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Free Indices*)


lI::usage = "lI[i] represents the i'th Lorentz (spacetime) index."
lI::author = "Gustav Mogull"


(* spI::usage = "spI[i] represents the i'th spinor index."
spI::author = "Gustav Mogull" *)


contract::usage = "contract[expr] contracts all pairs of raised and lowered Lorentz indices in expr."
contract::author = "Gustav Mogull"


contractv2::usage = "contract[expr] contracts all pairs of raised and lowered Lorentz indices in expr."
contractv2::author = "Gregor Kaelin"


expandTensor::usage="expandTensor[dot[f]] expands all the tensor into the conponent"
expandTensor::author="Gang Chen"


dressIndex::usage = "dressIndex[expr,i] dresses the expression expr with the free Lorentz index i."
dressIndex::author = "Gustav Mogull"


freeIndices::usage = "freeIndices[expr] finds the overall free index structure in expr."
freeIndices::author = "Gustav Mogull"


exposeIndex::author = "Gustav Mogull"


indexCoefficient::author = "Gustav Mogull"


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Extrafuns*)


Diamond::usage="a distributive general product"


(* ::Section:: *)
(*Private*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Scalars*)


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


(* ::Subsection:: *)
(*Vectors*)


(* Todo: fix dimensionality of dot and eps *)
(* $declaredTensors = {{dot,$globalDim,2},{eps,4,4}}; *)
$declaredTensors = {{eta,4,2}};
$declaredTensorHeads = {{spOuter,4,1},{outer,4,2}};
(* Todo: there should be an (optional) function that declares these...? *)
$extMomLabel = "p";
$loopLabel = "l";


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
	vec /: MakeBoxes[vec[basicfun`lI[basicfun`Private`i_]],TraditionalForm] :=
		FormBox[SuperscriptBox[MakeBoxes[vec,TraditionalForm],SubscriptBox["\[Nu]",ToString[basicfun`Private`i]]],TraditionalForm];
	AppendTo[$declaredTensors,{vec,OptionValue["dim"],1}];
	If[OptionValue["verbose"],Print[ToString[vec]<>" declared as a "<>ToString[OptionValue["dim"]]<>"-dimensional vector."]];
)
declareVector[__] := $Failed


Options[declareVectorHead] := {"dim"->$globalDim,"verbose"->True}
declareVectorHead[expr_List,options___] := Scan[declareVectorHead[#,options]&,expr];
declareVectorHead[vec_Symbol,OptionsPattern[]] /; If[!protectedQ[vec],True,Message[declareVectorHead::safe]; False] := (
	If[tensorQ[vec[1]],undeclareTensorHead[vec];];
	vec /: MakeBoxes[vec[i_],TraditionalForm] := FormBox[SubscriptBox[ToString[vec],MakeBoxes[i,TraditionalForm]],TraditionalForm];
	MakeBoxes[vec[basicfun`Private`i_][basicfun`lI[basicfun`Private`j_]],TraditionalForm] :=
		FormBox[SuperscriptBox[MakeBoxes[vec[basicfun`Private`i],TraditionalForm],SubscriptBox["\[Nu]",ToString[basicfun`Private`j]]],TraditionalForm];
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


SetAttributes[Diamond,{NHoldAll}];
declareDistributive[Diamond,tensorQ];
declareDistributive[Diamond,vectorQ];


(* declareTensorPattern[expr_, pattern_] := Module[{indices = Sort[freeIndices[pattern]]},
	If[indices =!= Array[lI, Length[indices]], Message[declareTensorPattern::badpattern]; Break[];];
	If[!tensorQ[expr],declareTensor[expr,dim->4,rank->Length[indices]];];

  	expr[args__lI] := Evaluate[pattern //. lI[i_] :> Evaluate[{args}[[i]]]];
  	If[tensorRank[expr] === 1,
   		Evaluate[If[Head[expr] === Symbol, expr, Head[expr]]] /: dot[expr, vec2_?tensorQ] := contract[(expr[lI[#]]*vec2[lI[-#]] &)[Unique[]]];
   		Evaluate[If[Head[expr] === Symbol, expr, Head[expr]]] /: sp[p1__, expr, p2__] := contract[(expr[lI[#]] sp[p1, lI[-#], p2] &)[Unique[]]];
   	];
 ] *)


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


(* ::Subsection:: *)
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
spp[spA[p1_],args___,spA[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && OddQ[Length[{args}]] := 0
spp[spB[p1_],args___,spB[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && OddQ[Length[{args}]] := 0
spp[spA[p1_],args___,spB[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && EvenQ[Length[{args}]] := 0
spp[spB[p1_],args___,spA[p2_],OptionsPattern[]] /; patternFreeQ[{args}] && EvenQ[Length[{args}]] := 0

spp[___,p_?tensorQ,___] /; If[tensorRank[p]===1,False,Message[spp::badrank]; True] := $Failed

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


(* ::Subsection:: *)
(*Free Indices*)


SetAttributes[lI,NHoldAll]
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


contractv2[expr_] /; FreeQ[expr, lI, Infinity] := expr 
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
contractv2[expr_] := contractv2 /@ expr


contractTwo::usage = "Helper for contractv2: contracts same indices that appears in two different parts of Times"
contractTwo::author = "Gregor Kaelin"
contractTwo[expr_Times, positions_List : {}] := Module[
	{temp, repl = {}, tempExpr, tempScan},
	(* get top position of lI's inside expr *)
	If[positions === {},
		(* then *)
		temp = DeleteCases[topPosLI[expr], {_}];,
		(* else *)
		temp = positions;
	];
	
	(* Apply contractAtom *)
	tempExpr = List@@expr;
	Scan[
		(
			Print[tempExpr, " ", tempScan];
			tempScan = #//.repl;
			With[{contracted = contractAtom[tempExpr[[tempScan[[1]]]], tempExpr[[tempScan[[2]]]]]},
				(* handle composite return objects *)
				If[Head[contracted] === Times,
					(*then*)
					If[Length[contracted] === 2 && contracted[[1]] == -1,
						(*then*)
						tempExpr[[tempScan[[1]]]] = -contracted;
						tempExpr[[tempScan[[2]]]] = -1;
						AppendTo[repl, tempScan[[2]] -> tempScan[[1]]];,
						(*else*)
						(*restart*)
						tempExpr[[tempScan[[1]]]] = contracted;
						tempExpr[[tempScan[[2]]]] = 1;
						Return[contractTwo[Times@@tempExpr]];
					],
					(*else*)
					tempExpr[[tempScan[[1]]]] = contracted;
					tempExpr[[tempScan[[2]]]] = 1;
					AppendTo[repl, tempScan[[2]] -> tempScan[[1]]]
				];
			];
		)&,
		temp
	];
	Times @@ tempExpr
]
contractTwo[expr_] := expr


(* TODO: associate rules with individual functions via tags *)
contractAtom::usage = "Helper for contractv2: contracts 'atomic' expressions and expands sums if necessary"
contractAtom::author = "Gregor Kaelin"
SetAttributes[contractAtom, {Orderless}]
(* if we don't distribute objects automatically *)
contractAtom[vecs1_?vectorIndexSumQ, vecs2_?vectorIndexSumQ] /; ((vecs1[[1,-1]] === vecs2[[1, -1]]) && !(distributive /. Options[dot])) := dot[Head /@ vecs1, Head /@ vecs2] 
contractAtom[vecs_?vectorIndexSumQ, vec_?vectorQ[lI[i_]]] /; ((vecs[[1, -1, 1]] === i) && !(distributive /. Options[dot])) := dot[Head /@ vecs, vec] 
contractAtom[vecs_?vectorIndexSumQ, spp[p1___, lI[i_], p2___]] /; (!(distributive /. Options[spp])) := spp[p1, Head /@ vecs, p2]
contractAtom[vecs_?vectorIndexSumQ, dot[lI[i_], p__]] /; (!(distributive /. Options[dot])) := dot[Head /@ vecs, p]
contractAtom[vecs_?vectorIndexSumQ, dot[p__, lI[i_]]] /; (!(distributive /. Options[dot])) := dot[p, Head /@ vecs]
(* rest *)
contractAtom[a___, b_Plus, c___] := contractv2[Distribute[{a, b, c}, Plus, List, Plus, Times]]
contractAtom[vec1_?vectorQ[lI[i_]], vec2_?vectorQ[lI[i_]]] := dot[vec1, vec2]
contractAtom[vec_?vectorQ[lI[i_]], spp[p1___, lI[i_], p2___]] := spp[p1, vec, p2]
contractAtom[vec_?vectorQ[lI[i_]], dot[lI[i_], p__]]  := dot[vec, p]
contractAtom[vec_?vectorQ[lI[i_]], dot[p__, lI[i_]]]  := dot[p, vec]
contractAtom[eta[lI[i_], lI[j_]], f_?tensorQ[p2___, lI[i_], p3___]] := contractv2[f[p2, lI[j], p3]] (* might contain more indices that have to be contracted *)
contractAtom[mu[x_,lI[i_]],p_?vectorQ[lI[i_]]] := mu[x,p]
(*contractAtom[vec_?vectorQ[lI[i_]], f_?tensorQ[p2___, lI[i_], p3___]] := f[p2, vec, p3]*)
contractAtom[spp[sp1_spA, lI[i_], sp2_spB], spp[sp3_spA, lI[i_], sp4_spB]] := 2 spp[sp1, sp3] spp[sp4, sp2]
contractAtom[spp[sp1_spA, lI[i_], sp2_spB], spp[sp4_spB, lI[i_], sp3_spA]] := 2 spp[sp1, sp3] spp[sp4, sp2]
contractAtom[spp[sp2_spB, lI[i_], sp1_spA], spp[sp3_spA, lI[i_], sp4_spB]] := 2 spp[sp1, sp3] spp[sp4, sp2]
contractAtom[spp[sp2_spB, lI[i_], sp1_spA], spp[sp4_spB, lI[i_], sp3_spA]] := 2 spp[sp1, sp3] spp[sp4, sp2]
contractAtom[eta[x_, lI[i_]], spp[sp1__, lI[i_], sp2__]] := contractv2[spp[sp1, lI[i], sp2]] (* might contain more indices that have to be contracted *)
contractAtom[eta[x_,lI[i_]],mu[lI[i_],y_]] := contractv2[-mu[x,y]] (* might contain more indices that have to be contracted *)
contractAtom[mu[x_,lI[i_]],mu[lI[i_],y_]] := contractv2[mu[x,y]] (* might contain more indices that have to be contracted *)
(* contractAtom[eps[p1___, lI[i_], p2___], sp[p3___, lI[i_], p4___]] := contractv2[toSpinors[eps[p1, lI[i], p2], {p1, p2}] sp[p3, lI[-i], p4]](* ??? *) *)
contractAtom[expr1_, expr2_] := expr1 * expr2 (* everything else is just returned as is *)


topPosLI::usage = "Helper for contractv2: gets the top positions of lI objects (sorted and deleted duplicates)"
topPosLI::author = "Gregor Kaelin"
topPosLI[expr_] := DeleteDuplicates[
	Function[
		x, 
		Sort[DeleteDuplicates[#[[2, 1]] & /@ x]]
	] /@ GatherBy[{Part[expr, Sequence @@ #], #} & /@ Position[expr, lI[_], Infinity], First]
];


vectorIndexSumQ::usage = "Helper for contractv2: checks if expression is a sum of vectors with open indices"
vectorIndexSumQ::author = "Gregor Kaelin"
vectorIndexSumQ[expr_Plus] := (And @@ (vectorQ[Head[#]]& /@ (List @@ expr))) && MatchQ[expr[[1, -1]], lI[_]]&& (SameQ[Sequence @@ (Last /@ (List @@ expr))])
vectorIndexSumQ[_] := False


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


(* ::Subsection:: *)
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


End[] (* End Private Context *)


Protect@@Names["basicfun`*"];


EndPackage[]
