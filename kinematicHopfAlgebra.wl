(* ::Package:: *)

(*Copyright 2022,Gang Chen AmplitudeGravity/kinematicHofpAlgebra is licensed under the
GNU General Public License v3.0 .  Permissions of this strong copyleft license are 
conditioned on making available complete source code of licensed works and modifications,
which include larger works using a licensed work,under the same license.
Copyright and license notices must be preserved.Contributors provide an express grant 
of patent rights.*)
(*  Updating version can be find in  
https://github.com/AmplitudeGravity/kinematicHofpAlgebra
  *)


nice={dot[f___]:> CenterDot[f],DOT[f_,g_]:>  CenterDot[f,g],CenterDot[p[i_],\[Epsilon][j_]]:>CenterDot[\[Epsilon][j],p[i]],\[Epsilon][p[i_]]:> Subscript[\[Epsilon], i],p[i_]:> Subscript[p, i],a_[i_]?vectorQ:> Subscript[a, i],a_[i_]?tensorQ:> Subscript[a, i],spBracket-> Diamond,spB[\[Beta]]:> \!\(\*OverscriptBox[\(u\), \(_\)]\),spA[\[Alpha]]:>v,k[f_]:> f,P[f_]:> f,Q[f_]:> f,J[f1_,f2_,f3_]:> Subscript[J, f2],J[f1_,f2_]:> Subscript[J, f2]};
nice2={CenterDot[f_Plus,f2_]:> CenterDot[T@@f,f2],CenterDot[f2_,f_Plus]:> CenterDot[f2,T@@f],CircleTimes[f___,g_Plus,f2___]:> CircleTimes[f,T@@g,f2],T[f___]:>TT[{f}/.Subscript[p, i_]:>i] ,TT[f_]:> Subscript[p, ToExpression[StringJoin@@ToString/@(f)]],s[f___]:> Subscript[s, ToExpression[StringJoin@@ToString/@(f)]]};
niceF={dot[f_]:> CenterDot[f,f],dot[f__]:> CenterDot[f],p[i__]:> Subscript[p, i],F[i_]:> Subscript[F, i],a[i_]:> Subscript[a, i]};


niceT={T[f___]:>(Subscript[T, f]/. List->L),L[f1___]:>"("<>ToString/@{f1}<>")"}
niceET={a[i_]:>t^Subscript[a, i],GT[{f1___},{f2___}]:>(\!\(\*SuperscriptBox[
SubscriptBox[\(T\), \(f2\)], \("\<(\>" <> ToString /@ {f1} <> "\<)\>"\)]\)/. List->L),L[f1___]:>"("<>ToString/@{f1}<>")",ET[f1_GT,{f2__}]:>f1 tr[CenterDot[f2]]}


rmzero={dot[a_,F[i_],a_]:> 0,dot[p[],f__,v]:> 0,dot[p[],g___]:> 0,dot[g___,p[]]:> 0};
rmzeroT={T[{i_},g___]:> 0,T[f__,{1,h___},g___]:>0};


Diamond[a_Plus,b_]:=Diamond[a,b]=Plus@@(Diamond[#,b]&/@a)//Expand
Diamond[a_,b_Plus]:=Diamond[a,b]=Plus@@(Diamond[a,#]&/@b)//Expand
Diamond[c___,a_Plus,b___]:=Plus@@(Diamond[c,#,b]&/@a)//Expand
Diamond[Times[-1,a__],b___]:=-Diamond[Times[a],b]
Diamond[b___,Times[-1,a__],c___]:=-Diamond[b,Times[a],c]
Diamond[a___,Times[-1,b__]]:=-Diamond[a,Times[b]]
Diamond[Times[s_,a_J],b___]:=s Diamond[Times[a],b]
Diamond[b___,Times[s_,a_J],c___]:=s Diamond[b,Times[a],c]
Diamond[a___,Times[s_,b_J]]:=s Diamond[a,Times[b]]
Diamond[b___,Times[a_p,s_],c___]:=s Diamond[b,Times[a],c]
Diamond[Times[a_p,s_],b___]:=s Diamond[Times[a],b]
Diamond[a___,Times[b_p,s_]]:=s Diamond[a,Times[b]]
Diamond[b___,Times[a_Diamond,s_],c___]:=s Diamond[b,Times[a],c]
Diamond[Times[a_Diamond,s_],b___]:=s Diamond[Times[a],b]
Diamond[a___,Times[b_Diamond,s_]]:=s Diamond[a,Times[b]]


(*Generating the binary product*)


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


(*Programms on the closed form and convolution map*)


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


KLTR[f_List]:=KLTR[f]=Which[Length[f]>2,Module[{imax=Max[f],idmax,fr,fl},idmax=Position[f,imax]//Flatten;fr=Drop[f,idmax];
fl=f[[1;;(idmax[[1]]-1)]];2dot[p[imax],(p/@fl)//Total]KLTR[fr]],Length[f]==2,2dot@@(p/@f),Length[f]<2,Print["Length of f should be large than 2"]]
KLTRM[n_]:=KLTR/@(Drop[#,-1]&/@ddmBasis[Range[n]])


Tt[od_]:=Module[{phat=Range[Length@od],leftv,num,den},phat[[1]]=v;Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[od[[1;;(i-1)]]],#<od[[i,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,2,Length@od}];leftv=phat;
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],v]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=(-1)^Length@od dot[v,p[1]]Product[dot[v,p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
num/den
]
Tp[{i_}]:=dot[\[Epsilon][i],v]
Tp[f__]:=Module[{od={f},phat,leftv,num,den},phat=Range[Length@od];
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
ET2F[f_ET]:=Module[{flavor=f[[2]],od=f[[1,2]],cod=f[[1,1]],phat,leftv,rightv,num,den,scalars},phat=Range[Length@od];
scalars=Complement[Flatten[cod],Flatten[od]];
(*If[od===T[{i_}]];*)
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)phat[[i]]=p@@Select[Flatten[Join[scalars,od[[1;;(i-1)]]]],Position[cod,#][[1,1]]<Position[cod,od[[i,1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,1,Length@od}];leftv=phat;
rightv=Range[Length@od];
Table[(*If[od[[i,1]]>Max[Flatten[od\[LeftDoubleBracket]1;;(i-1)\[RightDoubleBracket]]],*)rightv[[i]]=p@@Select[Flatten[Join[scalars,od[[1;;(i-1)]]]],Position[cod,#][[1,1]]>Position[cod,od[[i,-1]]][[1,1]]&](*,phat[[i]]=p@@Range[od[[i,1]]-1]]*),{i,1,Length@od}];
num=Times@@Table[dot[leftv[[i]],F/@od[[i]],rightv[[i]]]/.dot[gg_,List[gf___],hh_]:> dot[gg,gf,hh],{i,Length@od}];
den=1/2 dot[p@@scalars]Product[1/2 dot[(p@@scalars)+p@@Flatten[od[[1;;(i-1)]]]],{i,2,Length@od}];
den=den/.dot[gg_]:>dot[gg]-m^2;
flavor=CenterDot@@(flavor/.a[i_]:>t^a[i]);
 flavor num/den
]
ET2F2s[f_ET]:=Module[{flavor=f[[2]],od=f[[1,2]],cod=f[[1,1]],phat,leftv,rightv,num,den,refs,sc,refg},leftv=Range[Length@od];
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
den=den/.dot[gg_]:>dot[gg]-m^2;
flavor=CenterDot@@(flavor/.a[i_]:>t^a[i]);
 flavor num/den
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


(*Programs on the Hopf algebra*)


shuffle[s1_List,s2_List]:=shuffle[s1,s2]=Module[{p,tp,ord},p=Permutations@Join[1&/@s1,0&/@s2]\[Transpose];
tp=BitXor[p,1];
ord=Accumulate[p] p+(Accumulate[tp]+Length[s1]) tp;
Outer[Part,{Join[s1,s2]},ord,1][[1]]\[Transpose]]


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


\[ScriptCapitalK][i_,0]:=ET[GT[{i},{}],{a[i]}]
\[ScriptCapitalK][i_,1]:=ET[GT[{i},{{i}}],{}]


\[CapitalOmega][f_,g_]:=\[FivePointedStar][f,g]-\[FivePointedStar][g,f]
\[CapitalOmega][f1_,f2_,g__]:=\[CapitalOmega][\[CapitalOmega][f1,f2],g]


S[f_T]:=If[Length[f]>1,-f-Sum[\[FivePointedStar][S[f[[1;;i]]],f[[i+1;;Length[f]]]],{i,1,Length[f]-1}],-f]/.T[gg___]:>T@@(Sort/@{gg})
S[f_Plus]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> S
S[f_Times]:=(NonCommutativeMultiply[f]//ExpandNCM)/.NonCommutativeMultiply-> S
S[\[DoubleStruckCapitalI]]:=\[DoubleStruckCapitalI]


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
