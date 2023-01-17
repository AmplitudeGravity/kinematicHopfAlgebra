(* MultivariateApart: Generalized Partial Fractions
 *
 * By Matthias Heller, Johannes Gutenberg University Mainz,
 * and Andreas von Manteuffel, Michigan State University.
 *
 * Original release: January 2021
 *
 * This software is provided for academic purposes, but without any warranty;  without
 * even the implied warranty of merchantability or fitness for a particular purpose.
 * If you use this package for your research, please cite our paper.
 *)

(* ::Package:: *)

MultivariateApart`Private`Version = "2021-01-18";
BeginPackage["MultivariateApart`"]

If[$Notebooks,
  CellPrint[Cell[#,"Print",CellFrame->0.5,FontColor->RGBColor[0.09375,0.26953125,0.23046875]]]&,
  Print]["\
MultivariateApart -- Multivariate partial fractions. \
By Matthias Heller (maheller@students.uni-mainz.de) and Andreas von Manteuffel (vmante@msu.edu). \
Version "<>MultivariateApart`Private`Version<>". \
Please see ?MultivariateApart and ?MultivariateApart`* for help."]

(*** options: ***)

MaxDegree::usage="\
MaxDegree is an option for the function MultivariateApart. Denominator factors with total degree \
larger than this value will be left untouched in the partial fraction decomposition.";

Verbose::usage="\
Verbose is an option for different functions, which controls the display of auxiliary information.";

PreTransformation::usage="\
PreTransformation is an option for the function MultivariateApart, which specifies the function \
to be applied initially. Please note the default call of Together ensures a unique result, \
independ of the initial representation.";

InverseSymbol::usage="\
InverseSymbol is an option to the function AbbreviateDenominators. It specifies a string, which is \
used as a prefix to generate symbols representing inverse denominator factors."

PartitionSize::usage="\
PartitionSize is an option to MultivariateApart, MultivariateAbbreviatedApart and ApartReduce, \
which partitions the terms into smaller pieces. Using a number like 1000 can sometimes help \
to workaround a possible performance degradation in Mathematica for a large number of terms.";

Iterate::usage="\
Iterate is an option to MultivariateApart, MultivariateAbbreviatedApart and ApartReduce. \
For a factorized input expression, the option Iterate->True can be used to iteratively take \
factors into account. This can lead to a reduction in runtime if many cancellations between \
different terms of a sum occur.";

UseFormProgram::usage="\
UseFormProgram Iterate is an option to ApartReduce. If enabled, the external symbolic manipulation \
program Form will be used. Please make sure to check the options of Write"

(*** the all-in-one functions: ***)

MultivariateApart::struct="Encountered unknown structure, can't handle this expression.";
MultivariateApart::usage="\
MultivariateApart[expr] calculates a partial fraction decomposition of a rational expression of several variables. \
Several options are available to improve the performance of the function, please see the documentation for \
the entries in Options[MultivarateApart] for further information.";

MultivariateAbbreviatedApart::usage="\
MultivariateAbbreviatedApart[expr] has the same syntax as MultivariateApart[expr], but returns a list {res,abbr} where \
the result res is given in terms of abbreviated denominator factors and abbr defines the abbreviations. ";

(*** step-by-step approach within Mathematica: ***)

AbbreviateDenominators::collision="A new symbol in `1` clashes with an existing symbol. Please change option InverseSymbol.";
AbbreviateDenominators::usage="\
AbbreviateDenominators[expr] analyzes the denominators of a rational function and abbreviates \
the denominator factors by symbols. The argument should be a ratio of two polynomials or a sum \
of such expressions. The return value is {exprsym, denfac, invsym}, where exprsym is the original \
expression with denominators replaced by symbols and denfac is a list of denominator factors. \
The symbols invsym represent the inverses of the denominator factors such that invsym == 1/denfac. \
The function can be called with a second argument which is a list of denominator factors that are in \
use already. In that way the function will replace denominators in the corresponding order.";

ApartOrder::usage="\
ApartOrder[denfac,invsym,othersym] generates an optimized monomial ordering for a partial fraction \
decomposition. The argument dens is a list of denominator factors, invsym is a list of symbols \
denoting the inverse powers of the factors in den. The optional argument othersym is a list of \
additional symbols which occur as numerator factors."

ApartBasis::usage="\
ApartBasis[denfac,invsym,ord] calculates a Groebner basis for a partial fraction decomposition. \
The list denfac specifies the denominator factors and invsym symbols for their inverses. \
The argument ord is a monomial ordering specified in the format as returned by ApartOrder.";

ApartReduce::usage="\
ApartReduce[expr,basis,ord] performs a partial fraction decomposition of the expression expr. \
The argument basis is the Groebner basis as calculated by ApartBasis and ord is the monomial \
ordering as calculated by ApartOrder. Several options are available for performance tuning, \
please see the documentation of the entries in Options[ApartReduce].";

Options[ApartReduce] = { PartitionSize->Infinity, Iterate->False, UseFormProgram->False };
Options[AbbreviateDenominators]={InverseSymbol->"q", Verbose->False};
Options[MultivariateApart]=Join[{MaxDegree->Infinity, Verbose->False, PreTransformation->Together, PartitionSize->Infinity, Iterate->False, UseFormProgram->False },
  Options[AbbreviateDenominators]];
Options[MultivariateAbbreviatedApart]=Options[MultivariateApart];

(*** interface to Groebner basis calculation with Singular or rule applications with Form ***)

WriteSingularBasisInput::usage="\
WriteSingularBasisInput[denfac,invsym,ord,dir] generates a Singular input file \"apartbasisin.sing\" in the directory dir to replace the call to ApartBasis by an external computation. \
This file can be executed in Singular with the command execute(read(\"apartbasisin.sing\")). Singular calculates the Groebner basis and writes it to the \
file \"apartbasisout.m\". That basis can be read into Mathematica with Get[\"apartbasisout.m\"] and used in the same way as calculated with ApartBasis[denfac,invsym,ord]";

WriteFormReduceInput::usage="\
WriteFormReduceInput[gb,ord,\"dir\"] generates three files called \"multivariateapart.h\", \"apartpattern.inc\" and \"apartrules.inc\" for the Form symbolic manipulation program. \
The first file contains a Form procedure, that implements partial fractioning in Form. It can be included in a Form script via \"#include multivariateapart.h\" in the beginning of the script. \
Note that the auxillary files \"apartrules.inc\" and \"apartpattern.inc\" have to be in the same directory. The Form procedure can be invoked by \"#call partialfraction(expr)\", \
where expr is a local expression in Form.
";

(*** auxiliary functions ***)

CollectFlat::usage = "\
CollectFlat[expr,x] returns Collect[expr,x] with nested powers of x \
expanded.
CollectFlat[expr,{x1,x2,...}] and CollectFlat[expr,var,h] perform the \
same for the respective Collect[] invocations.
CollectFlat[expr,x,hc,hx] applies hc to the coefficients and hx to the \
powers products.";

(* at least the message should be public: *)
FindRawDenominators::struct="Not a sum of simple rational functions, can't handle this expression.";
FindRawDenominators::usage="\
FindRawDenominators[expr,vars] scans for all denominators in a rational function written over a single denominator or \
in a sum of such expressions. The optional argument vars specifies the list of indeterminants."

(* should be public for further applications like MonomialList etc. *)
DegRevLexWeights::usage="\
DegRevLexWeights[{{x1,x2,..},{y1,y2,..},..}] returns the weight matrix corresponding to a block ordering, \
where the ordering within a block is degree reverse lexicographic."

TotalDegree::usage="\
TotalDegree[poly] returns the total degree of a polynomial poly in all of its variables."

$ApartTemporaryDirectory::usage="\
$ApartTemporaryDirectory sets the directory name to be used for temporary files when running external programs like Form."



(*** some documentation for functions which are currently private ***)

(*

FactorDenominators::usage="\
FactorDenominators[dens] factorizes a list of polynomials.
The return value is a two-component list, where the first component is the list of factors and the second \
component gives the factorization of the input polynomials."

OrderPolynomials::usage="\
OrderPolynomials[{poly1,poly2,..}] orders a list of polynomials into groups according to the appearance of variables \
and sorts within each group according to total degree."

*)



Begin["`Private`"]

AbbreviateDenominators[expr_,den0_List:{},OptionsPattern[]]:=Module[{den,fac,fac0,inv,den2sym,fac2sym,res,sym,new},
  Check[den=FindRawDenominators[expr],Return[$Failed,Module]];
  {fac, den2sym}=FactorDenominators[den];
  fac0=FactorDenominators[den0][[1]];
  If[Complement[den0,fac0]=!={}||Complement[fac0,den0]=!={},
    Print["List of denominators is reducible. Please normalize with FactorDenominators[]."];
    Return[$Failed,Module]];
  new=Complement[fac,den0];
  If[den0!={}&&new=!={},Print["found ",Length@new," new denominator(s): ",InputForm[new]]];
  fac=Join[den0,new]; (* never change user specified den0 defs ! *)
  sym=ToString[OptionValue[InverseSymbol]];
  inv=ToExpression[sym<>ToString[#]]&/@Range[Length[fac]];
  If[!FreeQ[expr,Alternatives@@inv],Message[AbbreviateDenominator::collision,inv];Return[$Failed,Module]];
  fac2sym=Dispatch[Thread[fac->(1/inv)]];
  den2sym=Dispatch@Table[r[[1]]->Times@@(Replace[#[[1]],fac2sym]^#[[2]]&/@r[[2]]),{r,den2sym}];
  res=expr /. Power[x_,n_] /; (n<0) :> Replace[x,den2sym]^n;
  {res,fac,inv}
];

Clear[MultivariateApart];
MultivariateApart[expr_,opts:OptionsPattern[]]:=Module[{res=MultivariateAbbreviatedApart[expr,{},opts,InverseSymbol->Unique["q"]]},
  res[[1]]/.res[[2]]];

Clear[MultivariateAbbreviatedApart];
MultivariateAbbreviatedApart[expr_,den0_List:{},opts:OptionsPattern[]]:=Module[{nf,fac,bnd,inv,allfac0,allinv0,pos,allfac,allinv,fac2sym,sym2fac,gb,ord,verb,pre,res},
  verb=OptionValue[Verbose];
  nf=OptionValue[PreTransformation]@expr; (* could also be used without Together, but not unique in that case *)
  If[verb,Print["applied pre transformation"]];
  {nf,allfac0,allinv0}=AbbreviateDenominators[nf,den0,FilterRules[Flatten[{opts}],Options[AbbreviateDenominators]]];
  (* select only qis that appear in expr *)
  allinv=Select[allinv0,MemberQ[{nf},#,Infinity]&];
  pos=Table[Position[allinv0,allinv[[i]]],{i,1,Length@allinv}]//Flatten;
  allfac=allfac0[[pos]];
  If[verb, Print["all invariants: ", allinv]];
  {fac,inv}=If[#==={},{{},{}},Transpose[#]]&@Select[Transpose@{allfac,allinv}, (TotalDegree[#[[1]]]<=OptionValue[MaxDegree])&];
  ord=ApartOrder[fac, inv, Sort@Variables@nf];
  If[!PolynomialQ[nf,Flatten@ord],Message[MultivariateApart::struct];Return[$Failed,Module]];
  sym2fac=Dispatch[Thread[allinv->(1/allfac)]]; (* before filtering *)
  If[verb, Print["all factors: ",Normal@sym2fac, "\nordering:  ",ord]];
  If[fac==={}, Return[{CollectFlat[nf, allinv], sym2fac}], Module]; (* avoid errors in GroebnerBasis, e.g. for pure numbers *)
  gb=ApartBasis[fac, inv, ord];
  res=ApartReduce[nf, gb, ord, FilterRules[Flatten[{opts}],Options[ApartReduce]]];
  {CollectFlat[res, allinv], sym2fac}
];

Clear[ApartOrder];
ApartOrder[dens_List,inv_List,othervars_:{}]:=Module[{dengrps},
  dengrps=OrderPolynomials[dens, inv];
  Append[dengrps,ComplementPreserveOrder[DeleteDuplicates[Join[othervars,Variables[dens]]],inv]]
]

ApartBasis[fac_, inv_, ord_]:=GroebnerBasis[PartialFractionGenerators[fac,inv], ord//Flatten, MonomialOrder->DegRevLexWeights[ord], Method->"Buchberger"];


$ApartTemporaryDirectory=FileNameJoin[{$TemporaryDirectory,"multivariateapart"}]

CreateDirectoryIfNeeded[dirname_String]:=Switch[FileType[dirname],
  None,Return[CreateDirectory[dirname,CreateIntermediateDirectories->True]] (*create dir*),
  Directory,Null;Return[dirname] (*do nothing*),
  File,Print["ERROR: File with name chosen for the directory already exists !";Return[$Failed]] (*error!*)
]

Clear[WriteSingularBasisInput];
WriteSingularBasisInput[fac_, inv_, ord_, dir_]:=Module[{str,infile,outfile},
  If[CreateDirectoryIfNeeded[dir]==$Failed,Return[$Failed,Module]];
  infile=FileNameJoin[{dir,"apartbasisin.sing"}];
  outfile=FileNameJoin[{dir,"apartbasisout.m"}];
  str=OpenWrite[infile];
  WriteLine[str, "ring myring=0,"];
  WriteLine[str, "("<>Riffle[ToString[CForm[#]]&/@Flatten[ord],","]<>"),"];
  WriteLine[str, "("<>Riffle[("dp("<>ToString[Length[#]]<>")")&/@ord,","]<>");"];
  Do[WriteLine[str, "poly mygen"<>ToString[i]<>"=1-"<>ToString[InputForm[fac[[i]]*inv[[i]]]]<>";"],{i,Length[fac]}];
  WriteLine[str, "ideal myideal="<>Riffle[Table["mygen"<>ToString[i],{i,Length[fac]}],","]<>";"];
  WriteLine[str,"ideal mygb = slimgb(myideal);"];
  WriteLine[str,"link outfile = \"ASCII:w "<>outfile<>"\";"];
  WriteLine[str,"write(outfile, \"{\", mygb, \"}\");"];
  WriteLine[str,"close(outfile);"];
  WriteLine[str,"quit;"];
  Close[str];
]

Clear[formtemplate];
formtemplate="\n
#procedure apartreduce(expr)
  .sort
  CF fac;
  S x;
  Skip;
  NSkip `expr';
  AntiBracket dens;
  .sort
  Skip;
  NSkip `expr';

  Collect fac;
  factarg fac;
  id fac(x?number_,?a) = x*fac(?a);
*  repeat id fac(?a,x?,x?,?b) = fac(?a,x*x,?b);
  .sort
  Skip;
  NSkip `expr';
  $repeat = 0;
  #do i=1,1
   id once fac(x?,?a) = x*fac(?a);
   id fac = 1;
    #do k=1,1
     #include apartpattern.inc
     #write \"%$\",$repeat
     #include apartrules.inc

     if ($repeat == 1);
      redefine k \"0\";
      $repeat = 0;
     endif;
     .sort
     Skip;
     NSkip `expr';
    #enddo
  if (match(fac(?a)));
   redefine i \"0\";
  endif;
  .sort
  Skip;
  NSkip `expr';
  #enddo

#endprocedure
";

Clear[WriteFormReduceInput];
WriteFormReduceInput[gb_,ord_,formdirectory_]:=Module[{monlist,formrules,patternmatch,pre,str,qstring,varstring,template},
  If[CreateDirectoryIfNeeded[formdirectory]==$Failed,Return[$Failed,Module]];
  monlist=MonomialList[SortBy[gb,Length],Flatten@ord,DegRevLexWeights[ord]];
  patternmatch=MonomialList[SortBy[gb,Length],Flatten@ord,DegRevLexWeights[ord]];
  formrules=StringJoin@@Table[
    pre=FactorList[monlist[[i,1]]][[1]];
    pre=pre[[1]]^pre[[2]];
    "#if `rule"<>ToString[i]<>"'==1\n"<>
    "#write \"calc module nr "<>ToString[i]<>"\"\n#do j=1,1\n"<>
    "if (match("<>ToString@InputForm@Factor@(monlist[[i,1]]/pre)<>"));\n"<>
    "id "<>ToString@InputForm@Factor@(monlist[[i,1]]/pre)<>" = "<>ToString@InputForm@(Plus@@Factor/@(-monlist[[i,2;;All]]/pre))<>";\n"<>
    "redefine j \"0\";\nendif;\n.sort\nSkip;\nNSkip `expr';\n#enddo\n#endif\n"
  ,{i,1,Length@monlist}];

  str=OpenWrite[formdirectory<>"/apartrules.inc"];
  WriteString[str,formrules];
  Close[str];  
  
  patternmatch=StringJoin@@Table[
    pre=FactorList[patternmatch[[i,1]]][[1]];
    pre=pre[[1]]^pre[[2]];
    "if (match("<>ToString@InputForm@Factor@(monlist[[i,1]]/pre)<>"));\n"<>
    " redefine rule"<>ToString[i]<>" \"1\";\n"<>
    " $repeat = 1;\nendif;\n"
  ,{i,1,Length@monlist}];
  
  str=OpenWrite[formdirectory<>"/apartpattern.inc"];
  Do[WriteLine[str,"#define rule"<>ToString[i]<>" \"0\""],{i,1,Length@monlist}];
  WriteString[str,patternmatch];
  WriteLine[str,".sort\nSkip;\nNSkip `expr';"];
  Close[str];
  
  qstring=Reverse@Flatten@ord[[1;;-2]];
  varstring=ord[[-1]];
  qstring=StringDrop[StringJoin@@Table[ToString[qstring[[i]]]<>",",{i,1,Length@qstring}],-1];  
  varstring=StringDrop[StringJoin@@Table[ToString[varstring[[i]]]<>",",{i,1,Length@varstring}],-1];  

  str=OpenWrite[formdirectory<>"/apartreduce.h"];
  WriteLine[str,"S I;"]; (* imaginary unit *)
  WriteLine[str,"S "<>qstring<>";"];
  WriteLine[str,"S "<>varstring<>";"];
  WriteLine[str,"set dens: "<>qstring<>";"];
  WriteString[str,formtemplate];
  Close[str];
]


PlusToList[expr_Plus]:=List@@expr;
PlusToList[expr_]:={expr};

TimesToList[expr_Times]:=List@@expr;
TimesToList[expr_]:={expr};

Clear[formreduce];
formreduce="
#:workspace 100M
*Off Statistics;
#:TermsInSmall 500000
#:SmallSize 1000M
#:largesize 1000M
#:smallextension 2000M
#:scratchsize 1000M
#include apartreduce.h
#include apartreducein.inc
#call apartreduce(M);
.sort
#write <apartreduceout.m> \"%e\",M;
.end"

Clear[ApartReduce];
ApartReduce[expr_List, gb_, ord_, opts:OptionsPattern[]]:=ApartReduce[#,gb,ord,opts]&/@expr;
ApartReduce[expr_, gb_, ord_, opts:OptionsPattern[]]:=Module[{n, p, var=Flatten[ord], t, p1, p2, p3, ordm,str,ret,outfile,dirname,runfilename},
  If[NumberQ[expr],Return[expr]];
  If[OptionValue[UseFormProgram],
    dirname=$ApartTemporaryDirectory;
    WriteFormReduceInput[gb,ord,dirname];
    str=OpenWrite[FileNameJoin[{dirname,"apartreducein.inc"}]];
    WriteLine[str,"L M="];
    WriteString[str,(ToString@InputForm@expr)<>";"]; 
    Close[str];
    runfilename="apartreducerun.frm";
    str=OpenWrite[FileNameJoin[{dirname,runfilename}]];
    WriteLine[str,formreduce];
    Close[str];
    outfile=FileNameJoin[{dirname,"apartreduceout.m"}];
    If[FileExistsQ[outfile],DeleteFile[outfile]];
    ret=RunProcess[{"form", runfilename}, ProcessDirectory->dirname];
    If[("ExitCode"/.ret)=!=0,Print["ERROR: running form failed"];Print["StandardOutput"/.ret];Print["StandardError"/.ret];Return[$Failed]];
    Return@ToExpression@StringDelete[Import[outfile,"String"],{";"," ","\n"}];
  ];
  If[OptionValue[Iterate],
    p=TimesToList[expr];
    p={#,Which[FreeQ[#,Alternatives@@var],1,MatchQ[#,d_^_./;MemberQ[var,d]],2,True,3]}&/@p;
    (* p1=numerical factors, p2=denominators, p3=numerators *)
    {p1,p2,p3}={Cases[p,{t_,1}:>t],Cases[p,{d_^n_.,2}:>{d,n}],Cases[p,{t_,3}:>t]};
    If[Length[p1]+Length[p2]+Length[p3]=!=Length[p],Message[MultivariateApart::struct];Return[$Failed,Module]];
    (*p=ApartReduce[Times@@p3,gb,ord,Iterate->False,opts];*)
    p=Times@@p3;
    (* now reduce denominator factors one at a time *)
    p2=#[[1]]^#[[2]]&/@Reverse[SortBy[p2, Position[var, #[[1]]][[1]]&]];
	(*Monitor[*)
	Do[p=Expand@ApartReduce[p*p2[[t]], gb, ord, Iterate->False, opts], {t, Length@p2}];
    (*, {ToString[t]<>"/"<>ToString[Length[p2]],ByteCount[p],p2[[t]], Drop[p2,t]}];*)
	Return[(Times@@p1)*p, Module];
  ];
  n=OptionValue[PartitionSize];
  p=If[n===Infinity, {expr}, Plus@@@Partition[PlusToList@Expand@expr,n,n,{1,1},{}]];
  ordm=DegRevLexWeights[ord];
  (*Monitor[*)
  p=Table[PolynomialReduce[p[[t]],gb,var,MonomialOrder->ordm][[-1]],{t,Length@p}];
  (*, ToString[t]<>"/"<>ToString[Length[p]]];*)
  Plus@@p
]

Clear[CollectFlat];
CollectFlat[e_, t_, hc_: Identity, ht_: Identity] := Module[{h},
 Expand@Collect[e, t, h] /. h[x_] y_. :> hc[x] ht[y]]

Clear[ComplementPreserveOrder];
ComplementPreserveOrder[l_List, othersets__]:=l[[Sort[Complement[l, othersets] /. Thread[l -> Range@Length@l]]]]

TotalDegree[poly_]:=Max[Plus@@@(First/@CoefficientRules[poly])]


Clear[RawFactors];
RawFactors[expr_Times]:=Flatten[RawFactors/@(List@@expr)];
RawFactors[Power[x_,n_?IntegerQ]]:=RawFactors[x];
RawFactors[_?NumericQ]:={};
RawFactors[x_]:={x};

Clear[FindRawDenominators];
FindRawDenominators[expr_]:=FindRawDenominators[expr, Variables[expr]];
FindRawDenominators[expr_List,vars_List]:=Union@Flatten[FindRawDenominators[#,vars]&/@expr];
FindRawDenominators[expr_Plus,vars_List]:=FindRawDenominators[List@@expr,vars];
FindRawDenominators[expr_,vars_List]:=Module[{nd},
  nd=If[TrueQ[$VersionNumber >= 12.0], NumeratorDenominator[expr],{Numerator[expr],Denominator[expr]}];
  If[!And@@(PolynomialQ[#,vars]&/@nd),
    Message[FindRawDenominators::struct]; Return[$Failed,Module]];
  RawFactors[nd[[2]]]
];

Clear[FactorDenominators];
FactorDenominators[dens_List]:=Module[{faclst},
  faclst=FactorList/@dens;
  {Union@Cases[First/@(Join@@faclst), Except[_?NumericQ]], Thread[dens->faclst]}];

OrderPolynomials[{},_]={};
OrderPolynomials[polys_List,syms_]:=Module[{comb=Transpose[{polys,syms}]},
  (* we combine polys (explicit polynomials) and syms (just a list of symbols) to two-element lists,
     sort according to the explicit poly, but return the symbol for it *)
  Map[Last,Reverse@SortBy[
    Reverse/@SortBy[TotalDegree[#[[1]]]&]/@GatherBy[comb,Sort@Variables[#[[1]]]&],
    Length@Variables[#[[1]]]&], {2}]
]

Clear[DegRevLexWeights];
DegRevLexWeights[varblocks_]:=Module[{l=Length@Flatten@varblocks,n=1},
  Join@@Table[Join[{Table[If[Length@Flatten@varblocks[[1;;i-1]]<k<=Length@Flatten@varblocks[[1;;i]],1,0],{k,1,l}]},
         Table[Table[If[k==Length@Flatten@varblocks[[1;;i]]-j+1,-1,0],{k,1,l}],{j,1,Length[varblocks[[i]]]-1}]],
  {i,1,Length[varblocks]}] (*//SparseArray would be nice but conflicts with Singular interface*)
];

PartialFractionGenerators[dens_,inv_]:=1-dens*inv;

End[]

EndPackage[]

