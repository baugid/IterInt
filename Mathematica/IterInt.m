(* ::Package:: *)

Print["****************** IterInt v1.0 ******************\n* Authors: Gideon Baur, Claude Duhr              *\n* Email: gbaur@uni-bonn.de, cduhr@uni-bonn.de    *\n* arXiv: 2606.02744                              *\n**************************************************"];


BeginPackage["IterInt`"]

(*regulate::usage="regulate[{f, \!\(\*SubscriptBox[\(Res\), \(0\)]\)[f]},\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\), \!\(\*SubscriptBox[\(z\), \(max\)]\), v] gives a regulated expression for the iterated integral \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] f\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)] with regulator v in terms of integral objects.";*)
(*execIntegration::usage="execIntegration[f,\[Ellipsis], \!\(\*SubscriptBox[\(z\), \(max\)]\)] performs the convergent iterated integral \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] f\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)]";*)
plainIntegrate::usage="plainIntegrate is a value for the integrationMethod option. It performs both the regularization and integration of each occuring iterated integral independently.";
(*"plainIntegrate[{f, \!\(\*SubscriptBox[\(Res\), \(0\)]\)[f]},\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\), \!\(\*SubscriptBox[\(z\), \(max\)]\)] performs both the regularization and the integration of \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] f\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)].";*)
splittingPlainIntegrate::usage="splittingPlainIntegrate is a value for the integrationMethod option. It performs all iterated integrals independently, but only evaluates the regulated expressions along a short segment.";
(*"splittingPlainIntegrate[{f, \!\(\*SubscriptBox[\(Res\), \(0\)]\)[f]},\!\(\*
StyleBox[\"\[Ellipsis]\", \"TR\"]\), \!\(\*SubscriptBox[\(z\), \(max\)]\)] performs both the regularization and the integration of \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] f\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)] using a split approach.";*)
treeIntegrate::usage="treeIntegrate is the default value for the integrationMethod option. It combines all appearing iterated integrals into trees and evaluates them together.";
(*"treeIntegrate[{{f, \!\(\*SubscriptBox[\(Res\), \(0\)]\)[f]},\[Ellipsis]},{{\!\(\*SubscriptBox[\(i\), \(1\)]\), \!\(\*SubscriptBox[\(i\), \(2\)]\),\[Ellipsis]}, \[Ellipsis]}, \!\(\*SubscriptBox[\(z\), \(max\)]\)] performs the regularization and integration of all iterated integrals {\!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] \*SubscriptBox[\(f\), SubscriptBox[\(i\), \(1\)]]\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)],\[Ellipsis]} at once.";*)
splittingTreeIntegrate::usage="splittingTreeIntegrate is a value for the integrationMethod option. It combines the appearing iterated integrals into trees, but only evaluates the regulated expressions along a short segment. This combines essentially treeIntegrate with splittingPlainIntegrate.";
(*"splittingTreeIntegrate[{{f, \!\(\*SubscriptBox[\(Res\), \(0\)]\)[f]},\[Ellipsis]},{{\!\(\*SubscriptBox[\(i\), \(1\)]\), \!\(\*SubscriptBox[\(i\), \(2\)]\),\[Ellipsis]}, \[Ellipsis]}, \!\(\*SubscriptBox[\(z\), \(max\)]\)] performs the regularization and integration of all iterated integrals {\!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] \*SubscriptBox[\(f\), SubscriptBox[\(i\), \(1\)]]\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)],\[Ellipsis]} at once using a split approach.";*)
integral::usage="integral[f,\[Ellipsis], \!\(\*SubscriptBox[\(z\), \(max\)]\)] describes the convergent iterated integral \!\(\*FormBox[\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), SubscriptBox[\(z\), \(max\)]]\[DifferentialD]\*TemplateBox[{\"n\"},\n\"C\"] \[Ellipsis] \(\*SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\n\"C\"]]\[DifferentialD]\*TemplateBox[{\"1\"},\n\"C\"] f\)\),
TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\n\"C\"]\)]. This head is only used for checkFunction.";

regulated::usage="regulated[name] denotes the integration kernel \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\) with the pole removed.";
pole::usage="pole[name] denotes the pole term of the integration kernel \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\).";
regulateExpression::usage="regulateExpression[expr] regulates all IIntegral objects in \!\(\*StyleBox[\"expr\",FontSlant->\"Italic\"]\) using regulator 1.\nregulateExpression[expr,regulator] performs the regulation with the given regulator.";
kernelResidue::usage="kernelResidue[name] gives the residue of the integration kernel \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\) in zero.";
kernelExpression::usage="kernelExpression[name] gives the defining functions for the integration kernel \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\).";
defineKernel::usage="defineKernel[name, expr] adds the integration kernel \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\) with callable expression \!\(\*StyleBox[\"expr\",FontSlant->\"Italic\"]\)\!\(\*StyleBox[\".\",FontSlant->\"Italic\"]\)\ndefineKernel[name, expr, residue] also specifies the residue in zero.\ndefineKernel[assoc] adds the kernel specified by an association with at least the keys \!\(\*StyleBox[\"name\",FontSlant->\"Italic\"]\)\!\(\*StyleBox[\" \",FontSlant->\"Italic\"]\)and \!\(\*StyleBox[\"expression\",FontSlant->\"Italic\"]\)\!\(\*StyleBox[\".\",FontSlant->\"Italic\"]\)\!\(\*StyleBox[\" \",FontSlant->\"Italic\"]\)The key \!\(\*StyleBox[\"residue\",FontSlant->\"Italic\"]\)\!\(\*StyleBox[\" \",FontSlant->\"Italic\"]\)is optional.";
listKernels::usage="listKernels[] returns a list of all defined integration kernels.";
exportKernel::usage="exportKernel[name] returns an association describing the kernel \!\(\*
StyleBox[\"name\",\nFontSlant->\"Italic\"]\).";
exportKernels::usage="exportKernels[] returns a list of associations describing all defined kernels.";
deleteKernel::usage="deleteKernel[name] deletes the kernel given by \!\(\*
StyleBox[\"name\",\nFontSlant->\"Italic\"]\).";
name::usage="Key for integration kernels represented by associations. See defineKernel.";
expression::usage="Key for integration kernels represented by associations. See defineKernel.";
residue::usage="Key for integration kernels represented by associations. See defineKernel.";
IIntegral::usage="IIntegral[\!\(\*SubscriptBox[\(f\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(z\), \(max\)]\)] represents the iterated integral given by \!\(\*FormBox[RowBox[{SubsuperscriptBox[\"\[Integral]\", \"0\", SubscriptBox[\"z\", \"max\"]], RowBox[{RowBox[{\"\[DifferentialD]\", TemplateBox[{\"n\"},\"C\"]}], \"\[Ellipsis]\", RowBox[{SubsuperscriptBox[\"\[Integral]\", \"0\", TemplateBox[{\"2\"},\"C\"]], RowBox[{RowBox[{\"\[DifferentialD]\", TemplateBox[{\"1\"},\"C\"]}], SubscriptBox[\"f\", \"1\"]}]}]}]}],TraditionalForm]\)[\!\(\*TemplateBox[{\"1\"},\"C\"]\)].";
computeIIntegrals::usage="computeIIntegrals[expr] numerically evaluates all IIntegral objects in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).\ncomputeIIntegrals[expr, rules] numerically evaluates all IIntegral objects in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\) and beforehand performs the replacements in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\) specified by \!\(\*
StyleBox[\"rules\",\nFontSlant->\"Italic\"]\).";

(*Options*)
startRegulator::usage="startRegulator is an option to the integration methods moving the starting point away from zero. Possible values are between zero and one.";
regulator::usage="regulator is used to change the tangent used for the tangential base-point regularization. It can take any complex value except zero.";
checkFunction::usage="checkFunction is used by the integration functions to quickly evaluate some integrals by the user.";
switchPoint::usage="switchPoint is used by the splitting algorithm to determine where to transition between the different methods. Valid values are between zero and one.";
integrationAlgorithm::usage="integrationAlgorithm is an option value to choose the integration algorithm to use.";

Global`$NativeIteratedIntegral::usage="If set to True before loading this library, native functions for iterated integrals will be provided.";

If[Global`$NativeIteratedIntegral=!=True, Begin["Ununsed`"]];
plainIntegrateC::usage="plainIntegrateC is, up to different options, a drop in replacement for plainIntegrate using the C++/GSL backend.";
splittingPlainIntegrateC::usage="splittingPlainIntegrateC is, up to different options, a drop in replacement for splittingPlainIntegrateC using the C++/GSL backend.";
treeIntegrateC::usage="treeIntegrateC is, up to different options, a drop in replacement for treeIntegrateC using the C++/GSL backend.";
splittingTreeIntegrateC::usage="splittingTreeIntegrateC is, up to different options, a drop in replacement for splittingTreeIntegrateC using the C++/GSL backend.";
computeIIntegralsC::usage="computeIIntegralsC is, up to different options, a replacement for computeIIntegrals, which uses the C++/GSL backend.";
epsAbs::usage="epsAbs is an option for native integration functions specifying a target on the absolute deviation of the result.";
epsRel::usage="epsRel is an option for native integration functions specifying a target on the relative deviation of the result.";
hStart::usage="hStart is an option for native integration functions and is used to specify the initial step size of the solver.";
hMax::usage="hMax is an option for native integration functions and can be used to set an upper limit on the step size used by the solver.";
baseRange::usage="baseRange is an option for native integration functions corresponding to switchPoint. It sets the transition point between different methods.";
If[Global`$NativeIteratedIntegral=!=True, End[]];


Begin["Private`"];

If[Global`$NativeIteratedIntegral===True,
(*load the C functions*)
libraryFile = FindLibrary[DirectoryName[$InputFileName]<>"libIntegrate"];
evalIntegrF = ForeignFunctionLoad[libraryFile,"plainIntegrate",{{"CDouble","CDouble"},"OpaqueRawPointer","RawPointer"::[{"CDouble","CDouble"}],"UnsignedInteger64","CDouble","CDouble","CDouble","CDouble","CDouble",{"CDouble","CDouble"}}->{"CDouble","CDouble"}];
evalIntegrSplitF = ForeignFunctionLoad[libraryFile,"splittingPlainIntegrate",{{"CDouble","CDouble"},"OpaqueRawPointer","RawPointer"::[{"CDouble","CDouble"}],"UnsignedInteger64","CDouble","CDouble","CDouble","CDouble","CDouble","CDouble",{"CDouble","CDouble"}}->{"CDouble","CDouble"}];
treeIntegrateF = ForeignFunctionLoad[libraryFile,"treeIntegrate",{{"CDouble","CDouble"},"OpaqueRawPointer","RawPointer"::[{"CDouble","CDouble"}],"UnsignedInteger64",
"RawPointer"::["RawPointer"::["UnsignedInteger64"]],"UnsignedInteger64","CDouble","CDouble","CDouble","CDouble","CDouble",{"CDouble","CDouble"}}->"RawPointer"::[{"CDouble","CDouble"}]];
treeIntegrateSplitF = ForeignFunctionLoad[libraryFile,"splittingTreeIntegrate",{{"CDouble","CDouble"},"OpaqueRawPointer","RawPointer"::[{"CDouble","CDouble"}],"UnsignedInteger64",
"RawPointer"::["RawPointer"::["UnsignedInteger64"]],"UnsignedInteger64","CDouble","CDouble","CDouble","CDouble","CDouble","CDouble",{"CDouble","CDouble"}}->"RawPointer"::[{"CDouble","CDouble"}]];
arrayDelF = ForeignFunctionLoad[libraryFile,"deleteCmplxArray",{"RawPointer"::[{"CDouble","CDouble"}]}->"Void"];

Options[plainIntegrateC]={epsAbs->10^-12,epsRel->10^-12,hStart->1*10^-5,hMax->1,startRegulator->0,regulator->1};
plainIntegrateC[integrands__List,zVal_,OptionsPattern[]]:=Module[{functions,funcMem,residues,resMem},
functions=If[MemberQ[{ManagedObject,ForeignCallback,OpaqueRawPointer},Head[#]],#,CreateForeignCallback[val|->ReIm[#[Complex@@val]],{{"CDouble","CDouble"}}->{"CDouble","CDouble"}]]&/@{integrands}[[All,1]];
residues=ReIm[{integrands}[[All,2]]]//N;
funcMem=RawMemoryExport[functions,"OpaqueRawPointer"];
resMem=RawMemoryExport[residues,{"CDouble","CDouble"}];
Complex@@evalIntegrF[N[ReIm[zVal]],funcMem,resMem,Length[residues],OptionValue[epsAbs]//N,OptionValue[epsRel]//N,OptionValue[hStart]//N,OptionValue[hMax]//N,OptionValue[startRegulator]//N,ReIm[OptionValue[regulator]]//N]
];

Options[splittingPlainIntegrateC]={epsAbs->10^-12,epsRel->10^-12,hStart->1*10^-5,hMax->1,startRegulator->0,baseRange->2*10^-5,regulator->1};
splittingPlainIntegrateC[integrands__List,zVal_,OptionsPattern[]]:=Module[{functions,funcMem,residues,resMem},
functions=If[MemberQ[{ManagedObject,ForeignCallback,OpaqueRawPointer},Head[#]],#,CreateForeignCallback[val|->ReIm[#[Complex@@val]],{{"CDouble","CDouble"}}->{"CDouble","CDouble"}]]&/@{integrands}[[All,1]];
residues=ReIm[{integrands}[[All,2]]]//N;
funcMem=RawMemoryExport[functions,"OpaqueRawPointer"];
resMem=RawMemoryExport[residues,{"CDouble","CDouble"}];
Complex@@evalIntegrSplitF[N[ReIm[zVal]],funcMem,resMem,Length[residues],OptionValue[epsAbs]//N,OptionValue[epsRel]//N,OptionValue[hStart]//N,OptionValue[hMax]//N,OptionValue[startRegulator]//N,OptionValue[baseRange]//N,ReIm[OptionValue[regulator]]//N]
];

Options[treeIntegrateC]={epsAbs->10^-12,epsRel->10^-12,hStart->1*10^-5,hMax->1,startRegulator->0,regulator->1};
treeIntegrateC[integrands_List,integrals_List,zVal_,OptionsPattern[]]:=Module[{functions,funcMem,residues,resMem,integralLists,integralMemberMem,integralMem,res,output},
functions=If[MemberQ[{ManagedObject,ForeignCallback,OpaqueRawPointer},Head[#]],#,CreateForeignCallback[val|->ReIm[#[Complex@@val]],{{"CDouble","CDouble"}}->{"CDouble","CDouble"}]]&/@integrands[[All,1]];
residues=ReIm[integrands[[All,2]]]//N;
funcMem=RawMemoryExport[functions,"OpaqueRawPointer"];
resMem=RawMemoryExport[residues,{"CDouble","CDouble"}];
integralLists=Prepend[#-1,Length[#]]&/@integrals;
integralMemberMem=RawMemoryExport[#,"UnsignedInteger64"]&/@integralLists;
integralMem=RawMemoryExport[integralMemberMem,"RawPointer"::["UnsignedInteger64"]];
res=treeIntegrateF[N[ReIm[zVal]],funcMem,resMem,Length[residues],integralMem,Length[integrals],OptionValue[epsAbs]//N,OptionValue[epsRel]//N,OptionValue[hStart]//N,OptionValue[hMax]//N,OptionValue[startRegulator]//N,ReIm[OptionValue[regulator]]//N];
output=Table[Complex@@RawMemoryRead[res,i],{i,0,Length[integrals]-1}];
arrayDelF[res];
output
];

Options[splittingTreeIntegrateC]={epsAbs->10^-12,epsRel->10^-12,hStart->1*10^-5,hMax->1,startRegulator->0,baseRange->2*10^-5,regulator->1};
splittingTreeIntegrateC[integrands_List,integrals_List,zVal_,OptionsPattern[]]:=Module[{functions,funcMem,residues,resMem,integralLists,integralMemberMem,integralMem,res,output},
functions=If[MemberQ[{ManagedObject,ForeignCallback,OpaqueRawPointer},Head[#]],#,CreateForeignCallback[val|->ReIm[#[Complex@@val]],{{"CDouble","CDouble"}}->{"CDouble","CDouble"}]]&/@integrands[[All,1]];
residues=ReIm[integrands[[All,2]]]//N;
funcMem=RawMemoryExport[functions,"OpaqueRawPointer"];
resMem=RawMemoryExport[residues,{"CDouble","CDouble"}];
integralLists=Prepend[#-1,Length[#]]&/@integrals;
integralMemberMem=RawMemoryExport[#,"UnsignedInteger64"]&/@integralLists;
integralMem=RawMemoryExport[integralMemberMem,"RawPointer"::["UnsignedInteger64"]];
res=treeIntegrateSplitF[N[ReIm[zVal]],funcMem,resMem,Length[residues],integralMem,Length[integrals],OptionValue[epsAbs]//N,OptionValue[epsRel]//N,OptionValue[hStart]//N,OptionValue[hMax]//N,OptionValue[startRegulator]//N,OptionValue[baseRange]//N,ReIm[OptionValue[regulator]]//N];
output=Table[Complex@@RawMemoryRead[res,i],{i,0,Length[integrals]-1}];
arrayDelF[res];
output
];
];


(*calculate the shuffle broduct between two lists*)
shuffleProd[a_List, b_List] :=Module[{firstPositions, secondPositions, indexHelper, ord,combination},
If[a === {} === b, Return[{{}}]];
firstPositions = Permutations@Join[ConstantArray[1, Length[a]], ConstantArray[0, Length[b]]];
secondPositions = 1 - firstPositions;
indexHelper = Transpose[Accumulate[Transpose[#]]]*# &;
ord = indexHelper[firstPositions] + (indexHelper[secondPositions] + Length[a]) secondPositions;
combination = Join[a, b];
combination[[#]] & /@ ord
];

convWrap[fun_, res_] := (fun[#] -Quiet[res/#] &);
convWrap[fun_, 0] := fun;
rFun[asL_, zVal_] := Module[{},
If[Length[asL] == 0, Return[1]];
If[asL[[1, 2]] === 0,Return[integral[Sequence @@ asL[[All, 1]], zVal]]];
Sum[(-1)^(i + 1)*Sum[integral[convWrap @@ asL[[i]], Sequence @@ rem, zVal],
{rem,shuffleProd[(x |-> Quiet[#/x]) & /@ Reverse[asL[[;; i - 1, 2]]],asL[[i + 1 ;;, 1]]]}],
{i, 1,Min[Length[asL], LengthWhile[asL, #[[2]] =!= 0 &] + 1]}]
];

safePower[0, 0] := 1;
safePower[0., 0] := 1;
safePower[x_, n_] := x^n;

regulate[integrands__List, zVal_, regulator_] :=Module[{asL = {integrands}},
If[Length[asL] == 0, Return[1]];
If[zVal == 0, Return[0]];
Sum[1/k!*Product[asL[[j, 2]], {j, k}]*safePower[Log[zVal/regulator],k]*rFun[asL[[k + 1 ;;]], zVal], {k, 0, Length[asL]}]
];

Options[integCore] = Options[NDSolveValue];
integCore[integrands__, zVal_, startT_, stopT_, initialVals_, OptionsPattern[]] := Module[{intL = {integrands}, t, matr, functions, wrapper, args, f},
wrapper[f_, x_] /; NumericQ[x] :=f[x] /. {ComplexInfinity -> 0, Indeterminate -> 0};
matr =DiagonalMatrix[wrapper[#, zVal*t] & /@ intL, -1,TargetStructure -> "Sparse"][[2 ;;]];
functions = f[#][t] & /@ Range[Length[intL]];
args = Prepend[functions, 1];

NDSolveValue[{D[functions, t] == zVal*matr . args,initialVals == functions /. t -> startT}, functions /. t -> stopT, {t, startT, stopT},Sequence @@ (# -> OptionValue[#] & /@Options[NDSolveValue][[All, 1]])]
];

Options[execIntegration] =
  Join[Options[integCore], {startRegulator -> 0}];
execIntegration[integrands__, zVal_, OptionsPattern[]] := Module[{intL = {integrands}, t, matr, functions, wrapper, args, f},
If[Length[intL] == 0, Return[1]];
If[zVal == 0, Return[0]];
integCore[integrands, zVal, OptionValue[startRegulator], 1,ConstantArray[0, Length[intL]],
Sequence @@ (# -> OptionValue[#] & /@Options[integCore][[All, 1]])][[-1]]
];

Options[plainIntegrate] =
  Join[Options[execIntegration], {regulator -> 1, checkFunction -> (False &)}];
plainIntegrate[integrands__List, zVal_, OptionsPattern[]] := Module[{expr, ints, val, checkFct = OptionValue[checkFunction]},
If[Length[{integrands}] == 0, Return[1]];
If[zVal == 0, Return[0]];
expr=regulate[integrands, zVal, OptionValue[regulator]];
expr/.x:integral[integs__, zMax_] :>(val = checkFct[x];If[val=!=False, val, execIntegration[integs, zMax,Sequence @@ (# -> OptionValue[#] & /@Options[execIntegration][[All, 1]])]])
];

ClearAll[splittingPlainIntegrate];
Options[splittingPlainIntegrate] = Join[Options[execIntegration], {switchPoint -> 2*10^-5, regulator -> 1, startRegulator -> 0, checkFunction -> (False &)}];
splittingPlainIntegrate[integrands__List, zVal_, OptionsPattern[]] :=
 Module[{intL = {integrands}, exprs, integrals, evaluated = <||>,
    checkFct = OptionValue[checkFunction],val},
  If[Length[intL] == 0, Return[1]];
  If[zVal == 0, Return[0]];

  exprs =
   (regulate[Sequence @@ intL[[;; #]], zVal*OptionValue[switchPoint],
      OptionValue[regulator]] & /@ Range[Length[intL]]);
  exprs = exprs/.x_integral:>(val=checkFct[x];If[val=!=False,val,x]);
  integrals =
   Reverse[SortBy[DeleteDuplicates[Cases[exprs, integral[__], All]],
     Length]];

  Do[If[KeyExistsQ[evaluated, integ], Continue[]];
   AssociateTo[evaluated,
    Thread[(Join[integ[[;; #]], integ[[{-1}]]] & /@ Range[Length[integ] - 1]) ->
      integCore[Sequence @@ integ[[;; -2]], zVal,
       OptionValue[startRegulator], OptionValue[switchPoint], ConstantArray[0, Length[integ] - 1],
       Sequence @@ (# -> OptionValue[#] & /@ Options[integCore][[All, 1]])]]];
   , {integ, integrals}];

  integCore[Sequence @@ intL[[All, 1]], zVal,
    OptionValue[switchPoint], 1, exprs /. evaluated,
    Sequence @@ (# -> OptionValue[#] & /@ Options[integCore][[All, 1]])][[-1]]
  ]

convergent[x_, 0] := x;
rFunIdx[asL_, zVal_] := Module[{},
  If[Length[asL] == 0, Return[1]];
  If[asL[[1, 2]] == 0,
   Return[integral[Sequence @@ asL[[All, 1]], zVal]]];
  Sum[(-1)^(i + 1)*
    Sum[integral[convergent @@ asL[[i]], Sequence @@ rem, zVal], {rem,
       shuffleProd[poleTerm /@ Reverse[asL[[;; i - 1, 2]]],
       asL[[i + 1 ;;, 1]]]}], {i, 1,
    Min[Length[asL], LengthWhile[asL, #[[2]] != 0 &] + 1]}]
  ]

ClearAll[regulateIdx];
regulateIdx[integrands__List, zVal_, regulator_] :=
 Module[{asL = {integrands}},
  If[Length[asL] == 0, Return[1]];
  If[zVal == 0, Return[0]];
  Sum[1/k!*Product[asL[[j, 2]], {j, k}]*safePower[Log[zVal/regulator], k]*
    rFunIdx[asL[[k + 1 ;;]], zVal], {k, 0, Length[asL]}]
  ]


Options[treeIntegrateConvergent] =
  Join[Options[NDSolveValue], {startRegulator -> 0}];
treeIntegrateConvergent[integrandFunctions_List, integrals_List,
   zVal_, OptionsPattern[]] :=
  Module[{toCompute, computeGroups, integralSets, integralResults,
    fun, wrapper, t},
   wrapper[f_, x_] /; NumericQ[x] :=
    f[x] /. {ComplexInfinity -> 0, Indeterminate -> 0};
   (*gather all integrals including subintegrals necessary*)
   toCompute =
    Union @@
     Table[int[[;; #]] & /@ Range[1, Length[int]], {int, integrals}];

   (*group by innermost integrand*)
   computeGroups = GatherBy[toCompute, #[[1]] &];

   (*each integral is now only a pair of a function and the previous function (i.e. we have a tree)*)
   integralSets =
    Table[{integrandFunctions[[elem[[-1]]]], FirstPosition[grp, elem[[;; -2]], {Length[toCompute] + 1}, 1][[1]]},
      {grp, computeGroups}, {elem, grp}];
   (*root element is one*)
   fun[Length[toCompute] + 1] = 1 &;

   integralResults =
    Table[NDSolveValue[
      Flatten[{D[fun[#][t], t] ==
           wrapper[fctSet[[#, 1]], zVal*t]*zVal*
            fun[fctSet[[#, 2]]][t],
          fun[#][OptionValue[startRegulator]] == 0} & /@
        Range[Length[fctSet]]],
      fun[#][1] & /@ Range[Length[fctSet]], {t,
       OptionValue[startRegulator], 1},
      Sequence @@ (# -> OptionValue[#] & /@
         Options[NDSolveValue][[All, 1]])], {fctSet, integralSets}];

   integrals /. (Thread[(Join @@ computeGroups) -> (Join @@
         integralResults)])
   ];

Options[treeIntegrate] =
  Join[Options[treeIntegrateConvergent], {regulator -> 1,
    checkFunction -> (False &)}];
treeIntegrate[integrandFunctions_List, integrals_List, zVal_,
   OptionsPattern[]] :=
  Module[{regulatedExpressions, allIntegrals, checkFct = OptionValue[checkFunction], toFunctionRules, v, integrands, callIntegrals},
   toFunctionRules = {intFun[x_] :> (integrandFunctions[[x, 1]]),
     poleTerm[x_] :> (Quiet[x/#] &), convergent -> convWrap};
   (*regulate every integral keeping indices instead of function expressions*)
   regulatedExpressions =
    Table[regulateIdx[Sequence @@ ({intFun[#], integrandFunctions[[#, 2]]} & /@integ), zVal, OptionValue[regulator]],
    {integ, integrals}];
   (*Remove all easily calculated ones*)
   regulatedExpressions=regulatedExpressions/.x_integral:>(v = checkFct[x /. toFunctionRules];If[v=!=False, v, x]);
   
   (*List of all integrals*)
   allIntegrals = DeleteDuplicates[Cases[regulatedExpressions, _integral, All]];

   (*convert integrals to a call to the function for convergent integrals*)
   integrands = Union @@ List @@@ allIntegrals[[All, ;; -2]];
   callIntegrals =
    Table[FirstPosition[integrands, kernel, {}, 1][[1]], {int,
      allIntegrals[[All, ;; -2]]}, {kernel, List @@ int}];

   regulatedExpressions /.
    Thread[allIntegrals ->
      treeIntegrateConvergent[integrands /. toFunctionRules,
       callIntegrals, zVal,
       Sequence @@ (# -> OptionValue[#] & /@
          Options[treeIntegrateConvergent][[All, 1]])]]
   ];

Options[splittingTreeIntegrate] =
  Join[Options[treeIntegrate], {switchPoint -> 2*10^-5}];
splittingTreeIntegrate[integrands_List, integrals_List, zVal_,
  OptionsPattern[]] :=
 Module[{integVals, toCompute, computeGroups, integralSets, fun,
   integralResults, wrapper, t},
  If[zVal == 0, Return[ConstantArray[0, Length[integrals]]]];
  wrapper[f_, x_] /; NumericQ[x] :=
   f[x] /. {ComplexInfinity -> 0, Indeterminate -> 0};

  toCompute =
   Union @@ Table[int[[;; i]], {int, integrals}, {i, Length[int]}];
  integVals =
   treeIntegrate[integrands, toCompute, zVal*OptionValue[switchPoint],
     Sequence @@ (# -> OptionValue[#] & /@
       Options[treeIntegrate][[All, 1]])];

  (*group by innermost integrand*)
  computeGroups = GatherBy[toCompute, #[[1]] &];
  (*each integral is now only a pair of a function and the previous function (i.e. we have a tree)*)
  integralSets =
   Table[{integrands[[elem[[-1]], 1]],
     FirstPosition[grp, elem[[;; -2]], {Length[toCompute] + 1}, 1][[
      1]], integVals[[
      FirstPosition[toCompute, elem, {}, 1][[1]]]]}, {grp,
     computeGroups}, {elem, grp}];
  (*root element is one*)
  fun[Length[toCompute] + 1] = 1 &;

  integralResults =
   Table[NDSolveValue[
     Flatten[{D[fun[#][t], t] ==
          wrapper[fctSet[[#, 1]], zVal*t]*zVal*fun[fctSet[[#, 2]]][t],
          fun[#][OptionValue[switchPoint]] == fctSet[[#, 3]]} & /@
       Range[Length[fctSet]]],
     fun[#][1] & /@ Range[Length[fctSet]], {t,
      OptionValue[switchPoint], 1},
     Sequence @@ (# -> OptionValue[#] & /@
        Options[NDSolveValue][[All, 1]])], {fctSet, integralSets}];

  integrals /. (Thread[(Join @@ computeGroups) -> (Join @@
        integralResults)])
  ]


IIntegral[z_]:=1;
kernelResidue[regulated[name_]]:=0;
kernelResidue[pole[name_]]:=kernelResidue[name];
kernelExpression[regulated[name_]]:=kernelExpression[name][#] -Quiet[kernelResidue[name]/#] &;
kernelExpression[pole[name_]]:=Quiet[kernelResidue[name]/#]&;


convWrapExt[name_]:=Block[{},
If[kernelResidue[name]==0,Return[name]];
If[MemberQ[kernels,regulated[name]],Return[regulated[name]];];
Return[regulated[name]];
];
resKernel[name_]:=Block[{},
If[MemberQ[kernels,pole[name]],Return[pole[name]];];
Return[pole[name]];
];


rFunExt[asL_List, zVal_] := Module[{},
If[Length[asL] == 0, Return[1]];
If[kernelResidue[asL[[1]]] === 0,Return[IIntegral[Sequence @@ asL[[All]], zVal]]];
Sum[(-1)^(i + 1)*Sum[IIntegral[convWrapExt[asL[[i]]], Sequence @@ rem, zVal],
{rem,shuffleProd[resKernel /@ Reverse[asL[[;; i - 1]]],asL[[i + 1 ;;]]]}],
{i, 1,Min[Length[asL], LengthWhile[asL, kernelResidue[asL] =!= 0 &] + 1]}]
];

regulateIIntegral[funcs__, zVal_, regulator_] :=Module[{asL={funcs}},
If[Length[asL] == 0, Return[1]];
If[zVal == 0, Return[0]];
Sum[1/k!*Product[kernelResidue[asL[[j]]], {j, k}]*safePower[Log[zVal/regulator],k]*rFunExt[asL[[k + 1 ;;]], zVal], {k, 0, Length[asL]}]
];

regulateExpression[expr_,regulator_:1]:=Block[{},
expr/.IIntegral[funcs___,zMax_]:>regulateIIntegral[funcs,zMax,regulator]
];


kernels={};

defineKernel::AlreadyDefined="There already exists a kernel with name `1`, which will be overwritten.";
defineKernel::MissingEntry="The Association `1` is missing the key `2` to define a kernel!";

defineKernel[assoc_Association]:=Block[{},
If[!KeyMemberQ[assoc,name],Message[defineKernel::MissingEntry,assoc,name]; Return[$Failed]];
If[!KeyMemberQ[assoc,expression],Message[defineKernel::MissingEntry,assoc,expression]; Return[$Failed]];
If[MemberQ[kernels,assoc[name]],Message[defineKernel::AlreadyDefined,assoc[name]]];
AppendTo[kernels,assoc[name]];
kernelExpression[assoc[name]]=assoc[expression];
If[KeyMemberQ[assoc,residue],
kernelResidue[assoc[name]]=assoc[residue],
kernelResidue[assoc[name]]=0
];
];

defineKernel[namee_,callable_,residuee_:0]:=defineKernel[<|name->namee,expression->callable,residue->residuee|>];


listKernels[]:=kernels;
exportKernel[namee_]:=<|name->namee,expression->kernelExpression[namee],residue->kernelResidue[namee]|>;
exportKernels[]:=(exportKernel/@kernels);
deleteKernel[name_]:=Block[{},kernelResidue[name]=.;kernelExpression[name]=.;kernels=DeleteCases[kernels,name];]


performGroupIntegration[expr_,integrator_,externalRules_,integratorOptions_]:=Module[{integrals,groups,rules,kers,zMax,vals},
integrals=DeleteDuplicates[Cases[expr,_IIntegral,All]];
groups=GroupBy[integrals,Last];
rules={};
Do[
zMax=grp[[1,-1]];
kers=DeleteDuplicates[Flatten[Cases[grp,IIntegral[funcs__,z_]:>{funcs},All]]];
vals=integrator[({kernelExpression[#],kernelResidue[#]}/.externalRules)&/@kers,(int|->FirstPosition[kers,#,{Missing},{1}][[1]]&/@((List@@int)[[;;-2]]))/@grp, zMax, Sequence @@integratorOptions];
rules=Join[rules,Thread[grp->vals]];
,{grp,groups}];
Return[expr/.rules];
];


Options[computeIIntegrals]=Union[Options[plainIntegrate],Options[splittingPlainIntegrate],Options[treeIntegrate],Options[splittingTreeIntegrate],{integrationAlgorithm->treeIntegrate, Quiet->True}];
computeIIntegrals::unknownMethod="integrationAlgorithm `1` is not known. Possible values are: treeIntegrate, plainIntegrate, splittingPlainIntegrate, splittingTreeIntegrate.";
computeIIntegrals[expr_,rules_List:{},OptionsPattern[]]:=Module[{integrals,vals,temp,quietWrap=ReleaseHold},
If[OptionValue[Quiet],quietWrap=Quiet[ReleaseHold[#],{Power::infy,Infinity::indet}]&];
temp=expr/.rules;

If[OptionValue[integrationAlgorithm]===plainIntegrate,
integrals=DeleteDuplicates[Cases[temp,_IIntegral,All]];
quietWrap[Hold[vals=integrals/.IIntegral[funs___,limit_]:>plainIntegrate[Sequence@@(({kernelExpression[#],kernelResidue[#]}/.rules)&/@{funs}),limit/.rules,Sequence @@ (# -> OptionValue[#] & /@Options[plainIntegrate][[All, 1]])]]];
Return[temp/.Thread[integrals->vals]];];

If[OptionValue[integrationAlgorithm]===splittingPlainIntegrate,
integrals=DeleteDuplicates[Cases[temp,_IIntegral,All]];
quietWrap[Hold[vals=integrals/.IIntegral[funs___,limit_]:>splittingPlainIntegrate[Sequence@@(({kernelExpression[#],kernelResidue[#]}/.rules)&/@{funs}),limit/.rules,Sequence @@ (# -> OptionValue[#] & /@Options[splittingPlainIntegrate][[All, 1]])]]];
Return[temp/.Thread[integrals->vals]];];

If[OptionValue[integrationAlgorithm]===treeIntegrate,
Return[quietWrap[Hold[performGroupIntegration[temp/.IIntegral[funs___,limit_]:>IIntegral[funs,limit/.rules],treeIntegrate,rules,(# -> OptionValue[#] & /@Options[treeIntegrate][[All, 1]])]]]];
];

If[OptionValue[integrationAlgorithm]===splittingTreeIntegrate,
Return[quietWrap[Hold[performGroupIntegration[temp/.IIntegral[funs___,limit_]:>IIntegral[funs,limit/.rules],splittingTreeIntegrate,rules,(# -> OptionValue[#] & /@Options[splittingTreeIntegrate][[All, 1]])]]]];
];

Message[computeIIntegrals::unknownMethod,OptionValue[integrationAlgorithm]];
Return[$Failed];
];


If[Global`$NativeIteratedIntegral===True,
Options[computeIIntegralsC]=Union[Options[splittingTreeIntegrateC],{integrationAlgorithm->treeIntegrate}];
computeIIntegralsC[expr_,rules_List:{},OptionsPattern[]]:=Module[{integrals,vals,temp},
temp=expr/.rules;
If[OptionValue[integrationAlgorithm]===plainIntegrate,
integrals=DeleteDuplicates[Cases[temp,_IIntegral,All]];
vals=integrals/.IIntegral[funs___,limit_]:>plainIntegrateC[Sequence@@(({kernelExpression[#],kernelResidue[#]}/.rules)&/@{funs}),limit,Sequence @@ (# -> OptionValue[#] & /@Options[plainIntegrateC][[All, 1]])];
Return[temp/.Thread[integrals->vals]];];

If[OptionValue[integrationAlgorithm]===splittingPlainIntegrate,
integrals=DeleteDuplicates[Cases[temp,_IIntegral,All]];
vals=integrals/.IIntegral[funs___,limit_]:>splittingPlainIntegrateC[Sequence@@(({kernelExpression[#],kernelResidue[#]}/.rules)&/@{funs}),limit/.rules,Sequence @@ (# -> OptionValue[#] & /@Options[splittingPlainIntegrateC][[All, 1]])];
Return[temp/.Thread[integrals->vals]];];

If[OptionValue[integrationAlgorithm]===treeIntegrate,
Return[performGroupIntegration[temp/.IIntegral[funs___,limit_]:>IIntegral[funs,limit],treeIntegrateC,rules,(# -> OptionValue[#] & /@Options[treeIntegrateC][[All, 1]])]];
];

If[OptionValue[integrationAlgorithm]===splittingTreeIntegrate,
Return[performGroupIntegration[temp/.IIntegral[funs___,limit_]:>IIntegral[funs,limit],splittingTreeIntegrateC,rules,(# -> OptionValue[#] & /@Options[splittingTreeIntegrateC][[All, 1]])]];
];

Message[computeIIntegrals::unknownMethod,OptionValue[integrationAlgorithm]];
Return[$Failed];
];];


Format[pole[name_],TraditionalForm]:=Superscript[name,\[Infinity]];
Format[regulated[name_],TraditionalForm]:=Superscript[name,"R"];
Format[IIntegral,TraditionalForm]:="\[ScriptCapitalI]";


End[];
EndPackage[];
