(* ::Package:: *)

(* ::Title:: *)
(*GAMI package*)


(* ::Subtitle:: *)
(*tools for analysis of origami*)


BeginPackage["GAMI`"]


(* ::Chapter:: *)
(*usage*)


GAMI::usage=
"GAMI contains a number of tools for folding origami sheets."


(* ::Section::Closed:: *)
(*patterns*)


Network::usage=
"Network[vertices_,topology_] specifies a sheet as a network of sites located at vertices connected by bonds according to topology."


Origami::usage=
"Origami[sectors_,folds_,lengths_,topology_] specifies a sheet as a collection of vertices connected by creases according to topology."


Node::usage=
"Specifies a node on a sheet."


Crease::usage=
"Specifies a crease on a sheet."


Face::usage=
"Specifies a face on a sheet."


Unit::usage=
"Specifies a cell of the sheet."


Tessellation::usage=
"Specifies a tessellation of the sheet."


Triangulation::usage=
"Triangulation[connections_,periodicity_] specifies a triangulated topology."


Parallelogram::usage=
"Parallelogram[connections_,periodicity_] specifies a parallelogram topology."


Surface::usage=
"Specifies a continuous surface."


Plane::usage=
"Specifies a planar surface."


Uniaxial::usage=
"Specifies a uniaxial surface."


Spatial::usage=
"Spatial[translations_] specifies a spatially periodic topology."


Screw::usage=
"Screw[translations_,angles_,axis_] specifies a screw periodic topology."


FoldingMotion::usage=
"Describes an infinitesimal change to the fold angles on a sheet."


(* ::Section::Closed:: *)
(*queries*)


SheetQ::usage=
"SheetQ[obj_] tests whether obj is the characterization of a sheet."


NetworkQ::usage=
"NetworkQ[obj_] tests whether obj characterization a sheet as a network."


OrigamiQ::usage=
"OrigamiQ[obj_] tests whether obj characterization a sheet as origami."


NodeQ::usage=
"Tests whether obj is a node."


CreaseQ::usage=
"Tests whether obj is a crease."


FaceQ::usage=
"Tests whether obj is a face."


TopologyQ::usage=
"Tests whether obj defines a triangulated or parallelogram topology."


TriangulationQ::usage=
"TriangulationQ[obj_] tests whether obj is a triangulation."


ParallelogramQ::usage=
"ParallelogramQ[obj_] tests whether obj is a planar parallelogram."


PeriodicQ::usage=
"Tests whether obj defines spatial or screw periodicity."


SpatialQ::usage=
"SpatialQ[obj_] tests whether obj is a spatially periodic."


ScrewQ::usage=
"ScrewQ[obj_] tests whether obj is a screw periodic."


(* ::Section::Closed:: *)
(*computations*)


AxisTransform::usage=
"Constructs the similarity transform rotating the sheet's screw axis to the z-axis."


RotSum::usage=
"Sums rotation matrices about the z-axis."


TransverseVec::usage=
"Returns the project of a vector into the plane defined by an axis."


(* ::Chapter:: *)
(*errors*)


(* ::Chapter:: *)
(*routines*)


Begin["`Private`"]


(* ::Section:: *)
(*patterns*)


(* ::Subsection:: *)
(*Network*)


(* ::Subsubsection:: *)
(*Fundamentals*)


Network[vertices_,topology_]["vertices"[]]:=vertices
Network[vertices_,topology_]["topology"[]]:=topology
Network[vertices_,topology_]["connections"[]]:=topology["connections"[]]
Network[vertices_,topology_]["labels"[]]:=topology["labels"[]]
Network[vertices_,topology_]["faces"[]]:=topology["faces"[]]
Network[vertices_,topology_]["periodicity"[]]:=topology["periodicity"[]]
Network[vertices_,topology_]["translations"[]]:=topology["periodicity"[]]["translations"[]]/;PeriodicQ[topology["periodicity"[]]]
Network[vertices_,topology_]["angles"[]]:=topology["periodicity"[]]["angles"[]]/;ScrewQ[topology["periodicity"[]]]
Network[vertices_,topology_]["axis"[]]:=topology["periodicity"[]]["axis"[]]/;ScrewQ[topology["periodicity"[]]]


(* ::Subsubsection::Closed:: *)
(*Rigidity Matrix*)


Network[vertices_,topology_]["RigidityMatrix"]:=SparseArray@Flatten[ Last[ Reap[ 
Do[ 
With[
{
dir=Crease[Network[vertices,topology],i]["direction"[]]
},
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,1,1]]+#-3}->-dir[[#]]&/@Range[3]]; 
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,2,1]]+#-3}->dir[[#]]&/@Range[3]]], 
{i,Length[Network[vertices,topology]["connections"[]]]}]
] ] ]/;SpatialQ[topology["periodicity"[]]]


Network[vertices_,topology_]["RigidityMatrix"[z1_,z2_]]:=SparseArray@Flatten[ Last[ Reap[ 
Do[ 
With[
{
dir=Crease[Network[vertices,topology],i]["direction"[]], 
coefs=z1^First[Last[#]] z2^Last[Last[#]]&/@Network[vertices,topology]["connections"[]][[i]]
},
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,1,1]]+#-3}->-First[coefs]dir[[#]]&/@Range[3]]; 
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,2,1]]+#-3}->Last[coefs]dir[[#]]&/@Range[3]]], 
{i,Length[Network[vertices,topology]["connections"[]]]}]
] ] ]/;SpatialQ[topology["periodicity"[]]]


Network[vertices_,topology_]["RigidityMatrix"[]]:=SparseArray@Flatten[ Last[ Reap[ 
Do[ 
With[
{
dir=Crease[Network[vertices,topology],i]["direction"[]], 
rot=RotationMatrix[-topology["angles"[]].Flatten[Differences@(Last/@topology["connections"[]][[i]])],topology["axis"[]]]
}, 
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,1,1]]+#-3}->-dir[[#]]&/@Range[3]]; 
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,2,1]]+#-3}->(rot.dir)[[#]]&/@Range[3]]], 
{i,Length[Network[vertices,topology]["connections"[]]]}]
] ] ]/;ScrewQ[topology["periodicity"[]]]


Network[vertices_,topology_]["RigidityMatrix"[z1_,z2_]]:=SparseArray@Flatten[ Last[ Reap[ 
Do[ 
With[
{
dir=Crease[Network[vertices,topology],i]["direction"[]], 
rot=RotationMatrix[-topology["angles"[]].Flatten[Differences@(Last/@topology["connections"[]][[i]])],topology["axis"[]]],
coefs=z1^First[Last[#]] z2^Last[Last[#]]&/@Network[vertices,topology]["connections"[]][[i]]
},
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,1,1]]+#-3}->-First[coefs]dir[[#]]&/@Range[3]]; 
Sow[{i,3 Network[vertices,topology]["connections"[]][[i,2,1]]+#-3}->Last[coefs](rot.dir)[[#]]&/@Range[3]]], 
{i,Length[Network[vertices,topology]["connections"[]]]}]
] ] ]/;ScrewQ[topology["periodicity"[]]]


(* ::Subsubsection::Closed:: *)
(*Equilibrium Matrix*)


Network[vertices_,topology_]["EquilibriumMatrix"]:=Transpose[ Network[vertices,topology]["RigidityMatrix"] ]
Network[vertices_,topology_]["EquilibriumMatrix"[z1_,z2_]]:=Transpose[ Network[vertices,topology]["RigidityMatrix"[1/z1, 1/z2]] ]


(* ::Subsubsection::Closed:: *)
(*Angular Velocity*)


Network[vertices_,topology_]["CellAngularVelocities"[paths_,motion_]]:=FoldingMotion[Network[vertices,topology],motion]["CellAngularVelocity"[#]]&/@First[paths]
Network[vertices_,topology_]["CellAngularVelocities"[paths_,motion_,cell_]]:=FoldingMotion[Network[vertices,topology],motion]["CellAngularVelocity"[#,cell]]&/@First[paths]
Network[vertices_,topology_]["AngularVelocity"[paths_,index_,cell_,motion_]]:=Face[Network[vertices,topology],index,cell]["AngularVelocity"[paths,motion]]


Network[vertices_,topology_]["CellAngularVelocities"[paths_,motion_,z1_,z2_]]:=FoldingMotion[Network[vertices,topology],motion]["CellAngularVelocity"[#]]&/@First[paths]
Network[vertices_,topology_]["CellAngularVelocities"[paths_,motion_,cell_,z1_,z2_]]:=FoldingMotion[Network[vertices,topology],motion]["CellAngularVelocity"[#,cell]]&/@First[paths]
Network[vertices_,topology_]["AngularVelocity"[paths_,index_,cell_,motion_,z1_,z2_]]:=Face[Network[vertices,topology],index,cell]["AngularVelocity"[paths,motion,z1,z2]]


(* ::Subsubsection::Closed:: *)
(*Displacement*)


Network[vertices_,topology_]["CreaseRotation"[index_,cell_,paths_,face_,motion_]]:=Crease[Network[vertices,topology],index,cell]["Rotation"[paths,face,motion]]
Network[vertices_,topology_]["CellDisplacements"[facePaths_,paths_,motion_]]:=FoldingMotion[Network[vertices,topology],motion]["CellDisplacement"[facePaths,#]]&/@First[paths]
Network[vertices_,topology_]["Displacement"[paths_,index_,cell_,motion_]]:=Node[Network[vertices,topology],index,cell]["Displacement"[paths,motion]]


(* ::Subsubsection::Closed:: *)
(*Origami*)


Network[vertices_,topology_]["SectorAngles"[]]:=Module[
{
edgeDirs=Table[If[MatchQ[Last[Network[vertices,topology]["connections"[]][[#]]],{i,{_,_}}]==False,+1,-1]Crease[Network[vertices,topology],#]["direction"[]]&/@topology["labels"[]][[i]],{i,Length[topology["labels"[]]]}]
},
Table[VectorAngle[#[[1]],#[[2]]]&/@( Transpose[{#, Append[Delete[#, 1],First[#]]}&@edgeDirs[[i]]] ),{i,Length[topology["labels"[]]]}]
]


Network[vertices_,topology_]["FoldAngles"[]]:=Module[
{
edgeDirs=Table[If[MatchQ[Last[Network[vertices,topology]["connections"[]][[#]]],{i,{_,_}}]==False,+1,-1]Crease[Network[vertices,topology],#]["direction"[]]&/@topology["labels"[]][[i]],{i,Length[topology["labels"[]]]}],
faceDirs,
projDirs
},
faceDirs=Table[Normalize[Cross[#[[1]],#[[2]]]]&/@( Transpose[{Prepend[Delete[#, -1],Last[#]],#}&@edgeDirs[[i]]] ),{i,Length[topology["labels"[]]]}];
projDirs=Table[Normalize[Cross[#[[1]],#[[2]]]]&/@( Transpose[{faceDirs[[i]],edgeDirs[[i]]} ]),{i,Length[topology["labels"[]]]}];
Table[ArcTan[#[[1]].#[[2]],#[[2]].#[[3]]]&/@( Transpose[{#,Append[Delete[#, 1],First[#]],projDirs[[i]]}&@faceDirs[[i]]] ),{i,Length[topology["labels"[]]]}]
]


Network[vertices_,topology_]["Origami"[]]:=Origami[#[[1]]["SectorAngles"[]],#[[1]]["FoldAngles"[]],Table[Crease[#[[1]],#[[2,i,j]]]["length"[]],{i,Length[#[[2]]]},{j,Length[#[[2,i]]]}],topology]&@{Network[vertices,topology],topology["labels"[]]}




(* ::Subsection:: *)
(*Origami*)


(* ::Subsubsection::Closed:: *)
(*Fundamentals*)


Origami[sectors_,folds_,lengths_,topology_]["sectors"[]]:=sectors
Origami[sectors_,folds_,lengths_,topology_]["folds"[]]:=folds
Origami[sectors_,folds_,lengths_,topology_]["lengths"[]]:=lengths
Origami[sectors_,folds_,lengths_,topology_]["topology"[]]:=topology


(* ::Subsubsection::Closed:: *)
(*Constraints*)


Origami[sectors_,folds_,lengths_,topology_]["VertexConstraint"[]]:=Module[
{
edgeRotations=Table[
RotationMatrix[folds[[i,j]],{1,0,0}],
{i,Length[topology["labels"[]]]},{j,Length[topology["labels"[]][[i]]]}],

faceRotations=Table[
RotationMatrix[sectors[[i,j]],{0,0,1}],
{i,Length[topology["labels"[]]]},{j,Length[topology["labels"[]][[i]]]}]
},
Table[
Chop[Dot@@Table[
edgeRotations[[i,j]].faceRotations[[i,j]],
{j,Length[topology["labels"[]][[i]]]}]],
{i,Length[topology["labels"[]]]}]
]


Origami[sectors_,folds_,lengths_,topology_]["VertexEnergy"[]]:=Module[{
constraintViolation=Origami[sectors,folds,lengths,topology]["VertexConstraint"[]]
},
Chop[Sum[(KroneckerDelta[j,k]-constraintViolation[[i,j,k]])^2,{i,Length[topology["labels"[]]]},{j,3},{k,3}]^(1/2)]
]


Origami[sectors_,folds_,lengths_,topology_]["LinearConstraint"[]]:=Module[
{
tempVertexConstraintExpansion=Chop[SeriesCoefficient[Origami[sectors,folds+\[Epsilon] Table[ToExpression["\[Phi]"<>ToString[topology["labels"[]][[i,j]]]],{i,Length[topology["labels"[]]]},{j,Length[topology["labels"[]][[i]]]}],lengths,topology]["VertexConstraint"[]],{\[Epsilon],0,1}]]
},
Flatten[Table[
Transpose[Table[{
Coefficient[tempVertexConstraintExpansion[[j,3,2]],ToExpression[StringJoin["\[Phi]"<>ToString[i]]]],
Coefficient[tempVertexConstraintExpansion[[j,1,3]],ToExpression[StringJoin["\[Phi]"<>ToString[i]]]],
Coefficient[tempVertexConstraintExpansion[[j,2,1]],ToExpression[StringJoin["\[Phi]"<>ToString[i]]]]},
{i,Max[topology["labels"[]]]}]],
{j,Length[topology]}],1]
]


(* ::Subsubsection:: *)
(*Evolution*)


Origami[sectors_,folds_,lengths_,topology_]["RigidFolding"[foldDir_,stepSize_]]:={Origami[sectors,folds+stepSize Table[ToExpression["\[Phi]"<>ToString[topology["labels"[]][[i,j]]]],{i,Length[topology["labels"[]]]},{j,Length[topology["labels"[]][[i]]]}],lengths,topology],Table[ToExpression["\[Phi]"<>ToString[i]],{i,Max[topology["labels"[]]]}]}/.Last[Quiet[
FindMinimum[
Origami[sectors,folds+stepSize Table[ToExpression["\[Phi]"<>ToString[topology["labels"[]][[i,j]]]],{i,Length[topology["labels"[]]]},{j,Length[topology["labels"[]][[i]]]}],lengths,topology]["VertexEnergy"[]],
Transpose[{Table[ToExpression["\[Phi]"<>ToString[i]],{i,Max[topology["labels"[]]]}],foldDir}],
MaxIterations->10^5],
{FindMinimum::lstol}]]


(* ::Subsection::Closed:: *)
(*Simplexes*)


(* ::Subsubsection::Closed:: *)
(*node*)


Node[sheet_,index_]["sheet"[]]:=sheet
Node[sheet_,index_,cell_]["sheet"[]]:=sheet
Node[sheet_,index_]["index"[]]:=index
Node[sheet_,index_,cell_]["index"[]]:=index
Node[sheet_,index_,cell_]["cell"[]]:=cell
Node[sheet_,index_]["position"[]]:=#["vertices"[]][[index]]&@sheet


Node[sheet_,index_,cell_]["position"[]]:=#["vertices"[]][[index]]+cell.#["translations"[]]&@sheet/;SpatialQ[sheet["periodicity"[]]]
Node[sheet_,index_,cell_]["position"[]]:=With[{
trans=AxisTransform[sheet],
rots=RotationMatrix[First[#]Last[#],sheet["axis"[]]]&/@Transpose[{sheet["angles"[]],cell}]
},
(Dot@@rots).#["vertices"[]][[index]]+Inverse[trans].RotSum[First[#["angles"[]]],First[cell]].trans.First[#["translations"[]]]+First[rots].Inverse[trans].RotSum[Last[#["angles"[]]],Last[cell]].trans.Last[#["translations"[]]]&@sheet
]/;ScrewQ[sheet["periodicity"[]]]


Node[sheet_,index_,cell_]["complex"[]]:=Point[Node[sheet,index,cell]["position"[]]]
Node[sheet_,index_,cell_]["projection"[]]:=#[[1]]-#[[1]].#[[2]] #[[2]]&@{Node[sheet,index,cell]["position"[]],Plane[sheet]["orientation"[]]}/;SpatialQ[sheet["periodicity"[]]]
Node[sheet_,index_,cell_]["projection"[]]:={Uniaxial[sheet]["radius"[]],cell.sheet["angles"[]]+VectorAngle[TransverseVec[-Uniaxial[sheet]["center"[]],sheet["axis"[]]],TransverseVec[sheet["vertices"[]][[index]]-Uniaxial[sheet]["center"[]],sheet["axis"[]]]],sheet["axis"[]].Node[sheet,index,cell]["position"[]]}/;ScrewQ[sheet["periodicity"[]]]


(* ::Item:: *)
(*displacement*)


Node[sheet_,index_,cell_]["Displacement"[facePaths_,creasePaths_,motion_]]:=With[{
rotsum=Table[Inverse[#].RotSum[sheet["angles"[]][[i]],cell[[i]]].#,{i,2}]&@AxisTransform[sheet],
rots=RotationMatrix[First[#]Last[#],sheet["axis"[]]]&/@Transpose[{sheet["angles"[]],cell}],
celldisp=sheet["CellDisplacements"[facePaths,creasePaths,motion]],
cellvel=sheet["CellAngularVelocities"[facePaths,motion]]
},
If[index==1,0.,Total@(Crease[sheet,Last[#],cell]["Rotation"[facePaths,First[#],motion]]&/@creasePaths[[index]])]
+First[rotsum].First[celldisp]+First[rots].Last[rotsum].Last[celldisp]
+First[rotsum].Cross[First[cellvel],First[sheet["translations"[]]]]+First[rots].Last[rotsum].Cross[Last[cellvel],Last[sheet["translations"[]]]]+Cross[First[rotsum].First[cellvel],First[rots].Last[rotsum].Last[sheet["translations"[]]]]
]


(* ::Subsubsection::Closed:: *)
(*crease*)


Crease[sheet_,index_]["sheet"[]]:=sheet
Crease[sheet_,index_,cell_]["sheet"[]]:=sheet
Crease[sheet_,index_]["index"[]]:=index
Crease[sheet_,index_,cell_]["index"[]]:=index
Crease[sheet_,index_,cell_]["cell"[]]:=cell


Crease[sheet_,index_]["vector"[]]:=Flatten[Differences[Node[sheet,First[#],Last[#]]["position"[]]&/@Extract[sheet["connections"[]],index]]]
Crease[sheet_,index_,cell_]["vector"[]]:=Flatten[Differences[Node[sheet,First[#],cell+Last[#]]["position"[]]&/@Extract[sheet["connections"[]],index]]]
Crease[sheet_,index_]["direction"[]]:=Normalize[Crease[sheet,index]["vector"[]]]
Crease[sheet_,index_,cell_]["direction"[]]:=Normalize[Crease[sheet,index,cell]["vector"[]]]
Crease[sheet_,index_]["length"[]]:=Norm[Crease[sheet,index]["vector"[]]]
Crease[sheet_,index_,cell_]["length"[]]:=Norm[Crease[sheet,index,cell]["vector"[]]]


Crease[sheet_,index_,cell_]["complex"[]]:=Line[Node[sheet,First[#],cell+Last[#]]["position"[]]&/@sheet["connections"[]][[index]]]


Crease[sheet_,index_,cell_]["projection"[]]:=Node[sheet,sheet["connections"[]][[index,1,1]],cell]["projection"[]]+{0.,VectorAngle[TransverseVec[-Uniaxial[sheet]["center"[]],sheet["axis"[]]],TransverseVec[sheet["vertices"[]][[sheet["connections"[]][[index,1,1]]]]+1/2 Crease[sheet,index,cell]["vector"[]]-Uniaxial[sheet]["center"[]],sheet["axis"[]]]],1/2 Crease[sheet,index,cell]["vector"[]].sheet["axis"[]]}/;ScrewQ[sheet["periodicity"[]]]


(* ::Item:: *)
(*rotation*)


Crease[sheet_,index_,cell_]["Rotation"[paths_,face_,motion_]]:=Cross[ sheet["AngularVelocity"[paths,face,cell,motion]], Crease[sheet,index,cell]["vector"[]] ]
Crease[sheet_,index_,cell_]["Rotation"[paths_,face_,motion_,cell_]]:=Cross[ sheet["AngularVelocity"[paths,face,cell,motion]], Crease[sheet,index,cell]["vector"[]] ]


Crease[sheet_,index_,cell_]["Rotation"[paths_,face_,motion_,z1_,z2_]]:=Cross[ sheet["AngularVelocity"[paths,face,cell,motion,z1,z2]], Crease[sheet,index,cell]["vector"[]] ]



(* ::Subsubsection::Closed:: *)
(*face*)


Face[sheet_,index_]["sheet"[]]:=sheet
Face[sheet_,index_,cell_]["sheet"[]]:=sheet
Face[sheet_,index_]["index"[]]:=index
Face[sheet_,index_,cell_]["index"[]]:=index
Face[sheet_,index_,cell_]["cell"[]]:=cell
Face[sheet_,index_]["complex"[]]:=Polygon[Node[sheet,First[#],Last[#]]["position"[]]&/@sheet["faces"[]][[index]]]
Face[sheet_,index_,cell_]["complex"[]]:=Polygon[Node[sheet,First[#],cell+Last[#]]["position"[]]&/@sheet["faces"[]][[index]]]


(* ::Item:: *)
(*angular velocty*)


Face[sheet_,index_,cell_]["AngularVelocity"[paths_,motion_]]:=With[
{
trans=AxisTransform[sheet],
rots=RotationMatrix[First[#]Last[#],sheet["axis"[]]]&/@Transpose[{sheet["angles"[]],cell}],
cellvel=sheet["CellAngularVelocities"[paths,motion]]
},
If[index!=1,
(Dot@@rots).(Total@( Crease[sheet,First[#],Last[#]]["direction"[]]motion[[First[#]]]&/@paths[[index]]) )
+Inverse[trans].RotSum[First[#["angles"[]]],First[cell]].trans.First[cellvel]
+First[rots].Inverse[trans].RotSum[Last[#["angles"[]]],Last[cell]].trans.Last[cellvel],
Inverse[trans].RotSum[First[#["angles"[]]],First[cell]].trans.First[cellvel]
+First[rots].Inverse[trans].RotSum[Last[#["angles"[]]],Last[cell]].trans.Last[cellvel]
]&@sheet
]


Face[sheet_,index_,cell_]["AngularVelocity"[paths_,motion_,z1_,z2_]]:=With[
{
trans=AxisTransform[sheet],
rots=RotationMatrix[First[#]Last[#],sheet["axis"[]]]&/@Transpose[{sheet["angles"[]],cell}],
cellvel=sheet["CellAngularVelocities"[paths,motion,z1,z2]]
},
If[index!=1,
Re[
z1^First[cell] z2^Last[cell] (Dot@@rots).(Total@( Crease[sheet,First[#],Last[#]]["direction"[]]motion[[First[#]]]&/@paths[[index]]) )
+Inverse[trans].RotSum[First[#["angles"[]]],First[cell],z1].trans.First[cellvel]
+z1^First[cell] First[rots].Inverse[trans].RotSum[Last[#["angles"[]]],Last[cell],z2].trans.Last[cellvel]
],
Re[
Inverse[trans].RotSum[First[#["angles"[]]],First[cell],z2].trans.First[cellvel]
+z1^First[cell] First[rots].Inverse[trans].RotSum[Last[#["angles"[]]],Last[cell],z2].trans.Last[cellvel]
]
]&@sheet
]


(* ::Subsubsection::Closed:: *)
(*cell*)


Unit[sheet_]["sheet"[]]:=sheet
Unit[sheet_,cell_]["sheet"[]]:=sheet
Unit[sheet_,cell_]["cell"[]]:=cell
Unit[sheet_]["complex"[]]:=Table[Face[sheet,i]["complex"[]],{i,Length@sheet["faces"[]]}]
Unit[sheet_,cell_]["complex"[]]:=Table[Face[sheet,i,cell]["complex"[]],{i,Length@sheet["faces"[]]}]


(* ::Subsubsection::Closed:: *)
(*tessellation*)


Tessellation[sheet_,cells_]["sheet"[]]:=sheet
Tessellation[sheet_,cells_]["cells"[]]:=cells
Tessellation[sheet_,cells_]["complex"[]]:=Flatten[Table[Unit[sheet,{n1,n2}]["complex"[]],{n1,0,First[cells]-1},{n2,0,Last[cells]-1}]]


(* ::Subsection:: *)
(*Topology*)


(* ::Subsubsection:: *)
(*Triangulation*)


Triangulation[connections_,labels_,faces_,periodicity_]["connections"[]]:=connections
Triangulation[connections_,labels_,faces_,periodicity_]["labels"[]]:=labels
Triangulation[connections_,labels_,faces_,periodicity_]["faces"[]]:=faces
Triangulation[connections_,labels_,faces_,periodicity_]["periodicity"[]]:=periodicity
Triangulation[connections_,labels_,faces_,periodicity_]["translations"[]]:=periodicity["translations"[]]/;PeriodicQ[periodicity]
Triangulation[connections_,labels_,faces_,periodicity_]["angles"[]]:=periodicity["angles"[]]/;ScrewQ[periodicity]
Triangulation[connections_,labels_,faces_,periodicity_]["axis"[]]:=periodicity["axis"[]]/;ScrewQ[periodicity]


(* ::Subsubsection::Closed:: *)
(*Parallelogram*)


Parallelogram[connections_,faces_,periodicity_]["connections"[]]:=connections
Parallelogram[connections_,faces_,periodicity_]["faces"[]]:=faces
Parallelogram[connections_,faces_,periodicity_]["periodicity"[]]:=periodicity
Parallelogram[connections_,faces_,periodicity_]["translations"[]]:=periodicity["translations"[]]/;PeriodicQ[periodicity]
Parallelogram[connections_,faces_,periodicity_]["angles"[]]:=periodicity["angles"[]]/;ScrewQ[periodicity]
Parallelogram[connections_,faces_,periodicity_]["axis"[]]:=periodicity["axis"[]]/;ScrewQ[periodicity]


(* ::Subsection::Closed:: *)
(*Surfaces*)


(* ::Subsubsection:: *)
(*plane*)


Plane[sheet_]["sheet"[]]:=sheet/;SpatialQ[sheet["periodicity"[]]]
Plane[sheet_]["normal"[]]:=Cross@sheet["translations"[]]/;SpatialQ[sheet["periodicity"[]]]
Plane[sheet_]["height"[]]:=Mean@sheet["vertices"[]]/;SpatialQ[sheet["periodicity"[]]]


(* ::Subsubsection:: *)
(*uniaxial*)


Uniaxial[sheet_]["sheet"[]]:=sheet/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["axis"[]]:=sheet["axis"[]]/;ScrewQ[sheet["periodicity"[]]]


Uniaxial[sheet_]["basis"[1]]:=-Normalize[RotationMatrix[(\[Pi]-First[#["angles"[]]])/2,#["axis"[]]].TransverseVec[First[#["translations"[]]],#["axis"[]]]]&@sheet/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["basis"[2]]:=Cross[Uniaxial[sheet]["basis"[3]],Uniaxial[sheet]["basis"[1]]]/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["basis"[3]]:=sheet["axis"[]]/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["basis"[]]:={#["basis"[1]],#["basis"[2]],#["basis"[3]]}&@Uniaxial[sheet]/;ScrewQ[sheet["periodicity"[]]]


Uniaxial[sheet_]["radius"[]]:=Norm[TransverseVec[First[#["translations"[]]],#["axis"[]]]]/(2 Sin[First[#["angles"[]]]/2])&@sheet/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["curvature"[]]:=1/Uniaxial[sheet]["radius"[]]/;ScrewQ[sheet["periodicity"[]]]
Uniaxial[sheet_]["center"[]]:=Uniaxial[sheet]["radius"[]]Normalize[RotationMatrix[(\[Pi]-First[#["angles"[]]])/2,#["axis"[]]].TransverseVec[First[#["translations"[]]],#["axis"[]]]]&@sheet/;ScrewQ[sheet["periodicity"[]]]


(* ::Subsection::Closed:: *)
(*Periodicity*)


(* ::Subsubsection::Closed:: *)
(*Spatial*)


Spatial[translations_]["translations"[]]:=translations


(* ::Subsubsection::Closed:: *)
(*Screw*)


Screw[translations_,angles_,axis_]["translations"[]]:=translations
Screw[translations_,angles_,axis_]["angles"[]]:=angles
Screw[translations_,angles_,axis_]["axis"[]]:=axis


(* ::Subsection::Closed:: *)
(*Linear Modes*)


(* ::Subsubsection::Closed:: *)
(*Folding motions*)


FoldingMotion[sheet_,motion_]["sheet"[]]:=sheet
FoldingMotion[sheet_,motion_]["motion"[]]:=motion


FoldingMotion[sheet_,motion_]["CellAngularVelocity"[path_]]:=Total@(Crease[sheet,First[#],Last[#]]["direction"[]]motion[[First[#]]]&/@path)


FoldingMotion[sheet_,motion_]["CellAngularVelocity"[path_,z1_,z2_]]:=Total@(Crease[sheet,First[#],Last[#]]["direction"[]]z1^First[Last[#]] z2^Last[Last[#]] motion[[First[#]]]&/@path)


FoldingMotion[sheet_,motion_]["CellDisplacement"[facePaths_,path_]]:=Total@(Crease[sheet,Last[#],{0,0}]["Rotation"[facePaths,First[#],motion]]&/@path)


FoldingMotion[sheet_,motion_]["CellDisplacement"[facePaths_,path_,cell_,z1_,z2_]]:=Total@(Crease[sheet,Last[#],cell]["Rotation"[facePaths,First[#],motion,z1,z2]]&/@path)


(* ::Section::Closed:: *)
(*queries*)


SheetQ[_Lattice]:=True
SheetQ[_Origami]:=True
SheetQ[_]:=False


NetworkQ[_Network]:=True
NetworkQ[_]:=False


OrigamiQ[_Origami]:=True
OrigamiQ[_]:=False


NodeQ[_Node]:=True
NodeQ[_]:=False


CreaseQ[_Crease]:=True
CreaseQ[_]:=False


FaceQ[_Face]:=True
FaceQ[_]:=False


TopologyQ[_Triangulation]:=True
TopologyQ[_Parallelogram]:=True
TopologyQ[_]:=False


TriangulationQ[_Triangulation]:=True
TriangulationQ[_]:=False


ParallelogramQ[_Parallelogram]:=True
ParallelogramQ[_]:=False


PeriodicQ[_Spatial]:=True
PeriodicQ[_Screw]:=True
PeriodicQ[_]:=False


SpatialQ[_Spatial]:=True
SpatialQ[_]:=False


ScrewQ[_Screw]:=True
ScrewQ[_]:=False


(* ::Section::Closed:: *)
(*computations*)


AxisTransform[sheet_]:=RotationMatrix[VectorAngle@@#,Normalize[Cross@@#]]&@{sheet["axis"[]],UnitVector[3,3]}


RotSum[\[Theta]_,N_]:={{1/2 (1-Cos[N \[Theta]]+Cot[\[Theta]/2] Sin[N \[Theta]]),Csc[\[Theta]/2] Sin[(N \[Theta])/2] Sin[1/2 (\[Theta]-N \[Theta])],0},{-Csc[\[Theta]/2] Sin[(N \[Theta])/2] Sin[1/2 (\[Theta]-N \[Theta])],1/2 (1-Cos[N \[Theta]]+Cot[\[Theta]/2] Sin[N \[Theta]]),0},{0,0,N}}/;\[Theta]!=0.
RotSum[\[Theta]_,N_]:=N IdentityMatrix[3]/;\[Theta]==0.
RotSum[\[Theta]_,N_,z_]:={{(1-z Cos[\[Theta]]+z^N (-Cos[N \[Theta]]+z Cos[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),(-z Sin[\[Theta]]+z^N (Sin[N \[Theta]]+z Sin[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),0},{(z Sin[\[Theta]]-z^N (Sin[N \[Theta]]+z Sin[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),(1-z Cos[\[Theta]]+z^N (-Cos[N \[Theta]]+z Cos[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),0},{0,0,(-1+z^N)/(-1+z)}}/;\[Theta]!=0.
RotSum[\[Theta]_,N_,z_]:=(-1+z^N)/(-1+z) IdentityMatrix[3]/;\[Theta]==0.


TransverseVec[vec_,axis_]:=vec-vec.axis axis


(* ::Subtitle:: *)
(*end*)


End[]
EndPackage[]
