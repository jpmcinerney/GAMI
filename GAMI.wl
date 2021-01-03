(* ::Package:: *)

(* ::Title:: *)
(*GAMI Package*)


BeginPackage["GAMI`"]


(* ::Section::Closed:: *)
(*Usage*)


GAMI::usage="Tools for the analysis of origami sheets."


NetworkQ::usage="Function to check whether an object is a network."


SpatialQ::usage="Function to check whether an object is spatially periodic."


ScrewQ::usage="Function to check whether an object is screw periodic."


SpatialQ::usage="Function to check whether an object is spatially periodic."


TriangulationQ::usage="Function to check whether an object is a triangulated crease pattern."


QuadrilateralQ::usage="Function to check whether an object is a quadrilateral crease pattern."


Network::usage="An object defined by vertex positions of the unit cell, spatial or screw periodicity, and bond topology."


Spatial::usage="A periodicity object defined by lattice vectors."


Screw::usage="A periodicity object defined by lattice vectors, lattice rotation angles, and a rotation axis."


Triangulation::usage="A topology object defined by edge connections and triangular faces."


Quadrilateral::usage="A topology object defined by edge connections and quadrilateral faces."


Vertex::usage="An object defined by a particular vertex in a sheet."


Edge::usage="An object defined by a particular edge in a sheet."


Face::usage="An object defined by a particular face in a sheet."


Unit::usage="An object defined by the unit cell of a sheet."


Tessellation::usage="The tessellation of a sheet with a chosen number of cells."


RigidityMatrix::usage="The rigidity matrix for a network."


EquilibriumMatrix::usage="The equilibrium matrix for a network."


AxisTransform::usage=
"Constructs the similarity transform rotating the sheet's screw axis to the z-axis."


RotSum::usage=
"Sums rotation matrices about the z-axis."


(* ::Section:: *)
(*Errors*)


(* ::Section:: *)
(*Routines*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Network*)


NetworkQ[_Network]:=True
NetworkQ[_]:=False


(* ::Text:: *)
(* *)


Network[vertices_,periodicity_,topology_]["Vertices"]:=vertices
Network[vertices_,periodicity_,topology_]["Vertex"[index_]]:=vertices[[index]]


(* ::Text:: *)
(* *)


Network[vertices_,periodicity_,topology_]["Periodicity"]:=periodicity
Network[vertices_,periodicity_,topology_]["Translations"]:=periodicity["Translations"]
Network[vertices_,periodicity_,topology_]["Translation"[index_]]:=periodicity["Translation"[index]]
Network[vertices_,periodicity_,topology_]["Angles"]:=periodicity["Angles"]/;ScrewQ[periodicity]
Network[vertices_,periodicity_,topology_]["Angle"[index_]]:=periodicity["Angle"[index]]/;ScrewQ[periodicity]
Network[vertices_,periodicity_,topology_]["Axis"]:=periodicity["Axis"]/;ScrewQ[periodicity]


(* ::Text:: *)
(* *)


Network[vertices_,periodicity_,topology_]["Topology"]:=topology
Network[vertices_,periodicity_,topology_]["Edges"]:=topology["Edges"]
Network[vertices_,periodicity_,topology_]["Edge"[index_]]:=topology["Edge"[index]]
Network[vertices_,periodicity_,topology_]["Faces"]:=topology["Faces"]
Network[vertices_,periodicity_,topology_]["Face"[index_]]:=topology["Face"[index]]
Network[vertices_,periodicity_,topology_]["NumEdges"]:=Length[topology["Edges"]]
Network[vertices_,periodicity_,topology_]["NumFaces"]:=Length[topology["Faces"]]


(* ::Subsection::Closed:: *)
(*Periodicity*)


SpatialQ[_Spatial]:=True
SpatialQ[network_Network]:=SpatialQ@network["Periodicity"]
SpatialQ[_]:=False
ScrewQ[_Screw]:=True
ScrewQ[network_Network]:=ScrewQ@network["Periodicity"]
ScrewQ[_]:=False


(* ::Text:: *)
(* *)


Spatial[translations_]["Translations"]:=translations
Spatial[translations_]["Translation"[index_]]:=translations[[index]]


(* ::Text:: *)
(* *)


Screw[translations_,angles_,axis_]["Translations"]:=translations
Screw[translations_,angles_,axis_]["Translation"[index_]]:=translations[[index]]
Screw[translations_,angles_,axis_]["Angles"]:=angles
Screw[translations_,angles_,axis_]["Angle"[index_]]:=angles[[index]]
Screw[translations_,angles_,axis_]["Axis"]:=axis


(* ::Subsection::Closed:: *)
(*Topology*)


TriangulationQ[_Triangulation]:=True
TriangulationQ[network_Network]:=TriangulationQ@Network["Topology"]
TriangulationQ[_]:=False
QuadrilateralQ[_Quadrilateral]:=True
QuadrilateralQ[network_Network]:=QuadrilateralQ@Network["Topology"]
QuadrilateralQ[_]:=False


(* ::Text:: *)
(* *)


Triangulation[edges_,faces_]["Edges"]:=edges
Triangulation[edges_,faces_]["Edge"[index_]]:=edges[[index]]
Triangulation[edges_,faces_]["Faces"]:=faces
Triangulation[edges_,faces_]["Face"[index_]]:=faces[[index]]


(* ::Text:: *)
(* *)


Quadrilateral[edges_,faces_]["Edges"]:=edges
Quadrilateral[edges_,faces_]["Edge"[index_]]:=edges[[index]]
Quadrilateral[edges_,faces_]["Faces"]:=faces
Quadrilateral[edges_,faces_]["Face"[index_]]:=faces[[index]]


(* ::Subsection::Closed:: *)
(*Vertex*)


Vertex[network_,index_,cell_:{0,0}]["Position"]:=#["Vertex"[index]]+cell.#["Translations"]&@network/;SpatialQ[network["Periodicity"]]


AxisTransform[sheet_]:=RotationMatrix[VectorAngle@@#,Normalize[Cross@@#]]&@{sheet["Axis"],UnitVector[3,3]}


RotSum[\[Theta]_,N_]:={{1/2 (1-Cos[N \[Theta]]+Cot[\[Theta]/2] Sin[N \[Theta]]),Csc[\[Theta]/2] Sin[(N \[Theta])/2] Sin[1/2 (\[Theta]-N \[Theta])],0},{-Csc[\[Theta]/2] Sin[(N \[Theta])/2] Sin[1/2 (\[Theta]-N \[Theta])],1/2 (1-Cos[N \[Theta]]+Cot[\[Theta]/2] Sin[N \[Theta]]),0},{0,0,N}}/;\[Theta]!=0.
RotSum[\[Theta]_,N_]:=N IdentityMatrix[3]/;\[Theta]==0.
RotSum[\[Theta]_,N_,z_]:={{(1-z Cos[\[Theta]]+z^N (-Cos[N \[Theta]]+z Cos[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),(-z Sin[\[Theta]]+z^N (Sin[N \[Theta]]+z Sin[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),0},{(z Sin[\[Theta]]-z^N (Sin[N \[Theta]]+z Sin[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),(1-z Cos[\[Theta]]+z^N (-Cos[N \[Theta]]+z Cos[\[Theta]-N \[Theta]]))/(1+z^2-2 z Cos[\[Theta]]),0},{0,0,(-1+z^N)/(-1+z)}}/;\[Theta]!=0.
RotSum[\[Theta]_,N_,z_]:=(-1+z^N)/(-1+z) IdentityMatrix[3]/;\[Theta]==0.


(*Vertex[network_,index_,cell_:{0,0}]["Position"]:=rotationSum@{#["Angle"[1]],#["Axis"],cell[[1]]}.#["Translation"[1]]+
RotationMatrix[cell[[1]]#["Angle"[1]],#["Axis"]].rotationSum@{#["Angle"[2]],#["Axis"],cell[[2]]}.#["Translation"[2]]+
RotationMatrix[cell. #["Angles"],#["Axis"]].#["Vertex"[index]]&@network/;ScrewQ[network["Periodicity"]]*)


Vertex[network_,index_,cell_:{0,0}]["Position"]:=With[{
trans=AxisTransform[network],
rots=RotationMatrix[First[#]Last[#],network["Axis"]]&/@Transpose[{network["Angles"],cell}]
},
(Dot@@rots).#["Vertex"[index]]+Inverse[trans].RotSum[First[#["Angles"]],First[cell]].trans.First[#["Translations"]]+First[rots].Inverse[trans].RotSum[Last[#["Angles"]],Last[cell]].trans.Last[#["Translations"]]&@network
]/;ScrewQ[network["Periodicity"]]


Vertex[network_,index_,cell_:{0,0}]["Complex"]:=Point[Vertex[network,index,cell]["Position"]]


(*rotationSum={{-(1/(2 tempX^2 (tempX^2+tempY^2)))(-1+tempZ^2) (-1+tempY^2+tempZ^2) (-2 tempCellIndex+(-1+2 tempCellIndex) tempY^2+(-1+2 tempCellIndex) tempZ^2+(tempY^2+tempZ^2) Cos[tempCellIndex tempLatticeRotationAngle]-(tempY^2+tempZ^2) Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle]),1/(2 tempX) ((-1+2 tempCellIndex) tempY+(1-2 tempCellIndex) tempY^3+(1-2 tempCellIndex) tempY tempZ^2+Cos[tempCellIndex tempLatticeRotationAngle] (-tempY (-1+tempY^2+tempZ^2)+tempX tempZ Cot[tempLatticeRotationAngle/2])+tempX tempZ Sin[tempCellIndex tempLatticeRotationAngle]+Cot[tempLatticeRotationAngle/2] (-tempX tempZ+tempY (-1+tempY^2+tempZ^2) Sin[tempCellIndex tempLatticeRotationAngle])),1/2 (-tempY (-1+Cos[tempCellIndex tempLatticeRotationAngle]) Cot[tempLatticeRotationAngle/2]-tempY Sin[tempCellIndex tempLatticeRotationAngle]+tempX tempZ (-1+2 tempCellIndex+Cos[tempCellIndex tempLatticeRotationAngle]-Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle]))},{-(1/(2 tempX))(tempY-2 tempCellIndex tempY+(-1+2 tempCellIndex) tempY^3+(-1+2 tempCellIndex) tempY tempZ^2+tempY (-1+tempY^2+tempZ^2) Cos[tempCellIndex tempLatticeRotationAngle]-tempY (-1+tempY^2+tempZ^2) Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle]+2 tempX tempZ Csc[tempLatticeRotationAngle/2] Sin[(tempCellIndex tempLatticeRotationAngle)/2] Sin[1/2 (tempLatticeRotationAngle-tempCellIndex tempLatticeRotationAngle)]),1/2 (1+(-1+2 tempCellIndex) tempY^2+(-1+tempY^2) Cos[tempCellIndex tempLatticeRotationAngle]-(-1+tempY^2) Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle]),1/2 (tempX ((-1+Cos[tempCellIndex tempLatticeRotationAngle]) Cot[tempLatticeRotationAngle/2]+Sin[tempCellIndex tempLatticeRotationAngle])+tempY tempZ (-1+2 tempCellIndex+Cos[tempCellIndex tempLatticeRotationAngle]-Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle]))},{-(1/(2 tempX^2))(-1+tempY^2+tempZ^2) (tempY ((-1+Cos[tempCellIndex tempLatticeRotationAngle]) Cot[tempLatticeRotationAngle/2]+Sin[tempCellIndex tempLatticeRotationAngle])+tempX tempZ (-1+2 tempCellIndex+Cos[tempCellIndex tempLatticeRotationAngle]-Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle])),1/(2 tempX) (Cot[tempLatticeRotationAngle/2]-Cos[tempLatticeRotationAngle/2-tempCellIndex tempLatticeRotationAngle] Csc[tempLatticeRotationAngle/2]+tempX tempY tempZ (-1+2 tempCellIndex+Cos[tempCellIndex tempLatticeRotationAngle]-Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle])+2 tempY^2 Csc[tempLatticeRotationAngle/2] Sin[(tempCellIndex tempLatticeRotationAngle)/2] Sin[1/2 (tempLatticeRotationAngle-tempCellIndex tempLatticeRotationAngle)]+2 tempZ^2 Csc[tempLatticeRotationAngle/2] Sin[(tempCellIndex tempLatticeRotationAngle)/2] Sin[1/2 (tempLatticeRotationAngle-tempCellIndex tempLatticeRotationAngle)]),1/2 (1+(-1+2 tempCellIndex) tempZ^2+(-1+tempZ^2) Cos[tempCellIndex tempLatticeRotationAngle]-(-1+tempZ^2) Cot[tempLatticeRotationAngle/2] Sin[tempCellIndex tempLatticeRotationAngle])}}/.{tempLatticeRotationAngle->#[[1]],tempX->#[[2,1]],tempY->#[[2,2]],tempZ->#[[2,3]],tempCellIndex->#[[3]]}&;*)


(* ::Subsection::Closed:: *)
(*Edge*)


Edge[network_,index_,cell_:{0,0}]["Vector"]:=Flatten[Differences[Vertex[network,First[#],Last[#]+cell]["Position"]&/@network["Edge"[index]]]]


Edge[network_,index_,cell_:{0,0}]["Direction"]:=Normalize[Edge[network,index,cell]["Vector"]]


Edge[network_,index_,cell_:{0,0}]["Length"]:=Norm[Edge[network,index,cell]["Vector"]]


Edge[network_,index_,cell_:{0,0}]["Complex"]:=Line[Vertex[network,First[#],cell+Last[#]]["Position"]&/@network["Edge"[index]]]


(* ::Subsection::Closed:: *)
(*Face*)


Face[network_,index_,cell_:{0,0}]["Complex"]:=Polygon[Vertex[network,First[#],Last[#]+cell]["Position"]&/@network["Face"[index]]]


(* ::Subsection::Closed:: *)
(*Unit Cell*)


Unit[network_,cell_:{0,0}]["Complex"]:=Table[Face[network,index,cell]["Complex"],{index,Length@network["Faces"]}]


(* ::Subsection::Closed:: *)
(*Tessellation*)


Tessellation[network_,cells_]["Complex"]:=Flatten[Table[Unit[network,{n1,n2}]["Complex"],{n1,0,First[cells]-1},{n2,0,Last[cells]-1}]]


(* ::Subsection::Closed:: *)
(*Rigidity Matrix*)


(*RigidityMatrix[network_]:=Normal@SparseArray@Flatten@Table[
{
{index,3 network["Edge"[index]][[1,1]]-2}->#[[1]],
{index,3 network["Edge"[index]][[1,1]]-1}->#[[2]],
{index,3 network["Edge"[index]][[1,1]]}->#[[3]],
{index,3 network["Edge"[index]][[2,1]]-2}->-#[[1]],
{index,3 network["Edge"[index]][[2,1]]-1}->-#[[2]],
{index,3 network["Edge"[index]][[2,1]]}->-#[[3]]
}&@Edge[network,index]["Direction"],{index,network["NumEdges"]}]/;SpatialQ@network*)


RigidityMatrix[network_,bloch_:{1.,1.}]:=Normal@SparseArray@Flatten@Table[
{
{index,3 network["Edge"[index]][[1,1]]-2}->#[[1,1]] #[[2,2]],
{index,3 network["Edge"[index]][[1,1]]-1}->#[[1,2]]#[[2,2]],
{index,3 network["Edge"[index]][[1,1]]}->#[[1,3]]#[[2,2]],
{index,3 network["Edge"[index]][[2,1]]-2}->-#[[1,1]]#[[2,1]],
{index,3 network["Edge"[index]][[2,1]]-1}->-#[[1,2]]#[[2,1]],
{index,3 network["Edge"[index]][[2,1]]}->-#[[1,3]]#[[2,1]]
}&@{Edge[network,index]["Direction"],Times@@Thread[Power[bloch,Last[#]]]&/@network["Edges"][[index]]},{index,network["NumEdges"]}]/;SpatialQ@network


(* ::Subsection::Closed:: *)
(*Equilibrium Matrix*)


EquilibriumMatrix[network_,bloch_:{1.,1.}]:=Normal@SparseArray@Flatten@Table[
{
{3 network["Edge"[index]][[1,1]]-2,index}->#[[1,1]]#[[2,2]],
{3 network["Edge"[index]][[1,1]]-1,index}->#[[1,2]]#[[2,2]],
{3 network["Edge"[index]][[1,1]],index}->#[[1,3]]#[[2,2]],
{3 network["Edge"[index]][[2,1]]-2,index}->-#[[1,1]]#[[2,1]],
{3 network["Edge"[index]][[2,1]]-1,index}->-#[[1,2]]#[[2,1]],
{3 network["Edge"[index]][[2,1]],index}->-#[[1,3]]#[[2,1]]
}&@{Edge[network,index]["Direction"],Times@@Thread[Power[bloch,-Last[#]]]&/@network["Edges"][[index]]},{index,network["NumEdges"]}]/;SpatialQ@network


(* ::Title:: *)
(*End Package*)


End[]
EndPackage[]
