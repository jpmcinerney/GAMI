(* ::Package:: *)

(* ::Title:: *)
(*GAMI Package*)


BeginPackage["GAMI`"]


(* ::Section:: *)
(*Usage*)


(* ::Section:: *)
(*Errors*)


(* ::Section:: *)
(*Routines*)


Begin["`Private`"]


(* ::Subsection:: *)
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


Network[vertices_,periodiciity_,topology_]["Topology"]:=topology
Network[vertices_,periodiciity_,topology_]["Edges"]:=topology["Edges"]
Network[vertices_,periodiciity_,topology_]["Edge"[index_]]:=topology["Edge"[index]]
Network[vertices_,periodiciity_,topology_]["Faces"]:=topology["Faces"]


(* ::Subsection::Closed:: *)
(*Periodicity*)


SpatialQ[_Spatial]:=True
SpatialQ[_]:=False
ScrewQ[_Screw]:=True
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
TriangulationQ[_]:=False
QuadrilateralQ[_Quadrilateral]:=True
QuadrilateralQ[_]:=False


(* ::Text:: *)
(* *)


Triangulation[edges_,faces_]["Edges"]:=edges
Triangulation[edges_,faces_]["Edge"[index_]]:=edges[index]
Triangulation[edges_,faces_]["Faces"]:=faces
Triangulation[edges_,faces_]["Face"[index_]]:=faces[index]


(* ::Text:: *)
(* *)


Quadrilateral[edges_,faces_]["Edges"]:=edges
Quadrilateral[edges_,faces_]["Edge"[index_]]:=edges[index]
Quadrilateral[edges_,faces_]["Faces"]:=faces
Quadrilateral[edges_,faces_]["Face"[index_]]:=faces[index]


(* ::Subsection:: *)
(*Vertex*)


Vertex[network_,index_,cell_:{0,0}]["Position"]:=#["Vertex"[index]]+cell.#["Translations"]&@network/;SpatialQ[network["Periodicity"]]


Vertex[network_,index_,cell_:{0,0}]["Complex"]:=Point[Vertex[network,index,cell]["Complex"]]


(* ::Subsection:: *)
(*Edge*)


Edge[network_,index_,cell_:{0,0}]["Vector"]:=Flatten[Differences[Vertex[network,First[#],Last[#]+cell]["Position"]&/@Extract[network["Edge"[index]]]]]


Edge[network_,index_,cell_:{0,0}]["Direction"]:=Normalize[Edge[network,index,cell]["Vector"]]


Edge[network_,index_,cell_:{0,0}]["Length"]:=Norm[Edge[network,index,cell]["Vector"]]


Edge[network_,index_,cell_:{0,0}]["Complex"]:=Line[Vertex[network,First[#],cell+Last[#]]["Position"]&/@network["Edge"[index]]]


(* ::Subsection:: *)
(*Face*)


Face[network_,index_,cell_:{0,0}]["Complex"]:=Polygon[Vertex[network,First[#],Last[#]+cell]["Position"]&/@network["Face"[index]]]


(* ::Subsection:: *)
(*Unit Cell*)


Unit[network_,cell_:{0,0}]["Complex"]:=Table[Face[network,index]["Complex"],{index,Length@network["Faces"]}]


(* ::Subsection:: *)
(*Tessellation*)


Tessellation[network_,cells_]["Complex"]:=Flatten[Table[Unit[network,{n1,n2}]["Complex"],{n1,0,First[cells]-1},{n2,0,Last[cells]-1}]]


(* ::Title::Closed:: *)
(*End Package*)


End[]
EndPackage[]
