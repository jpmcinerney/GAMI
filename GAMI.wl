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


Vertex[network_,index_,cell_:{0,0}]["Complex"]:=Point[Vertex[network,index,cell]["Position"]]


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


(* ::Title:: *)
(*End Package*)


End[]
EndPackage[]
