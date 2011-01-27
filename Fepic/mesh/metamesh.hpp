#ifndef METAMESH_HPP
#define METAMESH_HPP

#include <fstream>
#include <iostream>

/*
	==================================================================
	========================== MeshMethods ===========================
	==================================================================
*/

template<class Mesh, int CellDim>
class MeshMethods;

template<class Mesh>
class MeshMethods<Mesh, 2>
{
public:
    template<class T>
	static void buildAdjacency(iMesh<T> & mesh)
	{
		mesh.buildAdjacency4face();
	}
    
    template<class T>
    static void remodelCellsNodes(iMesh<T> & mesh, int order)
    {
        mesh.remodelCellsNodes4face(order);
    }
    
    template<class T>
    static void buildCellLocalNodes(iMesh<T> & mesh)
    {
        mesh.edges_local_nodes = T::CellT::FaceT::getEdgesLocalNodes(mesh.getOrder());
    }
};

template<class Mesh>
class MeshMethods<Mesh, 3>
{
public:
    template<class Traits>
	static void buildAdjacency(iMesh<Traits> & mesh)
	{
		mesh.buildAdjacency4volume();
	}
    
    template<class T>
    static void remodelCellsNodes(iMesh<T> & mesh, int order)
    {
        mesh.remodelCellsNodes4volume(order);
    }
    
    template<class Traits>
    static void buildCellLocalNodes(iMesh<Traits> & mesh)
    {
        mesh.edges_local_nodes = Traits::CellT::getEdgesLocalNodes(mesh.getOrder());
        mesh.faces_local_nodes = Traits::CellT::getFacesLocalNodes(mesh.getOrder());
    }
};



/*
	==================================================================
	================== MeshReadMarkedElementsMSH =====================
	==================================================================
*/

template<class Mesh, int CellDim>
class MeshReadMarkedElementsMSH;

template<class Mesh>
class MeshReadMarkedElementsMSH<Mesh, 1>
{
public:
	static void buildAdjacency(std::ifstream &File)
	{
		Mesh::readMarkedElementsMSH4edge(File);
	}
};

template<class Mesh>
class MeshReadMarkedElementsMSH<Mesh, 2>
{
public:
	static void buildAdjacency(std::ifstream &File)
	{
		Mesh::readMarkedElementsMSH4face(File);
	}
};


template<class Mesh>
class MeshReadMarkedElementsMSH<Mesh, 3>
{
public:
	static void buildAdjacency(std::ifstream &File)
	{
		Mesh::readMarkedElementsMSH4volume(File);
	}
};


#endif // METAMESH_HPP

