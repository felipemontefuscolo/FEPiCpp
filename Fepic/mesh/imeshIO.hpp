#ifndef IMESHIO_HPP
#define IMESHIO_HPP

#include "imesh.hpp"
#include "../misc/misc.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>


/* Este arquivo contém iterações em células. Para melhorar o desempenho deve-se
 * refaze-los com iteradores.
 * 
 * - na função buildAdjacency()
 * 
 * Este arquivo contém:
 * - TERMINAR
 * - PENSAR EM ALGO MELHOR
 * 
 * 
 * */


/*
 * MUDANÇA: PARA MELHOR O DESEMPENHO, A FUNÇÃO PODE SER REFEITA USANDO ITERADORES
 * EM VEZ DE ACESSO RANDOMICO.
 */ 
/** Uma vez já lida uma tabela de conectividade, i.e., as células
 *  e seus nós, a função abaixo constrói as relações de adjacências
 *  (half-edfe, half-face ) ou nenhum se a célula for aresta.
 * @note células mortas não são cosideradas.
 */
template<class Traits>
void iMesh<Traits>::buildAdjacency4face()
{

	int n_borders = CellT::N_borders;
	Fepic::vectorui edge_vtx(2);
	Fepic::vectorui cell_ith(2);
	std::map<Fepic::vectorui, Fepic::vectorui> table; // Key: (n0, n1) ... atributos (elemento, ith)
	std::map<Fepic::vectorui, Fepic::vectorui>::iterator tab_it, tab_it2;
	CellIterator cell = iMesh::_cellL.begin();
	uint otherC, otherith, thisC, thisith;
	int thislabel;
	
	
	
	uint k=0;
	for (auto cellend=_cellL.end(); cell != cellend ; ++cell)
	{
		if (!cell->isDead())
		{
			for (int ith = 0; ith < n_borders; ++ith)
			{
				edge_vtx = cell->getBorderVertices(ith);
				cell_ith[0] = k;
				cell_ith[1] = ith;
				table.insert(std::pair<Fepic::vectorui, Fepic::vectorui>(edge_vtx, cell_ith));
			}
		}
		++k;
	}
	
	while (!table.empty())
	{
		tab_it = --table.end();
		reverse_copy(tab_it->first.begin(), tab_it->first.end(), edge_vtx.begin());
		tab_it2 = table.find(edge_vtx);
		if (tab_it2 != table.end()) // se econtrou uma edge em comum
		{
			otherC = tab_it2->second[0];
			otherith = tab_it2->second[1];
			thisC = tab_it->second[0];
			thisith =  tab_it->second[1];
			
			this->getCell(thisC)->getHalf(thisith)->setCompleteID(otherC, otherith);
			this->getCell(otherC)->getHalf(otherith)->setCompleteID(thisC, thisith);
			
			table.erase(tab_it2);
			table.erase(tab_it);
		}
		else
		{
			thisC = tab_it->second[0];
			thisith =  tab_it->second[1];
			thislabel = this->getCell(thisC)->getLabel();
			MHalfEdgeT temp(thisC, thisith, thislabel);
			uint mHE_universal_iD = addMHalf(temp);
			this->getCell(thisC)->getHalf(thisith)->setCompleteID(mHE_universal_iD, -1);
			table.erase(tab_it);	
		}
	}
	
}


/** @note células mortas não são consideradas.
 */ 
template<class Traits>
void iMesh<Traits>::buildAdjacency4volume()
{

	//const int n_borders = CellT::N_borders;
	//const int n_anch    = FaceT::N_borders;
    
    //typedef std::pair<std::array<uint, n_anch>, std::array<uint, 2> > PairT; // < (vértices da face) , (elemento, ith) >
    //typedef std::unordered_map<std::array<uint, n_anch>, std::array<uint, 2> >  MapT; 
    
	//std::array<uint, n_anch>  face_vtx;
	//std::array<uint, 2>       cell_ith;
	//MapT                      table; // < (vértices da face) , (elemento, ith) >
	//auto                      tab_it = table.begin(), tab_it2 = table.begin();  // iterators
	//CellIterator              cell = this->_cellL.begin();
	//uint                      otherC, otherith, thisC, thisith;
	//int                       thislabel;
    //int                       a, k=0;
    //bool                      found;
	
	//for (auto cellend=_cellL.end(); cell != cellend ; ++cell)
	//{
		//if (!cell->isDead())
		//{
			//for (int ith = 0; ith < n_borders; ++ith)
			//{
				////face_vtx = cell->getBorderVertices(ith);
                //copy_from(cell->getBorderVertices(ith).begin(), face_vtx.begin(), face_vtx.end());
				//cell_ith[0] = k;
				//cell_ith[1] = ith;
				//table.insert(PairT(face_vtx, cell_ith));
			//}
		//}
		//++k;
	//}
	
	//while (!table.empty())
	//{
        //tab_it = table.begin();
		//reverse_copy(tab_it->first.begin(), tab_it->first.end(), face_vtx.begin());
		//found = false;
		//for (a = 0; a != n_anch; ++a) // ancora
		//{
			//tab_it2 = table.find(face_vtx);
			//if (tab_it2 != table.end()) // se econtrou uma face em comum
			//{
				//otherC = tab_it2->second[0];
				//otherith = tab_it2->second[1];
				//thisC = tab_it->second[0];
				//thisith =  tab_it->second[1];
				
				//this->getCell(thisC)->getHalf(thisith)->setCompleteID(otherC, otherith, a);
				//this->getCell(otherC)->getHalf(otherith)->setCompleteID(thisC, thisith, a);
				
				//table.erase(tab_it2);
				//table.erase(tab_it);
				//found = true;
				//break;
			//}
			//rotate(face_vtx.begin(), face_vtx.begin()+1, face_vtx.end());
		//}
		//if(!found)
		//{
			//thisC = tab_it->second[0];
			//thisith =  tab_it->second[1];
			//thislabel = this->getCell(thisC)->getLabel();
			//MHalfFaceT temp(thisC, thisith, 0, thislabel);
			//uint mHE_universal_iD = addMHalf(temp);
			//this->getCell(thisC)->getHalf(thisith)->setCompleteID(mHE_universal_iD, -1, 0);
			//table.erase(tab_it);	
		//}
		
	//}


	const int n_borders = CellT::N_borders;
	const int n_anch    = FaceT::N_borders;
    
    typedef std::pair<Fepic::vectorui, Fepic::vectorui>   PairT; // < (vértices da face) , (elemento, ith) >
    typedef std::map<Fepic::vectorui, Fepic::vectorui>  MapT;
     
	Fepic::vectorui   face_vtx(FaceT::N_borders);
	Fepic::vectorui   cell_ith(2);
	MapT            table; // < (nós da face) , (elemento, ith) >
	MapT::iterator  tab_it, tab_it2;
	CellIterator    cell = this->_cellL.begin();
	uint            otherC, otherith, thisC, thisith;
	int             thislabel;
    int             a;
	
	uint k=0;
	for (auto cellend=_cellL.end(); cell != cellend ; ++cell)
	{
		if (!cell->isDead())
		{
			for (int ith = 0; ith != n_borders; ++ith)
			{
				face_vtx = std::move( cell->getBorderVertices(ith) );
				cell_ith[0] = k;
				cell_ith[1] = ith;
				table.insert(std::make_pair(face_vtx, cell_ith));
			}
		}
		++k;
	}
	
	bool found;
	while (!table.empty())
	{
        tab_it = table.begin();
		reverse_copy(tab_it->first.begin(), tab_it->first.end(), face_vtx.begin());
		found = false;
		for (a = 0; a != n_anch; ++a) // ancora
		{
			tab_it2 = table.find(face_vtx);
			if (tab_it2 != table.end()) // se econtrou uma face em comum
			{
				otherC = tab_it2->second[0];
				otherith = tab_it2->second[1];
				thisC = tab_it->second[0];
				thisith =  tab_it->second[1];
				
				_cellL[thisC].getHalf(thisith)->setCompleteID(otherC, otherith, a);
				_cellL[otherC].getHalf(otherith)->setCompleteID(thisC, thisith, a);
				
				table.erase(tab_it2);
				table.erase(tab_it);
				found = true;
				break;
			}
			rotate(face_vtx.begin(), face_vtx.begin()+1, face_vtx.end());
		}
		if(!found)
		{
			thisC = tab_it->second[0];
			thisith =  tab_it->second[1];
			thislabel = this->getCell(thisC)->getLabel();
			MHalfFaceT temp(thisC, thisith, 0, thislabel);
			uint mHE_universal_iD = addMHalf(temp);
			this->getCell(thisC)->getHalf(thisith)->setCompleteID(mHE_universal_iD, -1, 0);
			table.erase(tab_it);	
		}
		
	}

	
}

/** Lê a malha de um arquivo no formato .msh. A ordem da malha é detectada
 * altomaticamente.
 */ 
template<class Traits>
void iMesh<Traits>::readFileMSH(const char* filename)
{

	/*
	 * O que é feito:
	 * 
	 * 1) primeiramente lê-se apenas as células (conectividade clássica);
	 * 2) depois são construídos as half-edges e(ou) half-faces, lendo-se
	 * 	  os elementos de dimensões menores.
	 * 
	 * */

    std::ifstream File(filename);
    
    if (File.fail()) {
        std::cout << "I can't read the file "<< filename << " ... " << std::endl;
        throw;
    };   
    _meshfile = filename;
	int _final = _meshfile.find(".msh");
	_basename = _meshfile.substr(0,_final);
    
    double	  coord[3];
    int		  type_aux;
    char	  lixo[256];
    
    size_t    NODE_POS;  // $Nodes
	size_t    ELEM_POS;  // $Elements
    
    // nós
    NODE_POS = search_word(File, "$Nodes");

	uint num_pts;
	File >> num_pts;	
	_pointL.resize(num_pts);
	
	for (uint i=0; i< getNumNodes(); ++i) {
        File >> lixo;  // nº do nó
        File >> coord[0];
        File >> coord[1];
        File >> coord[2];
        _pointL[i].setCoord(coord);
	}
	// os pontos não estão completas: falta atribuir os labels
	
	// contagem de elementos e alocação
    ELEM_POS = search_word(File, "$Elements");
    
    uint num_cells=0;
    uint num_elms;
    File >> num_elms;

	/* ---------------------------------------
	 * Detectando a ordem da malha 
	 * --------------------------------------- */
	for (uint k = 0; k < num_elms; ++k)
	{
		File >> lixo; // numero do elemento
		File >> type_aux;
		int elm_dim = getDimForElementTypeMSH(type_aux);
		if (elm_dim == CellT::Dim)
		{
			this->setOrder( getOrderForElementTypeMSH(type_aux) );
            break;
		}
		
		File.getline(lixo, 256);
	}
	int order = this->_order;
    
    
	CellT	  Cell;
	/* --------------------------------------
	 * Lendo as células
	 * -------------------------------------- */
	File.seekg(ELEM_POS);
	File >> num_elms;
    Cell.setOrder(order);

	uint id_aux;
	int  num_tags;
	int  label_aux;
	uint counter; // contagem para verificação de erros.
	for (uint k=0; k < num_elms; ++k)
	{
		
		File >> counter; // numero do elemento .. verificação de sincronia
		if (counter != k+1)
		{
			std::cout << filename <<  ": invalid file format!" << std::endl;
			throw;
		}
		
		File >> type_aux;
		int elm_dim = getDimForElementTypeMSH(type_aux);
		
		File >> num_tags;
		
		File >> label_aux;
		
		for (int j=0; j<num_tags-1; ++j)
		{
			File >> lixo;
		}
		
		switch (elm_dim)
		{
			case 0: // caso ponto
				File >> id_aux;
				--id_aux;
				this->getNode(id_aux)->setLabel(label_aux);
				break;
			
			case 1: // caso aresta
				{
					if (M_Compare<EdgeT, CellT>::isEqual)
					{
						if (CellT::getMSHTag(order) != type_aux)	
						{
							std::cout << "Invalid edge! file has: " << getElementNameMSH(type_aux) << std::endl;
							std::cout << "But Mesh is configured to support: "
                                      << getElementNameMSH(EdgeT::getMSHTag(order)) << std::endl;
							throw;
						}
						for (int i=0; i<Cell.getNumNodes(); ++i)
						{
							uint nodeid;
							File >> nodeid;
							--nodeid;
							Cell.setNode(i, nodeid);
						}
						Cell.setLabel(label_aux);
						this->addCell(Cell);
						++num_cells;
					}
					else
						File.getline(lixo, 256);
				}
				break;
				
			case 2: // caso face
				{
					if (M_Compare<FaceT, CellT>::isEqual)
					{
						if (CellT::getMSHTag(order) != type_aux)	
						{
							std::cout << "Invalid face! file has: " << getElementNameMSH(type_aux) << std::endl;
							std::cout << "But Mesh is configured to support: "
                                      << getElementNameMSH(FaceT::getMSHTag(order)) << std::endl;
							throw;
						}
						for (int i=0; i<Cell.getNumNodes(); ++i)
						{
							uint nodeid;
							File >> nodeid;
							--nodeid;
							Cell.setNode(i, nodeid);
						}
						Cell.setLabel(label_aux);
						this->addCell(Cell);
						++num_cells;
					}
					else
						File.getline(lixo, 256);
				}
				break;

			case 3: // caso volume
				{
					if (M_Compare<VolumeT, CellT>::isEqual)
					{
						if (CellT::getMSHTag(order) != type_aux)	
						{
							std::cout << "Invalid volume! file has: " << getElementNameMSH(type_aux) << std::endl;
							std::cout << "But Mesh is configured to support: "
                                      << getElementNameMSH(FaceT::getMSHTag(order)) << std::endl;
							throw;
						}
						for (int i=0; i<Cell.getNumNodes(); ++i)
						{
							uint nodeid;
							File >> nodeid;
							--nodeid;
							Cell.setNode(i, nodeid);
						}
						Cell.setLabel(label_aux);
						this->addCell(Cell);
						++num_cells;
					}
					else
						File.getline(lixo, 256);
				}
				break;				
				
			default:
				std::cout << "invalid element ...\n";
				throw;
		}
		
	}
	
	if (num_cells == 0)
		std::cout << "WARNING: mesh file don't have any "
                  << getElementNameMSH(CellT::getMSHTag(order)) << std::endl;
	
	/* até aqui, apenas foi lido a conectividade
	 */ 
	
	/* constroi as halfs e Mhalfs */
	MeshMethods<iMesh<Traits>, CellT::Dim>::buildAdjacency(*this);

	
	/*
	___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__
	_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___| */
	
	/* Procura por elementos no contorno que coincidam com as halfs. Quando encontrados,
	 * os labels são propagados para as respectivas half que por sua vez propagam os
	 * labels para seus nós. Os labels só são propagados se eles não foram definidos
	 * nos nós anteriormente.
	 */
	File.seekg(ELEM_POS);
	File >> num_elms;
	
	
	for (uint k=0; k < num_elms; ++k)
	{
		
		File >> counter; // numero do elemento
		if (counter != k+1)
		{
			std::cout << filename <<  ": invalid file format!" << std::endl;
			throw;
		}
		
		File >> type_aux;
		int elm_dim = getDimForElementTypeMSH(type_aux);
		
		File >> num_tags;
		
		File >> label_aux;
		
		for (int j=0; j<num_tags-1; ++j)
		{
			File >> lixo;
		}
		
		uint id_aux;
		switch (elm_dim)
		{
			case 0: // caso ponto
				File >> id_aux;
				break;
			
			case 1: // caso aresta
				if (M_Compare<EdgeT, CellBT>::isEqual)
				{
					Fepic::vectorui nodes(CellBT::getNumNodesCell(order));
					Fepic::vectorui vtx(2);
					uint half_id;
					for (int i=0; i<CellBT::getNumNodesCell(order); ++i)
					{
						File >> nodes[i];
						--nodes[i];
						if (this->getNode(nodes[i])->getLabel() == 0)
						{
							this->getNode(nodes[i])->setLabel(label_aux);
						}
					}
					copy( nodes.begin(), nodes.begin()+2, vtx.begin() );
					if( this->theseVerticesFormAMHalf(vtx, half_id) )
					{
						this->getMHalf(half_id)->setLabel(label_aux);
					}
				}
				else
				{
					uint nodeid;
					for (int i=0; i<getNumVerticesForElementTypeMSH(type_aux); ++i)
					{
						File >> nodeid;
						--nodeid;
						if (this->getNode(nodeid)->getLabel() == 0)
						{
							this->getNode(nodeid)->setLabel(label_aux);
						}
					}
				}
				break;
				
			case 2: // caso face
				if (M_Compare<FaceT, CellBT>::isEqual)
				{
					Fepic::vectorui nodes(CellBT::getNumNodesCell(order));
					Fepic::vectorui vtx(CellBT::N_borders);
					uint half_id=0;
					for (int i=0; i<CellBT::getNumNodesCell(order); ++i)
					{
						File >> nodes[i];
						--nodes[i];
						if (this->getNode(nodes[i])->getLabel() == 0)
						{
							this->getNode(nodes[i])->setLabel(label_aux);
						}
					}
					copy( nodes.begin(), nodes.begin()+CellBT::N_borders, vtx.begin() );
					if( this->theseVerticesFormAMHalf(vtx, half_id) )
					{
						this->getMHalf(half_id)->setLabel(label_aux); //std::cout << (++TESTE) << std::endl;
					}
				}
				else
				{
					uint nodeid;
					for (int i=0; i<getNumVerticesForElementTypeMSH(type_aux); ++i)
					{
						File >> nodeid;
						--nodeid;
						
						if ((this->getNode(nodeid)->getLabel()) == 0)
						{
							this->getNode(nodeid)->setLabel(label_aux);
						}
					}
				}
				break;

			case 3: // caso volume
				uint nodeid;
				for (int i=0; i<getNumVerticesForElementTypeMSH(type_aux); ++i)
				{
					File >> nodeid;
					--nodeid;
					if (this->getNode(nodeid)->getLabel() == 0)
					{
						this->getNode(nodeid)->setLabel(label_aux);
					}
				}
				break;				
				
			default:
				std::cout << "invalid element ...\n";
		}
		
	}
	
    
	for (CellIterator cit = _cellL.begin(), cellend=_cellL.end() ; cit != cellend; ++cit)
	{
		if (!cit->isDead())
		{
			cit->propagateHalf(*this);
		}
	}
	
	for (MHalfIterator hit= _mhalfL.begin(), hend=_mhalfL.end(); hit != hend; ++hit)
	{
		if (!hit->isDead())
		{
			hit->propagateHalf(*this);
		}
	}
	
	
	File.close();
}


template<class Traits>
void iMesh<Traits>::writeFileState() {
	
	int order = this->_order;
    
	std::ofstream Fout( (_basename + ".state").data() );
	
	Fout << "BASENAME: " << _basename << std::endl << std::endl; 
	
	Fout << "NODES " << getNumNodesTotal() << std::endl;
	Fout << std::left << std::setw(22) << "x" << std::setw(22) << "y" << std::setw(22) << "z  label" << std::endl << std::endl;
	
	Fout.precision(15);
	
	for (uint i=0, tam=getNumNodesTotal(); i< tam; ++i)
	{
		this->getNode(i)->printSelfVTK(Fout); Fout << " " << this->getNode(i)->getLabel() << std::endl;
	}
	
	CellT *cell;
	Fout << std::endl;
	Fout << getElementNameMSH(CellT::getMSHTag(order)) << std::endl;
	for (uint i = 0; i < this->getNumCells(); i++)
	{
		cell = this->getCell(i);
		cell->printSelfState(Fout); Fout << " label:" << cell->getLabel();
		Fout << std::endl;
	}
	
	// print Halfs 
	Fout << std::endl << HalfT::getName() << std::endl;
	for (uint i = 0, tam=this->getNumCells(); i < tam; i++)
	{
		cell = this->getCell(i);
		
		HalfT *h_obj = cell->getHalf(0);
		
		for (int j = 1; j < CellT::N_borders; j++)
		{
			h_obj->printSelf(Fout);
			Fout << " | ";
			h_obj = cell->getHalf(j);
		}
		h_obj->printSelf(Fout);
		
		Fout << std::endl;
	}
	
	Fout << std::endl << "Marked " << HalfT::getName() << " " << this->getNumMHalfTotal() << std::endl;
	Fout << "CELL  POSIT  (ANC)  LABEL\n";
	for (uint i = 0, tam=this->getNumMHalfTotal(); i < tam; i++)
	{
		MHalfT *h_obj = this->getMHalf(i);
		h_obj->printSelf(Fout); Fout << " " << h_obj->getLabel();
		Fout << std::endl;
	}
	
	Fout << std::endl << std::endl;
	
	Fout.close();
}


/** Imprime os pontos e células no arquivo vtk.
 * @param flinear se true, força a malha ser impressa como linear, se false (default), a malha é impressa da ordem que ela é.
 * 
 * ESTA DESCRIÇÃO PRECISA SER MELHORADA
 */ 
template<class Traits>
void iMesh<Traits>::writeVTK(bool flinear)
{
	/* 
	 * 
	 * NOTA: TODOS os pontos são impressos. APENAS AS CÉLULAS VIVAS
	 * são impressas.
	 */ 

    int  order  = flinear ? 1 : this->_order;
    int  nc_mm  = CellT::getNumCellsMM(order);    // num de mini-células por mini-malha
    uint ncells = this->getNumCells() * nc_mm;
	std::stringstream ss;
	
    this->_add_scalar_vtk_n_calls=0;
    this->_add_vector_vtk_n_calls=0;
    
    if (this->_basename.size()==0)
    {
        _basename="untitled";
    }
    
	if (_is_family)
	{
		ss << (this->_basename) << std::setfill('0') << std::setw(5);
		ss << _family << ".vtk";
		++_family;
	}
	else
	{
		ss << (this->_basename) << ".vtk";
	}
	
	std::ofstream Fout( ss.str().data() );

	Fout << "# vtk DataFile Version 2.0" << std::endl
		 << "unstructured grid" << std::endl
		 << "ASCII" << std::endl
		 << "DATASET UNSTRUCTURED_GRID" << std::endl
		 << "POINTS " << _pointL.size() << " float" << std::endl;
	
	/* imprimindo pontos */
	PointIterator pit = _pointL.begin();
	for (auto pend=_pointL.end(); pit != pend ; ++pit)
	{
		pit->printSelfVTK(Fout);
		Fout << std::endl;
	}
	
	/* imprimindo células */
	Fout << std::endl;
	Fout << "CELLS " << ncells << " " << ncells*(1+CellT::N_vertices) << std::endl;
	CellIterator cit = _cellL.begin();
	for (auto cellend=_cellL.end(); cit != cellend ; ++cit)
	{
		if(!cit->isDead())
		{
			cit->printSelfVTK(Fout, order);
			Fout << std::endl;
		}
	}
    int type = CellT::getCellTypeVTK();
	Fout << std::endl; cit = _cellL.begin();
	Fout << "CELL_TYPES " << ncells << std::endl;
	for (auto cellend=_cellL.end(); cit != cellend ; ++cit)
	{
		if(!cit->isDead())
		{
            for (int i = 0; i < nc_mm; i++)
            {
                Fout << type << std::endl;
            }
		}
	}
	
	Fout << std::endl;

	
	Fout.close();
}

/** T: qualquer objeto chamável
 */ 
template<class Traits>
template<class T>
void iMesh<Traits>::addScalarVTK(const char* nome_var, T&& scalar, uint num_pts)
{
	
	std::stringstream ss;
	
	if (_is_family)
	{
		ss << (this->_basename) << std::setfill('0') << std::setw(5);
		ss << _family-1 << ".vtk";
	}
	else
	{
		ss << (this->_basename) << ".vtk";
	}
	
	std::ofstream Fout;
	Fout.open( ss.str().data(), std::ios::app);

	Fout << std::setprecision(8);
	Fout << std::setiosflags(std::ios::showpoint);
  
	if (_add_scalar_vtk_n_calls==0)
	{
		Fout << "POINT_DATA " << num_pts << std::endl;
		Fout << std::endl;
	}
	_add_scalar_vtk_n_calls++;
  	
  	Fout << "SCALARS " << nome_var << " float"   << std::endl;
	Fout << "LOOKUP_TABLE default"		        << std::endl;
  
	for (uint i=0; i<num_pts; ++i) {
		Fout << scalar[i] << std::endl;
	};
	
	Fout << std::endl << std::endl;
	
    Fout.close();
};


template<class Traits>
void iMesh<Traits>::addPointLabelVTK(const char* nome_var="node labels")
{
	
	std::stringstream ss;
	uint num_pts = this->getNumNodes();
    
	if (_is_family)
	{
		ss << (this->_basename) << std::setfill('0') << std::setw(5);
		ss << _family-1 << ".vtk";
	}
	else
	{
		ss << (this->_basename) << ".vtk";
	}
	
	std::ofstream Fout;
	Fout.open( ss.str().data(), std::ios::app);

	Fout << std::setprecision(8);
	Fout << std::setiosflags(std::ios::showpoint);
  
	if (_add_scalar_vtk_n_calls==0)
	{
		Fout << "POINT_DATA " << num_pts << std::endl;
		Fout << std::endl;
	}
	_add_scalar_vtk_n_calls++;
  	
  	Fout << "SCALARS " << nome_var << " int"   << std::endl;
	Fout << "LOOKUP_TABLE default"		        << std::endl;
  
	for (uint i=0; i<num_pts; ++i) {
		Fout << (getNode(i)->getLabel()) << std::endl;
	};
	
	Fout << std::endl << std::endl;
	
    Fout.close();
};


template<class Traits>
void iMesh<Traits>::addPointHalfVTK(const char* nome_var="node labels")
{
	
	std::stringstream ss;
	uint num_pts = this->getNumNodes();
    
	if (_is_family)
	{
		ss << (this->_basename) << std::setfill('0') << std::setw(5);
		ss << _family-1 << ".vtk";
	}
	else
	{
		ss << (this->_basename) << ".vtk";
	}
	
	std::ofstream Fout;
	Fout.open( ss.str().data(), std::ios::app);

	Fout << std::setprecision(8);
	Fout << std::setiosflags(std::ios::showpoint);
  
	if (_add_scalar_vtk_n_calls==0)
	{
		Fout << "POINT_DATA " << num_pts << std::endl;
		Fout << std::endl;
	}
	_add_scalar_vtk_n_calls++;
  	
  	Fout << "SCALARS " << nome_var << " int"   << std::endl;
	Fout << "LOOKUP_TABLE default"		        << std::endl;
  
	for (uint i=0; i<num_pts; ++i) {
		Fout << (this->getNode(i)->getHalf()->getIDCell()) << std::endl;
	};
	
	Fout << std::endl << std::endl;
	
    Fout.close();
};




#endif // IMESHIO_HPP