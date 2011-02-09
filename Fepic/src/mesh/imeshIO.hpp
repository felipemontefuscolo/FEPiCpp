// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.

#ifndef FEPIC_IMESHIO_HPP
#define FEPIC_IMESHIO_HPP





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
template<class _Traits>
void iMesh<_Traits>::buildAdjacency4face()
{

  int n_borders = CellT::n_borders;
  vectorui edge_vtx(2);
  vectorui cell_ith(2);
  std::map<vectorui, vectorui> table; // Key: (n0, n1) ... atributos (elemento, ith)
  std::map<vectorui, vectorui>::iterator tab_it, tab_it2;
  CellIterator cell = iMesh::_cellL.begin();
  uint otherC, otherith, thisC, thisith;
  int thistag;



  uint k=0;
  for (auto cellend=_cellL.end(); cell != cellend ; ++cell)
  {
    if (!cell->disabled())
    {
      for (int ith = 0; ith < n_borders; ++ith)
      {
        edge_vtx = cell->getBorderVertices(ith);
        cell_ith[0] = k;
        cell_ith[1] = ith;
        table.insert(std::pair<vectorui, vectorui>(edge_vtx, cell_ith));
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

      this->getCell(thisC)->getHalf(thisith)->setCompleteId(otherC, otherith);
      this->getCell(otherC)->getHalf(otherith)->setCompleteId(thisC, thisith);

      table.erase(tab_it2);
      table.erase(tab_it);
    }
    else
    {
      thisC = tab_it->second[0];
      thisith =  tab_it->second[1];
      thistag = this->getCell(thisC)->getTag();
      HalfLT temp(thisC, thisith, 0, thistag);
      uint mHE_universal_iD = addHalfl(temp);
      this->getCell(thisC)->getHalf(thisith)->setCompleteId(mHE_universal_iD, -1);
      table.erase(tab_it);
    }
  }

}


/** @note células mortas não são consideradas.
 */
template<class _Traits>
void iMesh<_Traits>::buildAdjacency4volume()
{

  const int n_borders = CellT::n_borders;
  const int n_anch    = CellT::n_vertices_per_border;

  typedef std::pair<vectorui, vectorui>   PairT; // < (vértices da face) , (elemento, ith) >
  typedef std::map<vectorui, vectorui>  MapT;

  vectorui face_vtx(CellT::n_vertices_per_border);
  vectorui cell_ith(2);
  MapT            table; // < (nós da face) , (elemento, ith) >
  MapT::iterator  tab_it, tab_it2;
  CellIterator    cell = this->_cellL.begin();
  uint            otherC, otherith, thisC, thisith;
  int             thistag;
  int             a;

  uint k=0;
  for (auto cellend=_cellL.end(); cell != cellend ; ++cell)
  {
    if (!cell->disabled())
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

        _cellL[thisC].getHalf(thisith)->setCompleteId(otherC, otherith, a);
        _cellL[otherC].getHalf(otherith)->setCompleteId(thisC, thisith, a);

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
      thistag = this->getCell(thisC)->getTag();
      HalfLT temp(thisC, thisith, 0, thistag);
      uint mHE_universal_iD = addHalfl(temp);
      this->getCell(thisC)->getHalf(thisith)->setCompleteId(mHE_universal_iD, -1, 0);
      table.erase(tab_it);
    }

  }


}

//template<class _Traits>
//void iMesh<_Traits>::writeFileState() {

  //int order = this->_order;

  //std::ofstream Fout( (_basename + ".state").data() );

  //Fout << "BASENAME: " << _basename << std::endl << std::endl;

  //Fout << "NODES " << getNumNodesTotal() << std::endl;
  //Fout << std::left << std::setw(22) << "x" << std::setw(22) << "y" << std::setw(22) << "z  label" << std::endl << std::endl;

  //Fout.precision(15);

  //for (uint k=0, tam=getNumNodesTotal(); k< tam; ++k)
  //{
    //this->getNode(k)->printSelfVtk(Fout);
    //Fout << " " << this->getNode(k)->getTag() << std::endl;
	//}

  //CellT *cell;
  //Fout << std::endl;
  //Fout << getElementNameMsh(CellT::getMshTag(order)) << std::endl;
  //for (uint i = 0; i < this->getNumCells(); i++)
  //{
    //cell = this->getCell(i);
    //cell->printSelfState(Fout); Fout << " label:" << cell->getTag();
    //Fout << std::endl;
  //}

  //// print Halfs
  //Fout << std::endl << HalfT::getName() << std::endl;
  //for (uint i = 0, tam=this->getNumCells(); i < tam; i++)
  //{
    //cell = this->getCell(i);

    //HalfT *h_obj = cell->getHalf(0);

    //for (int j = 1; j < CellT::n_borders; j++)
    //{
      //h_obj->printSelf(Fout);
      //Fout << " | ";
      //h_obj = cell->getHalf(j);
    //}
    //h_obj->printSelf(Fout);

    //Fout << std::endl;
  //}

  //Fout << std::endl << "Marked " << HalfT::getName() << " " << this->getNumHalfLTotal() << std::endl;
  //Fout << "CELL  POSIT  (ANC)  LABEL\n";
  //for (uint i = 0, tam=this->getNumHalfLTotal(); i < tam; i++)
  //{
    //HalfLT *h_obj = this->getHalfl(i);
    //h_obj->printSelf(Fout); Fout << " " << h_obj->getTag();
    //Fout << std::endl;
  //}

  //Fout << std::endl << std::endl;

  //Fout.close();
//}





#endif // IMESHIO_HPP
