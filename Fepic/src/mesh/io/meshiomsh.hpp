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

#ifndef FEPIC_MESHIOMSH_HPP
#define FEPIC_MESHIOMSH_HPP

#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<typename _Traits::MeshT*>(this)
  #define CONST_THIS static_cast<const typename _Traits::MeshT*>(this)
#endif  


template<class _Traits>
class _MeshIoMsh
{


public:  
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::MeshT  MeshT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoMsh() {};
  _MeshIoMsh(_MeshIoMsh const&) {};
  
public:
  
  void readFileMsh(const char* filename);
  
  /// return 1 if error
  template<class Cell_T>
  bool cellTypeCheckMsh(int tag_type, typename std::enable_if<std::is_same<Simplex<1>,  typename Cell_T::PolytopeT>::value ||
                                                               std::is_same<Hypercube<1>,typename Cell_T::PolytopeT>::value ||
                                                               std::is_same<Polytope<1>, typename Cell_T::PolytopeT>::value >::type* = NULL)
  {
    const int tags[] = {MSH_LIN_2, MSH_LIN_3, MSH_LIN_4, MSH_LIN_5, MSH_LIN_6};
    const int* p = tags;
    do
      if (*p==tag_type) return false;
    while (++p != tags+std::extent<decltype(tags)>::value);
    return true;
  }

  /// return 1 if error
  template<class Cell_T>
  bool cellTypeCheckMsh(int tag_type, typename std::enable_if<std::is_same<Simplex<2>,  typename Cell_T::PolytopeT>::value>::type* = NULL)
  {
    const int tags[] = {MSH_TRI_3, MSH_TRI_6, MSH_TRI_9, MSH_TRI_10, MSH_TRI_12, MSH_TRI_15, MSH_TRI_21};
    const int* p = tags;
    do
      if (*p==tag_type) return false;
    while (++p != tags+std::extent<decltype(tags)>::value);
    return true;
  }
  
  /// return 1 if error
  template<class Cell_T>
  bool cellTypeCheckMsh(int tag_type, typename std::enable_if<std::is_same<Simplex<3>,  typename Cell_T::PolytopeT>::value>::type* = NULL)
  {
    const int tags[] = {MSH_TET_4, MSH_TET_10, MSH_TET_20, MSH_TET_35, MSH_TET_56, MSH_TET_34};
    const int* p = tags;
    do
      if (*p==tag_type) return false;
    while (++p != tags+std::extent<decltype(tags)>::value);
    return true;
  }
  
    
protected:


};



/* Lê a malha de um arquivo no formato .msh. A ordem da malha é detectada
 * altomaticamente.
 * TODO
 */
template<class _Traits>
void _MeshIoMsh<_Traits>::readFileMsh(const char* filename)
{

  /*
   * O que é feito:
   *
   * 1) primeiramente lê-se apenas as células (conectividade clássica);
   * 2) depois são construídos as half-edges e(ou) half-faces, lendo-se
   *    os elementos de dimensões menores.
   *
   * */

  
  THIS->_registerFile(filename, ".msh");

  std::ifstream File(filename);

  double  coord[3];
  int     type_aux;
  char    lixo[256];

  size_t    NODE_POS;  // $Nodes
  size_t    ELEM_POS;  // $Elements

  // nós
  NODE_POS = search_word(File, "$Nodes");
  NODE_POS++;  // only to avoid gcc warning: variable ‘NODE_POS’ set but not used [-Wunused-but-set-variable]

  uint num_pts;
  File >> num_pts;
  THIS->_pointL.resize(num_pts);

  for (uint i=0; i< THIS->getNumNodes(); ++i) {
    File >> lixo;  // nº do nó
    File >> coord[0];
    File >> coord[1];
    File >> coord[2];
    THIS->_pointL[i].setCoord(coord);
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
  bool wrong_file_err=true;
  for (uint k = 0; k < num_elms; ++k)
  {
    File >> lixo; // numero do elemento
    File >> type_aux;
    int elm_dim = getDimForElementTypeMsh(type_aux);
    if (elm_dim == CellT::dim)
    {
      if (!_MeshIoMsh::cellTypeCheckMsh<CellT>(type_aux))
      {
        THIS->setOrder( getOrderForElementTypeMsh(type_aux) );
        wrong_file_err=false;
      }
      break;
    }
    File.getline(lixo, 256);
  }
  FEPIC_ASSERT(!wrong_file_err, "wrong file format", std::invalid_argument);
  int order = THIS->_order;

  CellT   Cell;
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
    FEPIC_ASSERT(counter==k+1, "invalid file format", std::invalid_argument);

    File >> type_aux;
    int elm_dim = getDimForElementTypeMsh(type_aux);

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
        THIS->getNode(id_aux)->setTag(label_aux);
        break;

      case 1: // caso aresta
        {
          if (CellT::dim==1)
          {
            if (CellT::getMshTag(order) != type_aux)
            {
              std::cout << "Invalid edge! file has: " << getElementNameMsh(type_aux) << std::endl;
              std::cout << "But Mesh is configured to support: "
                                      << getElementNameMsh(CellT::getMshTag(order)) << std::endl;
              throw;
            }
            for (int i=0; i<Cell.getNumNodes(); ++i)
            {
              uint nodeid;
              File >> nodeid;
              --nodeid;
              Cell.setNode(i, nodeid);
            }
            Cell.setTag(label_aux);
            THIS->addCell(Cell);
            ++num_cells;
          }
          else
            File.getline(lixo, 256);
        }
        break;

      case 2: // caso face
        {
          if (CellT::dim==2)
          {
            if (CellT::getMshTag(order) != type_aux)
            {
              std::cout << "Invalid face! file has: " << getElementNameMsh(type_aux) << std::endl;
              std::cout << "But Mesh is configured to support: "
                                      << getElementNameMsh(CellT::getMshTag(order)) << std::endl;
              throw;
            }
            for (int i=0; i<Cell.getNumNodes(); ++i)
            {
              uint nodeid;
              File >> nodeid;
              --nodeid;
              Cell.setNode(i, nodeid);
            }
            Cell.setTag(label_aux);
            THIS->addCell(Cell);
            ++num_cells;
          }
          else
            File.getline(lixo, 256);
        }
        break;

      case 3: // caso volume
        {
          if (CellT::dim==3)
          {
            if (CellT::getMshTag(order) != type_aux)
            {
              std::cout << "Invalid volume! file has: " << getElementNameMsh(type_aux) << std::endl;
              std::cout << "But Mesh is configured to support: "
                                      << getElementNameMsh(CellT::getMshTag(order)) << std::endl;
              throw;
            }
            for (int i=0; i<Cell.getNumNodes(); ++i)
            {
              uint nodeid;
              File >> nodeid;
              --nodeid;
              Cell.setNode(i, nodeid);
            }
            Cell.setTag(label_aux);
            THIS->addCell(Cell);
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
              << getElementNameMsh(CellT::getMshTag(order)) << std::endl;

  /* até aqui, apenas foi lido a conectividade
   */

  /* constroi as halfs e Mhalfs */
  _MeshMethods<MeshT, CellT::dim>::buildAdjacency(*THIS);


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
    int elm_dim = getDimForElementTypeMsh(type_aux);

    File >> num_tags;

    File >> label_aux;

    for (int j=0; j<num_tags-1; ++j)
    {
      File >> lixo;
    }

    uint id_aux;
    switch (elm_dim) {
      case 0: // caso ponto
        File >> id_aux;
        break;

      case 1: // caso aresta
        if (CellT::dim-1==1)
        {
          vectorui nodes(numNodes<CellDerivedPolytopeT>(order));
          vectorui vtx(2);
          uint half_id;
          for (int i=0; i<nodes.size(); ++i)
          {
            File >> nodes[i];
            --nodes[i];
            if (THIS->getNode(nodes[i])->getTag() == 0)
              THIS->getNode(nodes[i])->setTag(label_aux);
          }
          copy( nodes.begin(), nodes.begin()+2, vtx.begin() );
          if( THIS->theseVerticesFormAHalfl(vtx, half_id) )
            THIS->getHalfl(half_id)->setTag(label_aux);
        }
        else
        {
          uint nodeid;
          for (int i=0; i<getNumVerticesForElementTypeMsh(type_aux); ++i)
          {
            File >> nodeid;
            --nodeid;
            if (THIS->getNode(nodeid)->getTag() == 0)
              THIS->getNode(nodeid)->setTag(label_aux);
          }
        }
        break;

      case 2: // caso face
        if (CellT::dim - 1 == 2)
        {
          vectorui nodes(numNodes<CellDerivedPolytopeT>(order));
          vectorui vtx(CellT::n_vertices_per_border);
          uint half_id=0;
          for (int i=0; i<nodes.size(); ++i)
          {
            File >> nodes[i];
            --nodes[i];
            if (THIS->getNode(nodes[i])->getTag() == 0)
              THIS->getNode(nodes[i])->setTag(label_aux);
          }
          copy( nodes.begin(), nodes.begin()+CellT::n_vertices_per_border, vtx.begin() );
          if( THIS->theseVerticesFormAHalfl(vtx, half_id) )
            THIS->getHalfl(half_id)->setTag(label_aux); //std::cout << (++TESTE) << std::endl;
        }
        else
        {
          uint nodeid;
          for (int i=0; i<getNumVerticesForElementTypeMsh(type_aux); ++i)
          {
            File >> nodeid;
            --nodeid;

            if ((THIS->getNode(nodeid)->getTag()) == 0)
              THIS->getNode(nodeid)->setTag(label_aux);
          }
        }
        break;

      case 3: // caso volume
        uint nodeid;
        for (int i=0; i<getNumVerticesForElementTypeMsh(type_aux); ++i)
        {
          File >> nodeid;
      
          --nodeid;
          if (THIS->getNode(nodeid)->getTag() == 0)
            THIS->getNode(nodeid)->setTag(label_aux);
        }
        break;

      default:
        std::cout << "invalid element ...\n";
    }

  } // end for


  for (auto cit = THIS->_cellL.begin(), cellend=THIS->_cellL.end() ; cit != cellend; ++cit)
    if (!cit->disabled())
      cit->broadcastHalf2Nodes(*THIS);

  for (auto hit= THIS->_mhalfL.begin(), hend=THIS->_mhalfL.end(); hit != hend; ++hit)
    if (!hit->disabled())
      hit->broadcastHalf2Nodes(*THIS);


  File.close();
}

#undef THIS
#undef CONST_THIS

#endif
