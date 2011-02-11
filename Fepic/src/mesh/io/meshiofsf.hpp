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

#ifndef FEPIC_MESHIOFSF_HPP
#define FEPIC_MESHIOFSF_HPP


#if !defined(THIS) && !defined(CONST_THIS)
  #define THIS static_cast<typename _Traits::MeshT*>(this)
  #define CONST_THIS static_cast<const typename _Traits::MeshT*>(this)
#endif 


template<class _Traits>
class _MeshIoFsf
{
  
public:  
  typedef typename _Traits::CellT  CellT;
  typedef typename _Traits::MeshT  MeshT;
  typedef typename _Traits::PointT  PointT;
  typedef typename CellT::PolytopeT CellPolytopeT;
  typedef typename CellPolytopeT::Derived CellDerivedPolytopeT;

protected:
  _MeshIoFsf() : _filenumFsf(0), _add_scalar_fsf_n_calls(0) {};
  _MeshIoFsf(_MeshIoFsf const&) {};

public:

  void readFileFsf(const char* filename);
  void writeFsf(bool flinear=false);
  template<class T>
  void addScalarFsf(const char* nome_var, T&& scalar, uint num_pts);
  template<class T>
  void addVectorFsf(const char* nome_var, T&& arrayos, int dim, uint num_pts);
  void addPointTagFsf(const char* nome_var); // para debug
  void addPointHalfFsf(const char* nome_var);  // para debug     

  
protected:  
  uint _filenumFsf;
  uint _add_scalar_fsf_n_calls;
};



#undef THIS
#undef CONST_THIS

#endif


