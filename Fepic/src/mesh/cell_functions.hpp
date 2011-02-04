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

#ifndef FEPIC_CELL_FUNCTIONS_HPP
#define FEPIC_CELL_FUNCTIONS_HPP

// -------------------------------------------------------------------------

template<> inline uint numNodes<Polytope<1>>(uint order) {
  return order+1; }

template<> inline uint numNodes<Simplex<1>>(uint order) {
  return order+1; }

template<> inline uint numNodes<Simplex<2>>(uint order) {
  return (order+1)*(order+2)/2; }
  
template<> inline uint numNodes<Simplex<3>>(uint order) {
  return (order+1)*(order+2)*(order+3)/6; }

template<> inline uint numNodes<Hypercube<1>>(uint order) {
  return order+1; }

template<> inline uint numNodes<Hypercube<2>>(uint order) {
  return (order+1)*(order+1); }

template<> inline uint numNodes<Hypercube<3>>(uint order) {
  return (order+1)*(order+1)*(order+1); }

// -------------------------------------------------------------------------

template<> inline uint numBubbles<Polytope<1>>(uint order) {
  return (order-1); }

template<> inline uint numBubbles<Simplex<1>>(uint order) {
  return (order-1); }

template<> inline uint numBubbles<Simplex<2>>(uint order) {
  return (order-1)*(order-2)/2; }

template<> inline uint numBubbles<Simplex<3>>(uint order) {
  return (order-1)*(order-2)*(order-3)/6;; }

template<> inline uint numBubbles<Hypercube<1>>(uint order) {
  return order-1; }

template<> inline uint numBubbles<Hypercube<2>>(uint order) {
  return (order-1)*(order-1); }

template<> inline uint numBubbles<Hypercube<3>>(uint order) {
  return (order-1)*(order-1)*(order-1); }




#endif
