#include "mesh.hpp"


template<class Traits>
int iTriangle<Traits>::N_nodes = (MeshT::getOrder()+1)*(MeshT::getOrder()+2)/2;
