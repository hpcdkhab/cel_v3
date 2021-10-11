/*     
 * File:   cel_domain_3d.h
 * Project CRESTA (see details on https://cresta-project.eu) 
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
#ifndef CEL_DOMAIN_3D_H
#define CEL_DOMAIN_3D_H
#include "cel_types.h"

//! \struct 3D description of regular grid on domain [0;1]x[0;1]x[0;1]
typedef struct 
{
	INDEX block_indx; //! block_indx index of sub-domain
	INDEX nx;//! nx number of nodes in X-direction
	INDEX ny;//! ny number of nodes in Y-direction
	INDEX nz;//! nz number of nodes in Z-direction
	INDEX dx;//! dx number of sub-domains in X-direction
	INDEX dy;//! dy number of sub-domains in Y-direction
	INDEX dz;//! dz number of sub-domains in Z-direction
	INDEX blocks;//! blocks number blocks (sub-domains) in 3D domain
	INDEX bx;//! bx number nodes in X-direction in sub-domain
	INDEX by;//! by number nodes in Y-direction in sub-domain
	INDEX bz;//! bz number nodes in Z-direction in sub-domain
	INDEX bxyz;//! number nodes in sub-domain
} cel_domain_3D;



//! fill struct domain_3D with parameter for regular grid on domain  [0;1]x[0;1]x[0;1]
//! @param (in) nx number of nodes in X-direction
//! @param (in) ny number of nodes in Y-direction
//! @param (in) nz number of nodes in Z-direction
//! @param (in) dx number of sub-domains in X-direction
//! @param (in) dy number of sub-domains in Y-direction
//! @param (in) dz number of sub-domains in Z-direction
//! @param (out) domain_3D pointer to struct cel_domain_3D - description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) output_on > 1 - print error to sdterr; 0 - no output
//! @return  0 - no error ; -1 - wrong decomposition's parameters
int cel_domain_3D_init(INDEX nx, INDEX ny, INDEX nz, INDEX dx, INDEX dy, INDEX dz, cel_domain_3D* domain_3D, int output_on);

//! print parameter of regular grid on domain  [0;1]x[0;1]x[0;1]
//! @param (in) domain_3D pointer to struct cel_domain_3D - description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) output_on > 1 - print to sdterr; 0 - no output
void cel_domain_3D_print(cel_domain_3D* domain_3D, int output_on);

//! convert global coordinate indices of node to global and local indices of node in regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) domain_3D pointer to struct cel_domain_3D - description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) node_x_indx coordinate index of node in X-direction
//! @param (in) node_y_indx coordinate index of node in Y-direction
//! @param (in) node_z_indx coordinate index of node in Z-direction
//! @param (out) pointer to block index;
//! @param (out) pointer to local node index 
//! @return  global node index or INDEX_MAX if global coordinate indices of node are out of range (not in the domain)
inline INDEX cel_domain_3D_get_index(cel_domain_3D* domain_3D, INDEX node_x_indx, INDEX node_y_indx, INDEX node_z_indx, INDEX* block_idx, INDEX* cell_indx_local);

//! convert global node index to global coordinate indices in regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) domain_3D pointer to struct cel_domain_3D - description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in)  cell_indx index of node
//! @param (out) node_x_indx pointer to coordinate index of node in X-direction
//! @param (out) node_y_indx pointer to coordinate index of node in Y-direction
//! @param (out) node_z_indx pointer to coordinate index of node in Z-direction
inline void cel_domain_3D_get_node_coordinate_indicies(cel_domain_3D* domain_3D, INDEX cell_indx, INDEX* node_x_indx, INDEX* node_y_indx, INDEX* node_z_indx);

//! get boundaries coordinate indices of sub block in regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in) domain_3D pointer to struct cel_domain_3D - description of regular grid on domain [0;1]x[0;1]x[0;1]
//! @param (in)  block_indx index of sub block
//! @param (out) xx_start pointer to boundaries coordinate indices
//! @param (out) yy_start pointer to boundaries coordinate indices
//! @param (out) zz_start pointer to boundaries coordinate indices
//! @param (out) xx_end pointer to boundaries coordinate indices
//! @param (out) yy_end pointer to boundaries coordinate indices
//! @param (out) zz_start pointer to boundaries coordinate ndices
inline void cel_domain_3D_get_block_boundary(cel_domain_3D* domain_3D, INDEX block_indx, INDEX* xx_start, INDEX* yy_start, INDEX* zz_start, INDEX* xx_end, INDEX* yy_end, INDEX* zz_end);

#endif

