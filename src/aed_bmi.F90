!###############################################################################
!#                                                                             #
!# aed_bmi.F90                                                                 #
!#                                                                             #
!# A generic interface between model and libaed-xxx                            #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2024 - 2025 - The University of Western Australia                 #
!#                                                                             #
!#  This file is part of libaed (Library for AquaticEco Dynamics)              #
!#                                                                             #
!#  AED is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  AED is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!#  This module is not yet implemented ... DO NO USE!                          #
!#                                                                             #
!###############################################################################

#include "aed_api.h"

! The Basic Model Interface (BMI) Fortran specification.
!
! This language specification is derived from the Scientific
! Interface Definition Language (SIDL) file bmi.sidl located at
! https://github.com/csdms/bmi.

MODULE aed_bmif

  IMPLICIT NONE

  INTEGER,PARAMETER :: BMI_MAX_COMPONENT_NAME = 2048
  INTEGER,PARAMETER :: BMI_MAX_VAR_NAME = 2048
  INTEGER,PARAMETER :: BMI_MAX_TYPE_NAME = 2048
  INTEGER,PARAMETER :: BMI_MAX_UNITS_NAME = 2048

  INTEGER,PARAMETER :: BMI_FAILURE = 1
  INTEGER,PARAMETER :: BMI_SUCCESS = 0

CONTAINS

!###############################################################################
FUNCTION bmif_initialize(this, config_file) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Perform startup tasks for the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(out) :: this
   CHARACTER(len=*), INTENT(in) :: config_file
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_initialize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_update(this) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Advance the model one time step.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_update_until(this, time) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Advance the model until the given time.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
   DOUBLE PRECISION, INTENT(in) :: time
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_update_until
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_finalize(this) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Perform teardown tasks for the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_finalize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_component_name(this, name) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the name of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), pointer, INTENT(out) :: name
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_component_name
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_input_item_count(this, count) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Count a model's input variables.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(out) :: count
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_input_item_count
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_output_item_count(this, count) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Count a model's output variables.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(out) :: count
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_output_item_count
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_input_var_names(this, names) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! List a model's input variables.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), pointer, INTENT(out) :: names(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_input_var_names
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_output_var_names(this, names) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! List a model's output variables.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), pointer, INTENT(out) :: names(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_output_var_names
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_grid(this, name, grid) RESULT(bmi_status)
!-------------------------------------------------------------------------------
    ! Get the grid identifier for the given variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: grid
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_grid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_type(this, name, type) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the data type of the given variable as a string.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   CHARACTER(len=*), INTENT(out) :: type
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_units(this, name, units) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the units of the given variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   CHARACTER(len=*), INTENT(out) :: units
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_units
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_itemsize(this, name, size) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get memory use per array element, in bytes.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: size
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_itemsize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_nbytes(this, name, nbytes) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get size of the given variable, in bytes.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: nbytes
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_nbytes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_var_location(this, name, location) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Describe where a variable is located: node, edge, or face.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   CHARACTER(len=*), INTENT(out) :: location
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_var_location
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_current_time(this, time) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Current time of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   DOUBLE PRECISION, INTENT(out) :: time
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_current_time
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_start_time(this, time) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Start time of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   DOUBLE PRECISION, INTENT(out) :: time
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_start_time
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_end_time(this, time) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! End time of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   DOUBLE PRECISION, INTENT(out) :: time
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_end_time
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_time_units(this, units) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Time units of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(out) :: units
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_time_units
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_time_step(this, time_step) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Time step of the model.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   DOUBLE PRECISION, INTENT(out) :: time_step
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_time_step
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_int(this, name, dest) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a copy of values (flattened!) of the given INTEGER variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(inout) :: dest(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_int
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_float(this, name, dest) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a copy of values (flattened!) of the given real variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   real, INTENT(inout) :: dest(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_float
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_double(this, name, dest) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a copy of values (flattened!) of the given double variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   DOUBLE PRECISION, INTENT(inout) :: dest(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_ptr_int(this, name, dest_ptr) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a reference to the given INTEGER variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, pointer, INTENT(inout) :: dest_ptr(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_ptr_int
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_ptr_float(this, name, dest_ptr) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a reference to the given real variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   real, pointer, INTENT(inout) :: dest_ptr(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_ptr_float
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_ptr_double(this, name, dest_ptr) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get a reference to the given double variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   DOUBLE PRECISION, pointer, INTENT(inout) :: dest_ptr(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_ptr_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_at_indices_int(this, name, dest, inds) &
!-------------------------------------------------------------------------------
! Get INTEGER values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(inout) :: dest(:)
   INTEGER, INTENT(in) :: inds(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_at_indices_int
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_at_indices_float(this, name, dest, inds) &
!-------------------------------------------------------------------------------
! Get real values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   real, INTENT(inout) :: dest(:)
   INTEGER, INTENT(in) :: inds(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_at_indices_float
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_value_at_indices_double(this, name, dest, inds) &
!-------------------------------------------------------------------------------
! Get double values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   CHARACTER(len=*), INTENT(in) :: name
   DOUBLE PRECISION, INTENT(inout) :: dest(:)
   INTEGER, INTENT(in) :: inds(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_value_at_indices_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_int(this, name, src) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Set new values for an INTEGER model variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_int
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_float(this, name, src) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Set new values for a real model variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   real, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_float
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_double(this, name, src) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Set new values for a double model variable.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   DOUBLE PRECISION, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_at_indices_int(this, name, inds, src) &
!-------------------------------------------------------------------------------
! Set INTEGER values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(in) :: inds(:)
   INTEGER, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_at_indices_int
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_at_indices_float(this, name, inds, src) &
!-------------------------------------------------------------------------------
! Set real values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(in) :: inds(:)
   real, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_at_indices_float
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_set_value_at_indices_double(this, name, inds, src) &
!-------------------------------------------------------------------------------
! Set double values at particular (one-DIMENSIONal) indices.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(inout) :: this
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(in) :: inds(:)
   DOUBLE PRECISION, INTENT(in) :: src(:)
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_set_value_at_indices_double
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_rank(this, grid, rank) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get number of DIMENSIONs of the computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, INTENT(out) :: rank
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_rank
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_size(this, grid, size) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the total number of elements in the computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, INTENT(out) :: size
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_size
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_type(this, grid, type) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the grid type as a string.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   CHARACTER(len=*), INTENT(out) :: type
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_shape(this, grid, shape) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the DIMENSIONs of the computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, DIMENSION(:), INTENT(out) :: shape
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_shape
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_spacing(this, grid, spacing) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get distance between nodes of the computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: spacing
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_spacing
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_origin(this, grid, origin) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get coordinates of the origin of the computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: origin
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_origin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_x(this, grid, x) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the x-coordinates of the nodes of a computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: x
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_x
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_y(this, grid, y) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the y-coordinates of the nodes of a computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: y
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_y
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_z(this, grid, z) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the z-coordinates of the nodes of a computational grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: z
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_z
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_node_count(this, grid, count) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the number of nodes in an unstructured grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, INTENT(out) :: count
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_node_count
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_edge_count(this, grid, count) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the number of edges in an unstructured grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, INTENT(out) :: count
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_edge_count
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_face_count(this, grid, count) RESULT(bmi_status)
!-------------------------------------------------------------------------------
! Get the number of faces in an unstructured grid.
!-------------------------------------------------------------------------------
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, INTENT(out) :: count
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_face_count
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_edge_nodes(this, grid, edge_nodes) &
!-------------------------------------------------------------------------------
! Get the edge-node connectivity.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, DIMENSION(:), INTENT(out) :: edge_nodes
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_edge_nodes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_face_edges(this, grid, face_edges) &
!-------------------------------------------------------------------------------
! Get the face-edge connectivity.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, DIMENSION(:), INTENT(out) :: face_edges
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_face_edges
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_face_nodes(this, grid, face_nodes) &
!-------------------------------------------------------------------------------
! Get the face-node connectivity.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, DIMENSION(:), INTENT(out) :: face_nodes
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_face_nodes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION bmif_get_grid_nodes_per_face(this, grid, nodes_per_face) &
!-------------------------------------------------------------------------------
! Get the number of nodes for each face.
!-------------------------------------------------------------------------------
   RESULT(bmi_status)
   CLASS(bmi), INTENT(in) :: this
   INTEGER, INTENT(in) :: grid
   INTEGER, DIMENSION(:), INTENT(out) :: nodes_per_face
!
!LOCALS
   INTEGER :: bmi_status = BMI_FAILURE
!
!-------------------------------------------------------------------------------
!BEGIN
END FUNCTION bmif_get_grid_nodes_per_face
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_bmif
