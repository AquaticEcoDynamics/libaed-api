!###############################################################################
!#                                                                             #
!# aed_api_env.h                                                               #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2024 -  The University of Western Australia                      #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#                                                                             #
!# NB - This is currently a "work-in-progress" do not rely on it's function    #
!#      declarations/constants/whatever because it will almost certainly be    #
!#      changed.                                                               #
!#                                                                             #
!###############################################################################
#ifndef _AED_API_ENV_H_
#define _AED_API_ENV_H_

/* Environment Variables provided by driver model */

/* 3D vars */
      /* 'temperature'    'temperature'            'celsius'        */
#define A_E_TEMP       1
      /* 'salinity'       'salinity'               'g/kg'           */
#define A_E_SALT       2
      /* 'density'        'density'                'kg/m3'          */
#define A_E_DENS       3
      /* 'layer_ht'       'layer heights'          'm'              */
#define A_E_LAYER_HT   4
      /* 'extc_coef'      'extinction coefficient' '/m'             */
#define A_E_EXTC_COEF  5
      /* 'tss'            'tss'                    'g/m3'           */
#define A_E_TSS        6
      /* 'ss1'            'ss1'                    'g/m3'           */
#define A_E_SS1        7
      /* 'ss2'            'ss2'                    'g/m3'           */
#define A_E_SS2        8
      /* 'ss3'            'ss3'                    'g/m3'           */
#define A_E_SS3        9
      /* 'ss4'            'ss4'                    'g/m3'           */
#define A_E_SS4        10
      /* 'cell_vel'       'cell velocity'          'm/s'            */
#define A_E_CELL_VEL   11
      /* 'nir'            'nir'                    'W/m2'           */
#define A_E_NIR        12
      /* 'par'            'par'                    'W/m2'           */
#define A_E_PAR        13
      /* 'uva'            'uva'                    'W/m2'           */
#define A_E_UVA        14
      /* 'uvb'            'uvb'                    'W/m2'           */
#define A_E_UVB        15
      /* 'pressure'       'pressure'               ''               */
#define A_E_PRESSURE   16
      /* 'depth'          'depth'                  'm'              */
#define A_E_DEPTH      17

/* 2D vars */
      /* 'layer_area'     'layer area'             'm2'             */
#define A_E_LAYER_AREA 101
      /* 'rain'           'rainfall'               'm/s'            */
#define A_E_RAIN       102
      /* 'rainloss'       'rain loss'              'm/s'            */
#define A_E_RAINLOSS   103
      /* 'material'       'material'               '-'              */
#define A_E_MATERIAL   104
      /* 'bathy'          'bathy'                  'm above datum'  */
#define A_E_BATHY      105
      /* 'sed_zone'       'sediment zone'          '-'              */
#define A_E_SED_ZONE   106
      /* 'wind_speed'     'wind speed'             'm/s'            */
#define A_E_WIND_SPEED 107
      /* 'par_sf'         'par_sf'                 'W/m2'           */
#define A_E_PAR_SF     108
      /* 'taub'           'layer stress'           'N/m2'           */
#define A_E_TAUB       109
      /* 'air_temp'       'air temperature'        'celsius'        */
#define A_E_AIR_TEMP   110
      /* 'air_pres'       'air pressure'           'Pa'             */
#define A_E_AIR_PRES   111
      /* 'humidity'       'relative humidity'      '-'              */
#define A_E_HUMIDITY   112
      /* 'longwave'       'longwave'               'W/m2'           */
#define A_E_LONGWAVE   113
      /* 'col_num'        'column number'          '-'              */
#define A_E_COL_NUM    114
      /* 'col_depth'      'column water depth'     'm above bottom' */
#define A_E_COL_DEPTH  115
      /* 'longitude'      'longitude'              'radians'        */
#define A_E_LONGITUDE  116
      /* 'latitude'       'latitude'               'radians'        */
#define A_E_LATITUDE   117
      /* 'nearest_active' 'nearest active'         '-'              */
!#define A_E_NEAR_ACTV  118
      /* 'nearest_depth'  'nearest depth'          'm'              */
!#define A_E_NEAR_DEPTH 119

/* 0D vars */
      /* 'yearday'        'yearday'                'day'            */
#define A_E_YEARDAY    201
      /* 'timestep'       'timestep'               'seconds'        */
#define A_E_TIMESTEP   202


#endif
