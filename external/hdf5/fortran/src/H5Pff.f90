! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
!   access to either file, you may request a copy from help@hdfgroup.org.     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
!
! This file contains Fortran90 interfaces for H5P functions.
!
     MODULE H5P

       USE H5GLOBAL

       INTERFACE h5pset_fill_value_f
         MODULE PROCEDURE h5pset_fill_value_integer
         MODULE PROCEDURE h5pset_fill_value_real
! Comment if on Crays
         MODULE PROCEDURE h5pset_fill_value_double
! End comment if on Crays
         MODULE PROCEDURE h5pset_fill_value_char
       END INTERFACE
     
       INTERFACE h5pget_fill_value_f
         MODULE PROCEDURE h5pget_fill_value_integer
         MODULE PROCEDURE h5pget_fill_value_real
! Comment if on Crays
         MODULE PROCEDURE h5pget_fill_value_double
! End comment if on Crays
         MODULE PROCEDURE h5pget_fill_value_char
       END INTERFACE

       INTERFACE h5pset_f
         MODULE PROCEDURE h5pset_integer
         MODULE PROCEDURE h5pset_real
! Comment if on Crays
         MODULE PROCEDURE h5pset_double
! End comment if on Crays
         MODULE PROCEDURE h5pset_char
       END INTERFACE
     
       INTERFACE h5pget_f
         MODULE PROCEDURE h5pget_integer
         MODULE PROCEDURE h5pget_real
! Comment if on Crays
         MODULE PROCEDURE h5pget_double
! End comment if on Crays
         MODULE PROCEDURE h5pget_char
       END INTERFACE

       INTERFACE h5pregister_f
         MODULE PROCEDURE h5pregister_integer
         MODULE PROCEDURE h5pregister_real
! Comment if on Crays
         MODULE PROCEDURE h5pregister_double
! End comment if on Crays
         MODULE PROCEDURE h5pregister_char
       END INTERFACE
     
       INTERFACE h5pinsert_f
         MODULE PROCEDURE h5pinsert_integer
         MODULE PROCEDURE h5pinsert_real
! Comment if on Crays
         MODULE PROCEDURE h5pinsert_double
! End comment if on Crays
         MODULE PROCEDURE h5pinsert_char
       END INTERFACE
       
       INTERFACE h5pset_fapl_multi_f
         MODULE PROCEDURE h5pset_fapl_multi_l
         MODULE PROCEDURE h5pset_fapl_multi_s
       END INTERFACE        

     
     CONTAINS

!----------------------------------------------------------------------
! Name:		h5pcreate_f 
!
! Purpose: 	Creates a new property as an instance of a property 
!		list class.
!
! Inputs:  
!		class		- type of the property class to be created.
!				  Possible values are:
!				  H5P_FILE_CREATE_F
!				  H5P_FILE_ACCESS_F
!				  H5P_DATASET_CREATE_F
!				  H5P_DATASET_XFER_F
!				  H5P_FILE_MOUNT_F
! Outputs:  
!		prp_id		- property list identifier
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------
          SUBROUTINE h5pcreate_f(class, prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pcreate_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class      ! The type of the property list 
                                              ! to be created. Possible values
                                              ! are: 
                                              !  H5P_FILE_CREATE_F
                                              !  H5P_FILE_ACCESS_F
                                              !  H5P_DATASET_CREATE_F
                                              !  H5P_DATASET_XFER_F
                                              !  H5P_FILE_MOUNT_F
            INTEGER(HID_T), INTENT(OUT) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr        ! Error code

!            INTEGER, EXTERNAL :: h5pcreate_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pcreate_c(class, prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCREATE_C'::h5pcreate_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: class
              INTEGER(HID_T), INTENT(OUT) :: prp_id
              END FUNCTION h5pcreate_c
            END INTERFACE

            hdferr = h5pcreate_c(class, prp_id) 
          END SUBROUTINE h5pcreate_f

!----------------------------------------------------------------------
! Name:		h5pset_preserve_f 
!
! Purpose: 	Sets the dataset transfer property list status to 
!		TRUE or FALSE for initializing compound datatype
!		members during write/read operations.
!
! Inputs:  
!		prp_id		- property list identifier
!		flag		- status flag
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!                       Datatype of the flag parameter is changed from 
!                       INTEGER to LOGICAL 
!                               June 4, 2003 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_preserve_f(prp_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_preserve_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            LOGICAL, INTENT(IN) ::  flag ! TRUE/FALSE flag to set the dataset
                                         ! transfer property for partila writing/reading
                                         ! compound datatype 
            INTEGER, INTENT(OUT) :: hdferr    ! Error code
            INTEGER :: flag_c

!            INTEGER, EXTERNAL :: h5pset_preserve_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_preserve_c(prp_id, flag_c)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_PRESERVE_C'::h5pset_preserve_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER ::  flag_c
              END FUNCTION h5pset_preserve_c
            END INTERFACE
            flag_c = 0
            if(flag) flag_c = 1
            hdferr = h5pset_preserve_c(prp_id, flag_c) 
          END SUBROUTINE h5pset_preserve_f

!----------------------------------------------------------------------
! Name:		h5pget_preserve_f 
!
! Purpose: 	Checks status of the dataset transfer property list.
!
! Inputs:  
!		prp_id		- property list identifier
! Outputs:  
!		flag		- status flag
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!                       Datatype of the flag parameter is changed from 
!                       INTEGER to LOGICAL 
!                               June 4, 2003 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_preserve_f(prp_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_preserve_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            LOGICAL, INTENT(OUT) ::  flag ! TRUE/FALSE flag. Shows status of the dataset's
                                         ! transfer property for partial writing/reading
                                         ! compound datatype 
            INTEGER, INTENT(OUT) :: hdferr    ! Error code
            INTEGER :: flag_c

!            INTEGER, EXTERNAL :: h5pget_preserve_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_preserve_c(prp_id, flag_c)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_PRESERVE_C'::h5pget_preserve_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER ::  flag_c
              END FUNCTION h5pget_preserve_c
            END INTERFACE

            hdferr = h5pget_preserve_c(prp_id, flag_c) 
            flag = .FALSE.
            if(flag_c .eq. 1) flag = .TRUE.
          END SUBROUTINE h5pget_preserve_f

!----------------------------------------------------------------------
! Name:		h5pget_class_f 
!
! Purpose: 	Returns the property list class for a property list.
!
! Inputs:  
!		prp_id		- property list identifier
! Outputs:  
!		classtype	- property list class
!				  Possible values are:
!				  H5P_ROOT_F
!				  H5P_FILE_CREATE_F
!				  H5P_FILE_ACCESS_F
!				  H5PE_DATASET_CREATE_F
!				  H5P_DATASET_XFER_F
!				  H5P_FILE_MOUNT_F
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_class_f(prp_id, classtype, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_class_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: classtype  ! The type of the property list 
                                              ! to be created. Possible values
                                              ! are: 
                                              !  H5P_ROOT_F
                                              !  H5P_FILE_CREATE_F
                                              !  H5P_FILE_ACCESS_F
                                              !  H5PE_DATASET_CREATE_F 
                                              !  H5P_DATASET_XFER_F
                                              !  H5P_FILE_MOUNT_F
            INTEGER, INTENT(OUT) :: hdferr    ! Error code

!            INTEGER, EXTERNAL :: h5pget_class_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_class_c(prp_id, classtype)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_CLASS_C'::h5pget_class_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: classtype 
              END FUNCTION h5pget_class_c
            END INTERFACE

            hdferr = h5pget_class_c(prp_id, classtype) 
          END SUBROUTINE h5pget_class_f

!----------------------------------------------------------------------
! Name:		h5pcopy_f 
!
! Purpose: 	Copies an existing property list to create a new 
!		property list
!
! Inputs:  
!		prp_id		- property list identifier
! Outputs:  
!		new_prp_id	- new property list identifier
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pcopy_f(prp_id, new_prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pcopy_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(OUT) :: new_prp_id 
                                                ! Identifier  of property list
                                                ! copy  
            INTEGER, INTENT(OUT) :: hdferr      ! Error code

!            INTEGER, EXTERNAL :: h5pcopy_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pcopy_c(prp_id, new_prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCOPY_C'::h5pcopy_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(OUT) :: new_prp_id
              END FUNCTION h5pcopy_c
            END INTERFACE

            hdferr = h5pcopy_c(prp_id, new_prp_id)
          END SUBROUTINE h5pcopy_f

!----------------------------------------------------------------------
! Name:		h5pclose_f 
!
! Purpose: 	Terminates access to a property list. 
!
! Inputs:  
!		prp_id		- identifier of the property list to 
!				  terminate access to. 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pclose_f(prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pclose_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr        ! Error code

!            INTEGER, EXTERNAL :: h5pclose_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pclose_c(prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCLOSE_C'::h5pclose_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              END FUNCTION h5pclose_c
            END INTERFACE

            hdferr = h5pclose_c(prp_id)
          END SUBROUTINE h5pclose_f

!----------------------------------------------------------------------
! Name:		h5pset_chunk_f 
!
! Purpose: 	Sets the size of the chunks used to store 
!		a chunked layout dataset. 
!
! Inputs:  
!		prp_id		- datatset creation property list identifier
!		ndims		- number of dimensions for each chunk
!		dims		- array with dimension sizes for each chunk
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_chunk_f(prp_id, ndims, dims, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_chunk_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: ndims    ! Number of chunk dimensions
            INTEGER(HSIZE_T), DIMENSION(ndims), INTENT(IN) :: dims    
                                            ! Array containing sizes of
                                            ! chunk dimensions
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_chunk_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_chunk_c(prp_id, ndims, dims)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_CHUNK_C'::h5pset_chunk_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: ndims
              INTEGER(HSIZE_T), DIMENSION(ndims), INTENT(IN) :: dims
              END FUNCTION h5pset_chunk_c
            END INTERFACE

            hdferr =  h5pset_chunk_c(prp_id, ndims, dims)
          END SUBROUTINE h5pset_chunk_f

!----------------------------------------------------------------------
! Name:		h5pget_chunk_f 
!
! Purpose: 	Retrieves the size of chunks for the raw data of a 
!		chunked layout dataset
!
! Inputs:  
!		prp_id		- property list identifier
!		ndims		- size of dims array
! Outputs:  
!		dims		- array with dimension sizes for each chunk
!		hdferr:		- error code		
!				 	Success:  number of chunk dimensions
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_chunk_f(prp_id, ndims, dims, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_chunk_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: ndims    ! Number of chunk dimensions to
                                            ! to return
            INTEGER(HSIZE_T), DIMENSION(ndims), INTENT(OUT) :: dims    
                                            ! Array containing sizes of
                                            ! chunk dimensions
            INTEGER, INTENT(OUT) :: hdferr  ! Error code; number of
                                            ! chunk dimensions on success,
                                            ! -1 on failure

!            INTEGER, EXTERNAL :: h5pget_chunk_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_chunk_c(prp_id, ndims, dims)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_CHUNK_C'::h5pget_chunk_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER :: ndims
              INTEGER(HSIZE_T), DIMENSION(*), INTENT(OUT) :: dims
              END FUNCTION h5pget_chunk_c
            END INTERFACE

            hdferr =  h5pget_chunk_c(prp_id, ndims, dims)
          END SUBROUTINE h5pget_chunk_f

!----------------------------------------------------------------------
! Name:		h5pset_deflate_f 
!
! Purpose: 	Sets compression method and compression level. 
!
! Inputs:  
!		prp_id		- property list identifier
!		level		- compression level
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_deflate_f(prp_id, level, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_deflate_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: level        ! Compression level 
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_deflate_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_deflate_c(prp_id, level)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_DEFLATE_C'::h5pset_deflate_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: level
              END FUNCTION h5pset_deflate_c
            END INTERFACE
            hdferr = h5pset_deflate_c(prp_id, level)

          END SUBROUTINE h5pset_deflate_f

!----------------------------------------------------------------------
! Name:		h5pset(get)fill_value_f 
!
! Purpose: 	Sets(gets) fill value for a dataset creation property list
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
!		type_id		- datatype identifier for fill value
!		fillvalue	- fill value
! Outputs:  
!	(	type_id		- datatype identifier for fill value )
!	(		fillvalue	- fill value )
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:	h5pset(get)fill_value_f function is overloaded to support
!		INTEGER, REAL, DOUBLE PRECISION and CHARACTER dtatypes.		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_fill_value_integer(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_value_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype 
                                                  ! (in memory)
            INTEGER, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fill_value_integer_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fill_value_integer_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILL_VALUE_INTEGER_C'::h5pset_fill_value_integer_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              INTEGER, INTENT(IN) :: fillvalue
              END FUNCTION h5pset_fill_value_integer_c
            END INTERFACE

            hdferr = h5pset_fill_value_integer_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pset_fill_value_integer


          SUBROUTINE h5pget_fill_value_integer(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fill_value_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype
                                                  ! (in memory) 
            INTEGER, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_fill_value_integer_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fill_value_integer_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILL_VALUE_INTEGER_C'::h5pget_fill_value_integer_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              INTEGER :: fillvalue
              END FUNCTION h5pget_fill_value_integer_c
            END INTERFACE

            hdferr = h5pget_fill_value_integer_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pget_fill_value_integer


          SUBROUTINE h5pset_fill_value_real(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_value_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype 
                                                  ! (in memory)
            REAL, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fill_value_real_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fill_value_real_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILL_VALUE_REAL_C'::h5pset_fill_value_real_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              REAL, INTENT(IN) :: fillvalue
              END FUNCTION h5pset_fill_value_real_c
            END INTERFACE

            hdferr = h5pset_fill_value_real_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pset_fill_value_real


          SUBROUTINE h5pget_fill_value_real(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fill_value_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype
                                                  ! (in memory) 
            REAL, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_fill_value_real_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fill_value_real_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILL_VALUE_REAL_C'::h5pget_fill_value_real_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              REAL :: fillvalue
              END FUNCTION h5pget_fill_value_real_c
            END INTERFACE

            hdferr = h5pget_fill_value_real_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pget_fill_value_real


          SUBROUTINE h5pset_fill_value_double(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_value_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype 
                                                  ! (in memory)
            DOUBLE PRECISION, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fill_value_double_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fill_value_double_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILL_VALUE_DOUBLE_C'::h5pset_fill_value_double_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              DOUBLE PRECISION, INTENT(IN) :: fillvalue
              END FUNCTION h5pset_fill_value_double_c
            END INTERFACE

            hdferr = h5pset_fill_value_double_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pset_fill_value_double


          SUBROUTINE h5pget_fill_value_double(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fill_value_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype
                                                  ! (in memory) 
            DOUBLE PRECISION, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_fill_value_double_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fill_value_double_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILL_VALUE_DOUBLE_C'::h5pget_fill_value_double_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              DOUBLE PRECISION :: fillvalue
              END FUNCTION h5pget_fill_value_double_c
            END INTERFACE

            hdferr = h5pget_fill_value_double_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pget_fill_value_double

          SUBROUTINE h5pset_fill_value_char(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_value_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype 
                                                  ! (in memory)
            CHARACTER, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fill_valuec_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fill_valuec_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILL_VALUEC_C'::h5pset_fill_valuec_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: fillvalue 
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              CHARACTER, INTENT(IN) :: fillvalue
              END FUNCTION h5pset_fill_valuec_c
            END INTERFACE

            hdferr = h5pset_fill_valuec_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pset_fill_value_char

          SUBROUTINE h5pget_fill_value_char(prp_id, type_id, fillvalue, &
                                               hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fill_value_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: type_id ! Datatype identifier of 
                                                  ! of fillvalue datatype
                                                  ! (in memory) 
            CHARACTER, INTENT(IN) :: fillvalue   ! Fillvalue
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_fill_valuec_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fill_valuec_c(prp_id, type_id, fillvalue)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILL_VALUEC_C'::h5pget_fill_valuec_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: fillvalue 
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(IN) :: type_id
              CHARACTER :: fillvalue
              END FUNCTION h5pget_fill_valuec_c
            END INTERFACE

            hdferr = h5pget_fill_valuec_c(prp_id, type_id, fillvalue)
          END SUBROUTINE h5pget_fill_value_char

!----------------------------------------------------------------------
! Name:		h5pget_version_f 
!
! Purpose: 	Retrieves the version information of various objects 
!		for a file creation property list
!
! Inputs:  
!		prp_id		- file createion property list identifier
! Outputs:  
!		boot		- super block version number
!		freelist	- global freelist version number
!		stab		- symbol table version number
!		shhdr		- shared object header version number
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_version_f(prp_id, boot, freelist, &
                                    stab, shhdr, hdferr)
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_version_f
!DEC$endif
!

            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, DIMENSION(:), INTENT(OUT) :: boot  !array to put boot
                                                        !block version number
            INTEGER, DIMENSION(:), INTENT(OUT) :: freelist  !array to put global
                                                        !freelist version number
     
            INTEGER, DIMENSION(:), INTENT(OUT) :: stab  !array to put symbol
                                                        !table version number
            INTEGER, DIMENSION(:), INTENT(OUT) :: shhdr !array to put shared
                                                        !object header version number
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_version_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_version_c(prp_id, boot, freelist, stab, shhdr)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_VERSION_C'::h5pget_version_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, DIMENSION(:), INTENT(OUT) :: boot 
              INTEGER, DIMENSION(:), INTENT(OUT) :: freelist 
              INTEGER, DIMENSION(:), INTENT(OUT) :: stab
              INTEGER, DIMENSION(:), INTENT(OUT) :: shhdr
              END FUNCTION h5pget_version_c
            END INTERFACE

            hdferr = h5pget_version_c(prp_id, boot, freelist, stab, shhdr)
          END SUBROUTINE h5pget_version_f

!----------------------------------------------------------------------
! Name:		h5pset_userblock_f 
!
! Purpose: 	Sets user block size
!
! Inputs:  
!		prp_id		- file creation property list to modify
!		size		- size of the user-block in bytes
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------
 
          SUBROUTINE h5pset_userblock_f (prp_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_userblock_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: size !Size of the user-block in bytes 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_userblock_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_userblock_c(prp_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_USERBLOCK_C'::h5pset_userblock_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_userblock_c
            END INTERFACE

            hdferr = h5pset_userblock_c(prp_id, size)
          END SUBROUTINE h5pset_userblock_f

!----------------------------------------------------------------------
! Name:		h5pget_userblock_f 
!
! Purpose: 	Gets user block size.
!
! Inputs:  
!		prp_id		- file creation property list identifier
! Outputs:  
!		block_size	- size of the user block in bytes
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_userblock_f(prp_id, block_size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_userblock_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) ::  block_size !Size of the 
                                                               !user-block in bytes 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_userblock_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_userblock_c(prp_id, block_size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_USERBLOCK_C'::h5pget_userblock_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(OUT) :: block_size
              END FUNCTION h5pget_userblock_c
            END INTERFACE
            hdferr = h5pget_userblock_c(prp_id,  block_size)
          END SUBROUTINE h5pget_userblock_f

!----------------------------------------------------------------------
! Name:		h5pset_sizes_f 
!
! Purpose: 	Sets the byte size of the offsets and lengths used 
!		to address objects in an HDF5 file.
!
! Inputs:  
!		prp_id		- file creation property list identifier
!		sizeof_addr	- size of an object offset in bytes 
!		sizeof_size	- size of an object length in bytes
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_sizes_f (prp_id, sizeof_addr, sizeof_size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_sizes_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(SIZE_T), INTENT(IN) :: sizeof_addr !Size of an object 
                                                       !offset in bytes 
            INTEGER(SIZE_T), INTENT(IN) :: sizeof_size !Size of an object 
                                                       !length in bytes 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_sizes_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_sizes_c(prp_id, sizeof_addr, sizeof_size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SIZES_C'::h5pset_sizes_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER(SIZE_T), INTENT(IN) :: sizeof_addr
              INTEGER(SIZE_T), INTENT(IN) :: sizeof_size
              END FUNCTION h5pset_sizes_c
            END INTERFACE

            hdferr = h5pset_sizes_c(prp_id, sizeof_addr, sizeof_size)
          END SUBROUTINE h5pset_sizes_f

!----------------------------------------------------------------------
! Name:		h5pget_sizes_f 
!
! Purpose: 	Retrieves the size of the offsets and lengths used 
!		in an HDF5 file
!
! Inputs:  
!		prp_id		- file creation property list identifier
! Outputs:  
!		sizeof_addr	- size of an object offset in bytes 
!		sizeof_size	- size of an object length in bytes
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_sizes_f(prp_id, sizeof_addr, sizeof_size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_sizes_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(SIZE_T), INTENT(OUT) :: sizeof_addr !Size of an object
                                                                      !offset in bytes 
            INTEGER(SIZE_T), INTENT(OUT) :: sizeof_size !Size of an object
                                                                      !length in bytes 
  
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_sizes_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_sizes_c(prp_id, sizeof_addr, sizeof_size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_SIZES_C'::h5pget_sizes_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER(SIZE_T), INTENT(OUT) :: sizeof_addr
              INTEGER(SIZE_T), INTENT(OUT) :: sizeof_size
              END FUNCTION h5pget_sizes_c
            END INTERFACE

            hdferr = h5pget_sizes_c(prp_id, sizeof_addr, sizeof_size)
          END SUBROUTINE h5pget_sizes_f

!----------------------------------------------------------------------
! Name:		h5pset_sym_k_f 
!
! Purpose: 	Sets the size of parameters used to control the 
!		symbol table nodes
!
! Inputs:  
!		prp_id		- file creation property list identifier
!		ik		- symbol table tree rank
!		lk		- symbol table node size
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_sym_k_f (prp_id, ik, lk, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_sym_k_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: ik ! Symbol table tree rank 
            INTEGER, INTENT(IN) :: lk ! Symbol table node size 
                                                       
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_sym_k_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_sym_k_c(prp_id, ik, lk)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SYM_K_C'::h5pset_sym_k_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(IN) :: ik
              INTEGER, INTENT(IN) :: lk
              END FUNCTION h5pset_sym_k_c
            END INTERFACE

            hdferr = h5pset_sym_k_c(prp_id, ik, lk)
          END SUBROUTINE h5pset_sym_k_f

!----------------------------------------------------------------------
! Name:		h5pget_sym_k_f 
!
! Purpose: 	Retrieves the size of the symbol table B-tree 1/2 rank
!		 and the symbol table leaf node 1/2 size. 
!
! Inputs:  
!		prp_id		- file creation property list identifier
! Outputs:  
!		ik		- symbol table tree 1/2 rank
!		lk		- symbol table node 1/2 size
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_sym_k_f(prp_id, ik, lk, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_sym_k_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: ik !Symbol table tree rank
            INTEGER, INTENT(OUT) :: lk !Symbol table node size
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_sym_k_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_sym_k_c(prp_id, ik, lk)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_SYM_K_C'::h5pget_sym_k_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(OUT) :: ik
              INTEGER, INTENT(OUT) :: lk
              END FUNCTION h5pget_sym_k_c
            END INTERFACE

            hdferr = h5pget_sym_k_c(prp_id, ik, lk)
          END SUBROUTINE h5pget_sym_k_f

!----------------------------------------------------------------------
! Name:		h5pset_istore_k_f 
!
! Purpose: 	Sets the size of the parameter used to control the 
!		B-trees for indexing chunked datasets
!
! Inputs:  
!		prp_id		- file creation property list identifier
!		ik		- 1/2 rank of chunked storage B-tree
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_istore_k_f (prp_id, ik, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_istore_k_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: ik ! 1/2 rank of chunked storage B-tree
                                                       
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_istore_k_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_istore_k_c(prp_id, ik)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_ISTORE_K_C'::h5pset_istore_k_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: ik
              END FUNCTION h5pset_istore_k_c
            END INTERFACE

            hdferr = h5pset_istore_k_c(prp_id, ik)
          END SUBROUTINE h5pset_istore_k_f

!----------------------------------------------------------------------
! Name:		h5pget_istore_k_f 
!
! Purpose: 	Queries the 1/2 rank of an indexed storage B-tree. 
!
! Inputs:  
!		prp_id		- file creation property list identifier
! Outputs:  
!		ik		- 1/2 rank of chunked storage B-tree
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_istore_k_f(prp_id, ik, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_istore_k_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: ik !1/2 rank of chunked storage B-tree
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_istore_k_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_istore_k_c(prp_id, ik)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_ISTORE_K_C'::h5pget_istore_k_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: ik
              END FUNCTION h5pget_istore_k_c
            END INTERFACE

            hdferr = h5pget_istore_k_c(prp_id, ik)
          END SUBROUTINE h5pget_istore_k_f

!----------------------------------------------------------------------
! Name:		h5pget_driver_f 
!
! Purpose: 	Returns low-lever driver identifier. 
!
! Inputs:  
!		prp_id		- file access or data transfer property 
!				  list identifier. 
! Outputs:  
!		driver		- low-level driver identifier
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_driver_f(prp_id, driver, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_driver_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HID_T), INTENT(OUT) :: driver !low-level file driver identifier
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_driver_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_driver_c(prp_id, driver)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_DRIVER_C'::h5pget_driver_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(OUT) :: driver
              END FUNCTION h5pget_driver_c
            END INTERFACE

            hdferr = h5pget_driver_c(prp_id, driver)
          END SUBROUTINE h5pget_driver_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_stdio_f 
!
! Purpose: 	Sets the standard I/O driver. 
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fapl_stdio_f (prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_stdio_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fapl_stdio_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_stdio_c(prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_STDIO_C'::h5pset_fapl_stdio_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              END FUNCTION h5pset_fapl_stdio_c
            END INTERFACE

            hdferr = h5pset_fapl_stdio_c(prp_id)
          END SUBROUTINE h5pset_fapl_stdio_f

!----------------------------------------------------------------------
! Name:		h5pget_stdio_f 
!
! Purpose:  	NOT AVAILABLE
!
! Inputs:  
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

!          SUBROUTINE h5pget_stdio_f (prp_id, io, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_stdio_f
!DEC$endif
!
!            IMPLICIT NONE
!            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
!            INTEGER, INTENT(OUT) :: io   ! value indicates that the file 
                                         !access property list is set to 
                                         !the stdio driver
!            INTEGER, INTENT(OUT) :: hdferr  ! Error code
!            INTEGER, EXTERNAL :: h5pget_stdio_c
!            hdferr = h5pget_stdio_c(prp_id, io)
!          END SUBROUTINE h5pget_stdio_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_sec2_f 
!
! Purpose: 	Sets the sec2 driver. 
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fapl_sec2_f (prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_sec2_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fapl_sec2_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_sec2_c(prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_SEC2_C'::h5pset_fapl_sec2_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
              END FUNCTION h5pset_fapl_sec2_c
            END INTERFACE

            hdferr = h5pset_fapl_sec2_c(prp_id)
          END SUBROUTINE h5pset_fapl_sec2_f

!----------------------------------------------------------------------
! Name:		h5pget_sec2_f 
!
! Purpose: 	NOT AVAILABLE
!
! Inputs:  
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

!          SUBROUTINE h5pget_sec2_f (prp_id, sec2, hdferr) 
!            IMPLICIT NONE
!            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
!            INTEGER, INTENT(OUT) :: sec2   ! value indicates whether the file 
                                           !driver uses the functions declared
                                           !in the unistd.h file
!            INTEGER, INTENT(OUT) :: hdferr  ! Error code
!            INTEGER, EXTERNAL :: h5pget_sec2_c
!            hdferr = h5pget_sec2_c(prp_id, sec2)
!          END SUBROUTINE h5pget_sec2_f

!----------------------------------------------------------------------
! Name:		h5pset_alignment_f 
!
! Purpose: 	Sets alignment properties of a file access property list. 
!
! Inputs:  
!		prp_id		- file access property list identifier
!		threshold	- threshold value	
!		alignment	- alignment value
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_alignment_f(prp_id, threshold,  alignment, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_alignment_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: threshold ! Threshold value
            INTEGER(HSIZE_T), INTENT(IN) :: alignment ! alignment value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_alignment_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_alignment_c(prp_id, threshold, alignment)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_ALIGNMENT_C'::h5pset_alignment_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(IN) :: threshold
              INTEGER(HSIZE_T), INTENT(IN) :: alignment
              END FUNCTION h5pset_alignment_c
            END INTERFACE

            hdferr = h5pset_alignment_c(prp_id, threshold, alignment)
          END SUBROUTINE h5pset_alignment_f

!----------------------------------------------------------------------
! Name:		h5pget_alignment_f 
!
! Purpose: 	Retrieves the current settings for alignment 
!		properties from a file access property list. 
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		threshold	- threshold value	
!		alignment	- alignment value
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_alignment_f(prp_id, threshold,  alignment, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_alignment_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) :: threshold ! Threshold value
            INTEGER(HSIZE_T), INTENT(OUT) :: alignment ! alignment value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_alignment_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_alignment_c(prp_id, threshold, alignment)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_ALIGNMENT_C'::h5pget_alignment_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(OUT) :: threshold
              INTEGER(HSIZE_T), INTENT(OUT) :: alignment
              END FUNCTION h5pget_alignment_c
            END INTERFACE

            hdferr = h5pget_alignment_c(prp_id, threshold, alignment)
          END SUBROUTINE h5pget_alignment_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_core_f 
!
! Purpose: 	Modifies the file access property list to use the 
!		H5FD_CORE driver. 
!
! Inputs:  	prp_id		- file access property list identifier
!		increment	- size, in bytes, of memory increments 
!		backing_store	- boolean flag indicating whether to write 
!				  the file contents to disk when the file is closed. 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fapl_core_f(prp_id, increment, backing_store, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_core_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(SIZE_T), INTENT(IN) :: increment ! File block size in bytes.
            LOGICAL, INTENT(IN) :: backing_store ! flag to indicate that
                                    ! entire file contents are flushed to a file 
                                    ! with the same name as this core file.
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: backing_store_flag 

!            INTEGER, EXTERNAL :: h5pset_fapl_core_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_core_c(prp_id, increment, backing_store_flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_CORE_C'::h5pset_fapl_core_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER(SIZE_T), INTENT(IN) :: increment 
              INTEGER :: backing_store_flag 
              END FUNCTION h5pset_fapl_core_c
            END INTERFACE
            backing_store_flag = 0
            if(backing_store) backing_store_flag = 1
            hdferr = h5pset_fapl_core_c(prp_id, increment, backing_store_flag)
          END SUBROUTINE h5pset_fapl_core_f

!----------------------------------------------------------------------
! Name:		h5pget_fapl_core_f 
!
! Purpose: 	Queries core file driver properties. 
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		increment	- size, in bytes, of memory increments 
!		backing_store	- boolean flag indicating whether to write 
!				  the file contents to disk when the file is closed. 
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_fapl_core_f(prp_id, increment, backing_store, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fapl_core_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(SIZE_T), INTENT(OUT) :: increment ! File block size in bytes.
            LOGICAL, INTENT(OUT) :: backing_store ! flag to indicate that
                                    ! entire file contents are flushed to a file 
                                    ! with the same name as this core file.
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: backing_store_flag 

!            INTEGER, EXTERNAL :: h5pget_fapl_core_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fapl_core_c(prp_id, increment, backing_store_flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FAPL_CORE_C'::h5pget_fapl_core_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER(SIZE_T), INTENT(OUT) :: increment 
              INTEGER :: backing_store_flag 
              END FUNCTION h5pget_fapl_core_c
            END INTERFACE

            hdferr = h5pget_fapl_core_c(prp_id, increment, backing_store_flag)
            backing_store =.FALSE.
            if (backing_store_flag .eq. 1) backing_store =.TRUE.
          END SUBROUTINE h5pget_fapl_core_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_family_f 
!
! Purpose: 	Sets the file access property list to use the family driver. 
!
! Inputs:  
!		prp_id		- file access property list identifier
!		memb_size	- size in bytes of each file member 
!		memb_plist	- identifier of the file access property 
!				  list to be used for each family member
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fapl_family_f(prp_id, memb_size, memb_plist , hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_family_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: memb_size ! Logical size, in bytes,
                                                      !of each family member
            INTEGER(HID_T), INTENT(IN) :: memb_plist !Identifier of the file 
                                                     !access property list for 
                                                     !each member of the family
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_fapl_family_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_family_c(prp_id, memb_size, memb_plist)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_FAMILY_C'::h5pset_fapl_family_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(IN) :: memb_size
              INTEGER(HID_T), INTENT(IN) :: memb_plist
              END FUNCTION h5pset_fapl_family_c
            END INTERFACE

            hdferr = h5pset_fapl_family_c(prp_id, memb_size, memb_plist)
          END SUBROUTINE h5pset_fapl_family_f

!----------------------------------------------------------------------
! Name:		h5pget_fapl_family_f 
!
! Purpose:	Returns file access property list information.  	
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		memb_size	- size in bytes of each file member 
!		memb_plist	- identifier of the file access property 
!				  list to be used for each family member
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_fapl_family_f(prp_id, memb_size, memb_plist , hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fapl_family_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) :: memb_size ! Logical size, in bytes,
                                                      !of each family member
            INTEGER(HID_T), INTENT(OUT) :: memb_plist !Identifier of the file 
                                                     !access property list for 
                                                     !each member of the family
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_fapl_family_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fapl_family_c(prp_id, memb_size, memb_plist)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FAPL_FAMILY_C'::h5pget_fapl_family_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(OUT) :: memb_size
              INTEGER(HID_T), INTENT(OUT) :: memb_plist
              END FUNCTION h5pget_fapl_family_c
            END INTERFACE

            hdferr = h5pget_fapl_family_c(prp_id, memb_size, memb_plist)
          END SUBROUTINE h5pget_fapl_family_f

!----------------------------------------------------------------------
! Name:		h5pset_cache_f 
!
! Purpose: 	Sets the meta data cache and raw data chunk 
!		cache parameters
!
! Inputs:  
!		prp_id		- file access property list identifier
!		mdc_nelmts	- number of elements (objects) in the meta 
!				  data cache 
!		rdcc_nelmts	- number of elements (objects) in the raw 
!			          data chunk cache 
!		rdcc_nbytes	- total size of the raw data chunk cache, in bytes 
!		rdcc_w0		- preemption policy (0 or 1)
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_cache_f(prp_id, mdc_nelmts,rdcc_nelmts, rdcc_nbytes, rdcc_w0, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_cache_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: mdc_nelmts  !Number of elements (objects)
                                                        ! in the meta data cache
            INTEGER(SIZE_T), INTENT(IN) :: rdcc_nelmts  !Number of elements (objects)
                                                        ! in the meta data cache
            INTEGER(SIZE_T), INTENT(IN) :: rdcc_nbytes !Total size of the raw data 
                                                      !chunk cache, in bytes 
            REAL, INTENT(IN) :: rdcc_w0 !Preemption policy
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_cache_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_cache_c(prp_id,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_CACHE_C'::h5pset_cache_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: mdc_nelmts 
              INTEGER(SIZE_T), INTENT(IN) :: rdcc_nelmts 
              INTEGER(SIZE_T), INTENT(IN) :: rdcc_nbytes
              REAL, INTENT(IN) :: rdcc_w0
              END FUNCTION h5pset_cache_c
            END INTERFACE

            hdferr = h5pset_cache_c(prp_id, mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0 )
          END SUBROUTINE h5pset_cache_f

!----------------------------------------------------------------------
! Name:		h5pget_cache_f 
!
! Purpose: 	Queries the meta data cache and raw data chunk cache 
!		parameters.  
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		mdc_nelmts	- number of elements (objects) in the meta 
!				  data cache 
!		rdcc_nelmts	- number of elements (objects) in the raw 
!			          data chunk cache 
!		rdcc_nbytes	- total size of the raw data chunk cache, in bytes 
!		rdcc_w0		- preemption policy (0 or 1)
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
!                       Bug fix: type of the rdcc_nelmts parameter should be INTEGER
!                                instead of INTEGER(SIZE_T) October 10, 2003 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_cache_f(prp_id, mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_cache_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: mdc_nelmts  !Number of elements (objects)
                                                        ! in the meta data cache
            INTEGER(SIZE_T), INTENT(OUT) :: rdcc_nelmts  !Number of elements (objects)
                                                        ! in the meta data cache
            INTEGER(SIZE_T), INTENT(OUT) :: rdcc_nbytes !Total size of the raw data 
                                                      !chunk cache, in bytes 
            REAL, INTENT(OUT) :: rdcc_w0 !Preemption policy
            INTEGER, INTENT(OUT) :: hdferr  ! Error code


!            INTEGER, EXTERNAL :: h5pget_cache_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_cache_c(prp_id,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_CACHE_C'::h5pget_cache_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: mdc_nelmts 
              INTEGER(SIZE_T), INTENT(OUT) :: rdcc_nelmts 
              INTEGER(SIZE_T), INTENT(OUT) :: rdcc_nbytes
              REAL, INTENT(OUT) :: rdcc_w0
              END FUNCTION h5pget_cache_c
            END INTERFACE

            hdferr = h5pget_cache_c(prp_id, mdc_nelmts,rdcc_nelmts, rdcc_nbytes, rdcc_w0 )
          END SUBROUTINE h5pget_cache_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_split_f 
!
! Purpose: 	Emulates the old split file driver. 
!
! Inputs:  
!		prp_id		- file access property list identifier
!		meta_ext	- name of the extension for the metafile 
!				  filename
!		meta_plist	- identifier of the meta file access property 
!				  list
!		raw_ext 	- name extension for the raw file filename
!		raw_plist	- identifier of the raw file access property list
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fapl_split_f(prp_id, meta_ext, meta_plist, raw_ext, raw_plist, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_split_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: meta_ext  !Name of the extension for
                                                      !the metafile filename
            INTEGER(HID_T), INTENT(IN) :: meta_plist  ! Identifier of the meta file
                                                      ! access property list
            CHARACTER(LEN=*), INTENT(IN) :: raw_ext  !Name extension for the raw file filename
            INTEGER(HID_T), INTENT(IN) :: raw_plist  !Identifier of the raw file 
                                                     !access property list
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: meta_len, raw_len

!            INTEGER, EXTERNAL :: h5pset_fapl_split_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_split_c(prp_id,meta_len,meta_ext,meta_plist,raw_len,raw_ext,raw_plist)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_SPLIT_C'::h5pset_fapl_split_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: meta_ext
              !DEC$ATTRIBUTES reference :: raw_ext
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: meta_ext 
              INTEGER(HID_T), INTENT(IN) :: meta_plist
              CHARACTER(LEN=*), INTENT(IN) :: raw_ext
              INTEGER(HID_T), INTENT(IN) :: raw_plist 
              INTEGER :: meta_len, raw_len
              END FUNCTION h5pset_fapl_split_c
            END INTERFACE

            meta_len = LEN(meta_ext)
            raw_len = LEN(raw_ext)
        hdferr = h5pset_fapl_split_c(prp_id,meta_len,meta_ext,meta_plist,raw_len,raw_ext,raw_plist)
          END SUBROUTINE h5pset_fapl_split_f

!----------------------------------------------------------------------
! Name:		h5pget_split_f 
!
! Purpose: 	NOT AVAILABLE
!
! Inputs:  
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

!          SUBROUTINE h5pget_split_f(prp_id, meta_ext_size, meta_ext, meta_plist,raw_ext_size,&
!                                     raw_ext, raw_plist, hdferr) 
!            IMPLICIT NONE
!            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier
!            INTEGER(SIZE_T), INTENT(IN) :: meta_ext_size ! Number of characters of the meta
                                                         ! file extension to be copied to the
                                                         ! meta_ext buffer
      
!            CHARACTER(LEN=*), INTENT(OUT) :: meta_ext  !Name of the extension for
                                                      !the metafile filename
!            INTEGER(HID_T), INTENT(OUT) :: meta_plist  ! Identifier of the meta file
                                                      ! access property list
!            INTEGER(SIZE_T), INTENT(IN) :: raw_ext_size ! Number of characters of the raw
                                                         ! file extension to be copied to the
                                                         ! raw_ext buffer
!            CHARACTER(LEN=*), INTENT(OUT) :: raw_ext  !Name extension for the raw file filename
!            INTEGER(HID_T), INTENT(OUT) :: raw_plist  !Identifier of the raw file 
                                                     !access property list
!            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_split_c
!            hdferr = h5pget_split_c(prp_id, meta_ext_size, meta_ext, meta_plist, &
!                                    raw_ext_size, raw_ext, raw_plist )
!          END SUBROUTINE h5pget_split_f

!----------------------------------------------------------------------
! Name:		h5pset_gc_references_f 
!
! Purpose: 	Sets garbage collecting references flag. 
!
! Inputs:  
!		prp_id		- file access property list identifier
!		gc_reference	- flag for stting garbage collection on 
!				  and off (1 or 0)
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

 
          SUBROUTINE h5pset_gc_references_f (prp_id, gc_reference, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_gc_references_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: gc_reference !the flag for garbage collecting
                                                ! references for the file
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_gc_references_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_gc_references_c(prp_id, gc_reference)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_GC_REFERENCES_C'::h5pset_gc_references_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: gc_reference
              END FUNCTION h5pset_gc_references_c
            END INTERFACE

            hdferr = h5pset_gc_references_c(prp_id, gc_reference)
          END SUBROUTINE h5pset_gc_references_f

!----------------------------------------------------------------------
! Name:		h5pget_gc_references_f 
!
! Purpose: 	Returns garbage collecting references setting. 	
!
! Inputs:  
!		prp_id		- file access property list identifier
! Outputs:  
!		gc_reference	- flag for stting garbage collection on 
!				  and off (1 or 0)
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_gc_references_f (prp_id, gc_reference, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_gc_references_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: gc_reference !the flag for garbage collecting
                                                ! references for the file
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_gc_references_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_gc_references_c(prp_id, gc_reference)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_GC_REFERENCES_C'::h5pget_gc_references_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: gc_reference
              END FUNCTION h5pget_gc_references_c
            END INTERFACE

            hdferr = h5pget_gc_references_c(prp_id, gc_reference)
          END SUBROUTINE h5pget_gc_references_f

!----------------------------------------------------------------------
! Name:		h5pset_layout_f 
!
! Purpose: 	Sets the type of storage used store the raw data 
!		for a dataset. 
!
! Inputs:  
!		prp_id		- data creation property list identifier
!		layout		- type of storage layout for raw data
!				  possible values are:
!				  H5D_COMPACT_F
!				  H5D_CONTIGUOUS_F
!				  H5D_CHUNKED_F
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_layout_f (prp_id, layout, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_layout_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: layout !Type of storage layout for raw data
                                          !possible values are:
                                          !H5D_COMPACT_F
                                          !H5D_CONTIGUOUS_F
                                          !H5D_CHUNKED_F
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_layout_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_layout_c(prp_id, layout)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_LAYOUT_C'::h5pset_layout_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: layout
              END FUNCTION h5pset_layout_c
            END INTERFACE

            hdferr = h5pset_layout_c(prp_id, layout)
          END SUBROUTINE h5pset_layout_f

!----------------------------------------------------------------------
! Name:		h5pget_layout_f 
!
! Purpose: 	Returns the layout of the raw data for a dataset. 
!
! Inputs:  
!		prp_id		- data creation property list identifier
! Outputs:  
!		layout		- type of storage layout for raw data
!				  possible values are:
!				  H5D_COMPACT_F
!				  H5D_CONTIGUOUS_F
!				  H5D_CHUNKED_F
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_layout_f (prp_id, layout, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_layout_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: layout !Type of storage layout for raw data
                                          !possible values are:
                                          !H5D_COMPACT_F(0)
                                          !H5D_CONTIGUOUS_F(1)
                                          !H5D_CHUNKED_F(2)
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_layout_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_layout_c(prp_id, layout)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_LAYOUT_C'::h5pget_layout_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: layout
              END FUNCTION h5pget_layout_c
            END INTERFACE

            hdferr = h5pget_layout_c(prp_id, layout)
          END SUBROUTINE h5pget_layout_f

!----------------------------------------------------------------------
! Name:		h5pset_filter_f 
!
! Purpose: 	Adds a filter to the filter pipeline. 
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
!		filter		- filter to be added to the pipeline 
!		flags		- bit vector specifying certain general
!				  properties of the filter
!		cd_nelmts	- number of elements in cd_values
!		cd_values	- auxiliary data for the filter
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		February, 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_filter_f(prp_id, filter, flags, cd_nelmts, cd_values,  hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_filter_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: filter  !Filter to be added to the pipeline.
            INTEGER, INTENT(IN) :: flags  !Bit vector specifying certain general
                                          !properties of the filter.
            INTEGER(SIZE_T), INTENT(IN) :: cd_nelmts  !Number of elements in cd_values.
            INTEGER, DIMENSION(*), INTENT(IN) :: cd_values  !Auxiliary data for the filter.

            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_filter_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_filter_c(prp_id, filter, flags, cd_nelmts, cd_values)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILTER_C'::h5pset_filter_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(IN) :: filter 
              INTEGER, INTENT(IN) :: flags 
              INTEGER(SIZE_T), INTENT(IN) :: cd_nelmts 
              INTEGER, DIMENSION(*), INTENT(IN) :: cd_values 
              END FUNCTION h5pset_filter_c
            END INTERFACE

            hdferr = h5pset_filter_c(prp_id, filter, flags, cd_nelmts, cd_values )
          END SUBROUTINE h5pset_filter_f

!----------------------------------------------------------------------
! Name:		h5pget_nfilters_f 
!
! Purpose: 	Returns the number of filters in the pipeline. 
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
! Outputs:  
!		nfilters	- number of filters in the pipeline
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_nfilters_f (prp_id, nfilters, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_nfilters_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: nfilters !the number of filters in the pipeline
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_nfilters_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_nfilters_c(prp_id, nfilters)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_NFILTERS_C'::h5pget_nfilters_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: nfilters
              END FUNCTION h5pget_nfilters_c
            END INTERFACE

            hdferr = h5pget_nfilters_c(prp_id, nfilters)
          END SUBROUTINE h5pget_nfilters_f

!----------------------------------------------------------------------
! Name:		h5pget_filter_f 
!
! Purpose: 	Returns information about a filter in a pipeline
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
! Outputs:  
!				  identifier
!		filter		- filter to be added to the pipeline 
!		flags		- bit vector specifying certain general
!				  properties of the filter
!		cd_nelmts	- number of elements in cd_values
!		cd_values	- auxiliary data for the filter
!		namelen		- number of characters in the name buffer
!		name		- buffer to retrieve filter name
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_filter_f(prp_id, filter_number, flags, cd_nelmts, cd_values, namelen, name, filter_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_filter_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: filter_number  !Sequence number within the filter
                                                  !pipeline of the filter for which 
                                                  !information is sought
            INTEGER, DIMENSION(*), INTENT(OUT) :: cd_values  !Auxiliary data for the filter.
            INTEGER, INTENT(OUT) :: flags  !Bit vector specifying certain general
                                          !properties of the filter.
            INTEGER(SIZE_T), INTENT(INOUT) :: cd_nelmts  !Number of elements in cd_values.
            INTEGER(SIZE_T), INTENT(IN) :: namelen !Anticipated number of characters in name.
            CHARACTER(LEN=*), INTENT(OUT) :: name !Name of the filter
            INTEGER, INTENT(OUT) :: filter_id ! filter identification number  

            INTEGER, INTENT(OUT) :: hdferr  ! Error code


!            INTEGER, EXTERNAL :: h5pget_filter_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_filter_c(prp_id, filter_number, flags, cd_nelmts,  &
                                              cd_values, namelen, name, filter_id )
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILTER_C'::h5pget_filter_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: filter_number 
              INTEGER, DIMENSION(*), INTENT(OUT) :: cd_values
              INTEGER, INTENT(OUT) :: flags 
              INTEGER(SIZE_T), INTENT(INOUT) :: cd_nelmts
              INTEGER(SIZE_T), INTENT(IN) :: namelen
              CHARACTER(LEN=*), INTENT(OUT) :: name
              INTEGER, INTENT(OUT) :: filter_id
              END FUNCTION h5pget_filter_c
            END INTERFACE

            hdferr = h5pget_filter_c(prp_id, filter_number, flags, cd_nelmts,  & 
                                     cd_values, namelen, name, filter_id )
          END SUBROUTINE h5pget_filter_f

!----------------------------------------------------------------------
! Name:		h5pset_external_f 
!
! Purpose: 	Adds an external file to the list of external files. 
!
! Inputs:   
!		prp_id		- dataset creation property list identifier
!		name		- name of external file
!		offset		- offset in bytes from the beginning of the 
!				  file to the location in the file
!				  where the data starts
!		bytes		- size of the external file data. 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_external_f(prp_id, name, offset,bytes, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_external_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier
            CHARACTER(LEN=*), INTENT(IN) :: name !Name of an external file
            INTEGER, INTENT(IN) :: offset !Offset, in bytes, from the beginning 
                                          !of the file to the location in the file 
                                          !where the data starts.
            INTEGER(HSIZE_T), INTENT(IN) :: bytes ! Number of bytes reserved in the 
                                                 !file for the data
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

            INTEGER :: namelen

!            INTEGER, EXTERNAL :: h5pset_external_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_external_c(prp_id, name,namelen, offset, bytes)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_EXTERNAL_C'::h5pset_external_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER :: namelen
              INTEGER, INTENT(IN) :: offset
              INTEGER(HSIZE_T), INTENT(IN) :: bytes
              END FUNCTION h5pset_external_c
            END INTERFACE

            namelen = LEN(name)
            hdferr = h5pset_external_c(prp_id, name, namelen, offset, bytes)
          END SUBROUTINE h5pset_external_f

!----------------------------------------------------------------------
! Name:		h5pget_external_count_f 
!
! Purpose: 	Returns the number of external files for a dataset. 
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
! Outputs:  
!		count		- number of external files for the 
!				  specified dataset
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_external_count_f (prp_id, count, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_external_count_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: count !number of external files for the 
                                          !specified dataset
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
!            INTEGER, EXTERNAL :: h5pget_external_count_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_external_count_c(prp_id, count)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_EXTERNAL_COUNT_C'::h5pget_external_count_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(OUT) :: count
              END FUNCTION h5pget_external_count_c
            END INTERFACE

            hdferr = h5pget_external_count_c(prp_id, count)
          END SUBROUTINE h5pget_external_count_f

!----------------------------------------------------------------------
! Name:		h5pget_external_f 
!
! Purpose: 	Returns information about an external file. 
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
! Outputs:  
!		idx		- external file index 
!		name_size	- maximum size of name array
!		name		- name of the external file	
!		name		- name of external file
!	 	offset		- offset in bytes from the beginning of the 
!				  file to the location in the file
!				  where the data starts
!		bytes		- size of the external file data
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_external_f(prp_id, idx, name_size, name, offset,bytes, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_external_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier
            INTEGER, INTENT(IN) :: idx !External file index.
            INTEGER(SIZE_T), INTENT(IN) :: name_size !Maximum length of name array 
            CHARACTER(LEN=*), INTENT(OUT) :: name !Name of an external file
            INTEGER, INTENT(OUT) :: offset !Offset, in bytes, from the beginning 
                                          !of the file to the location in the file 
                                          !where the data starts.
            INTEGER(HSIZE_T), INTENT(OUT) :: bytes ! Number of bytes reserved in the 
                                                 !file for the data
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pget_external_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_external_c(prp_id, idx, name_size, name, offset, bytes)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_EXTERNAL_C'::h5pget_external_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: idx 
              INTEGER(SIZE_T), INTENT(IN) :: name_size
              CHARACTER(LEN=*), INTENT(OUT) :: name
              INTEGER, INTENT(OUT) :: offset
              INTEGER(HSIZE_T), INTENT(OUT) :: bytes
              END FUNCTION h5pget_external_c
            END INTERFACE

            hdferr = h5pget_external_c(prp_id, idx, name_size, name, offset, bytes)
          END SUBROUTINE h5pget_external_f

!----------------------------------------------------------------------
! Name:		h5pset_btree_ratios_f 
!
! Purpose: 	Sets B-tree split ratios for a dataset transfer 
!		property list. 
!
! Inputs:  	
!		prp_id		- the dataset transfer property list 
!				  identifier 
!		left		- the B-tree split ratio for left-most nodes 
!		middle		- the B-tree split ratio for all other nodes
!		right		- the B-tree split ratio for right-most nodes
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_btree_ratios_f(prp_id, left, middle, right, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_btree_ratios_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier
            REAL, INTENT(IN) :: left !The B-tree split ratio for left-most nodes.
            REAL, INTENT(IN) :: middle !The B-tree split ratio for all other nodes 
            REAL, INTENT(IN) :: right !The B-tree split ratio for right-most 
                                      !nodes and lone nodes. 

            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pset_btree_ratios_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION  h5pset_btree_ratios_c(prp_id, left, middle, right)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_BTREE_RATIOS_C'::h5pset_btree_ratios_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              REAL, INTENT(IN) :: left
              REAL, INTENT(IN) :: middle
              REAL, INTENT(IN) :: right
              END FUNCTION h5pset_btree_ratios_c
            END INTERFACE

            hdferr = h5pset_btree_ratios_c(prp_id, left, middle, right)
          END SUBROUTINE h5pset_btree_ratios_f

!----------------------------------------------------------------------
! Name:		h5pget_btree_ratios_f 
!
! Purpose: 	Gets B-tree split ratios for a dataset transfer property list
!
! Inputs:  
!		prp_id		- the dataset transfer property list 
!				  identifier 
! Outputs:  
!		left		- the B-tree split ratio for left-most nodes 
!		middle		- the B-tree split ratio for all other nodes
!		right		- the B-tree split ratio for right-most nodes
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		August 12, 1999	
!
! Modifications: 	Explicit Fortran interfaces were added for 
!			called C functions (it is needed for Windows
!			port).  March 14, 2001 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_btree_ratios_f(prp_id, left, middle, right, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_btree_ratios_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier
            REAL, INTENT(OUT) :: left !The B-tree split ratio for left-most nodes.
            REAL, INTENT(OUT) :: middle !The B-tree split ratio for all other nodes 
            REAL, INTENT(OUT) :: right !The B-tree split ratio for right-most 
                                      !nodes and lone nodes. 

            INTEGER, INTENT(OUT) :: hdferr  ! Error code


!            INTEGER, EXTERNAL :: h5pget_btree_ratios_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION  h5pget_btree_ratios_c(prp_id, left, middle, right)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_BTREE_RATIOS_C'::h5pget_btree_ratios_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              REAL, INTENT(OUT) :: left
              REAL, INTENT(OUT) :: middle
              REAL, INTENT(OUT) :: right
              END FUNCTION h5pget_btree_ratios_c
            END INTERFACE

            hdferr = h5pget_btree_ratios_c(prp_id, left, middle, right)
          END SUBROUTINE h5pget_btree_ratios_f

!----------------------------------------------------------------------
! Name:		h5pget_fclose_degree_f 
!
! Purpose: 	Returns the degree for the file close behavior.
!
! Inputs:  
!		fapl_id		- file access property list identifier
! Outputs:  
!		degree  	- one of the following:
!				  Possible values are:
!				  H5F_CLOSE_DEFAULT_F
!				  H5F_CLOSE_WEAK_F
!				  H5F_CLOSE_SEMI_F
!				  H5F_CLOSE_STRONG_F
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        September 26, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_fclose_degree_f(fapl_id, degree, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fclose_degree_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: fapl_id ! File Access Property list identifier 
            INTEGER, INTENT(OUT) :: degree     ! Possible values
                                              ! are: 
						!  H5F_CLOSE_DEFAULT_F
						!  H5F_CLOSE_WEAK_F
						!  H5F_CLOSE_SEMI_F
						!  H5F_CLOSE_STRONG_F

            INTEGER, INTENT(OUT) :: hdferr    ! Error code

!            INTEGER, EXTERNAL :: h5pget_fclose_degree_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fclose_degree_c(fapl_id, degree)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
         !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FCLOSE_DEGREE_C'::h5pget_fclose_degree_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: fapl_id
              INTEGER, INTENT(OUT) :: degree
              END FUNCTION h5pget_fclose_degree_c
            END INTERFACE

            hdferr = h5pget_fclose_degree_c(fapl_id, degree) 
          END SUBROUTINE h5pget_fclose_degree_f

!----------------------------------------------------------------------
! Name:		h5pset_fclose_degree_f 
!
! Purpose: 	Sets the degree for the file close behavior.
!
! Inputs:  
!		fapl_id		- file access property list identifier
!		degree  	- one of the following:
!				  Possible values are:
!				  H5F_CLOSE_DEFAULT_F
!				  H5F_CLOSE_WEAK_F
!				  H5F_CLOSE_SEMI_F
!				  H5F_CLOSE_STRONG_F
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        September 26, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fclose_degree_f(fapl_id, degree, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fclose_degree_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: fapl_id ! File Access Property list identifier 
            INTEGER, INTENT(IN) :: degree     ! Possible values
                                              ! are: 
						!  H5F_CLOSE_DEFAULT_F
						!  H5F_CLOSE_WEAK_F
						!  H5F_CLOSE_SEMI_F
						!  H5F_CLOSE_STRONG_F

            INTEGER, INTENT(OUT) :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_fclose_degree_c(fapl_id, degree)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FCLOSE_DEGREE_C'::h5pset_fclose_degree_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: fapl_id
              INTEGER, INTENT(IN) :: degree
              END FUNCTION h5pset_fclose_degree_c
            END INTERFACE

            hdferr = h5pset_fclose_degree_c(fapl_id, degree) 
          END SUBROUTINE h5pset_fclose_degree_f

!----------------------------------------------------------------------
! Name:		h5pequal_f 
!
! Purpose: 	Checks if two property lists are eqaul
!
! Inputs:  
!		plist1_id	- property list identifier
!		plist2_id	- property list identifier
! Outputs:  
!               flag		- flag, possible values
!				  .TRUE. or .FALSE.
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1, flag is set to .FALSE.   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        September 30, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pequal_f(plist1_id, plist2_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pequal_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist1_id ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: plist2_id ! Property list identifier 
            LOGICAL, INTENT(OUT)       :: flag      ! Flag
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code
            INTEGER                    :: c_flag

            INTERFACE
              INTEGER FUNCTION h5pequal_c(plist1_id, plist2_id, c_flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PEQUAL_C'::h5pequal_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist1_id
              INTEGER(HID_T), INTENT(IN) :: plist2_id
              INTEGER, INTENT(OUT) :: c_flag
              END FUNCTION h5pequal_c
            END INTERFACE

            flag = .FALSE.
            hdferr = h5pequal_c(plist1_id, plist2_id, c_flag) 
            if (c_flag .GT. 0) flag = .TRUE.
          END SUBROUTINE h5pequal_f

!----------------------------------------------------------------------
! Name:		h5pset_buffer_f 
!
! Purpose: 	Sets sixe for conversion buffer
!
! Inputs:  
!		plist_id	- data transfer property list identifier
!               size		- buffer size 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 2, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_buffer_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_buffer_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! Data transfer property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: size  ! Buffer size in bytes; 
                                                   ! buffer is allocated and freed by 
                                                   ! the library.
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_buffer_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_BUFFER_C'::h5pset_buffer_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_buffer_c
            END INTERFACE

            hdferr = h5pset_buffer_c(plist_id, size) 
          END SUBROUTINE h5pset_buffer_f

!----------------------------------------------------------------------
! Name:		h5pget_buffer_f 
!
! Purpose: 	Gets size for conversion buffer
!
! Inputs:  
!		plist_id	- data transfer property list identifier
! Outputs:  
!               size		- buffer size 
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 2, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_buffer_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_buffer_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! Data transfer property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) :: size ! Buffer size in bytes; 
                                                   ! buffer is allocated and freed by 
                                                   ! the library.
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_buffer_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_BUFFER_C'::h5pget_buffer_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_buffer_c
            END INTERFACE

            hdferr = h5pget_buffer_c(plist_id, size) 
          END SUBROUTINE h5pget_buffer_f

!----------------------------------------------------------------------
! Name:		h5pfill_value_defined_f
!
! Purpose: 	Check if fill value is defined.
!
! Inputs:  
!		plist_id	- dataset creation property list identifier
! Outputs:  
!               flag            - fill value status flag
!                                 Possible values are:
!				    H5D_FILL_VALUE_ERROR_F
!				    H5D_FILL_VALUE_UNDEFINED_F
!				    H5D_FILL_VALUE_DEFAULT_F
!				    H5D_FILL_VALUE_USER_DEFINED_F
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 4, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pfill_value_defined_f(plist_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pfill_value_defined_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id 
            INTEGER, INTENT(OUT) :: flag
            INTEGER, INTENT(OUT)       :: hdferr    

            INTERFACE
              INTEGER FUNCTION h5pfill_value_defined_c(plist_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PFILL_VALUE_DEFINED_C'::h5pfill_value_defined_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER, INTENT(OUT) :: flag
              END FUNCTION h5pfill_value_defined_c
            END INTERFACE

            hdferr = h5pfill_value_defined_c(plist_id, flag) 
          END SUBROUTINE h5pfill_value_defined_f

!----------------------------------------------------------------------
! Name:		h5pset_alloc_time_f
!
! Purpose: 	Set space allocation time for dataset during creation.
!
! Inputs:  
!		plist_id	- dataset creation property list identifier
!               flag            - allocation time flag
!                                 Possible values are:
!				    H5D_ALLOC_TIME_ERROR_F
!				    H5D_ALLOC_TIME_DEFAULT_F
!				    H5D_ALLOC_TIME_EARLY_F
!				    H5D_ALLOC_TIME_LATE_F
!				    H5D_ALLOC_TIME_INCR_F
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 4, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_alloc_time_f(plist_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_alloc_time_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id 
            INTEGER, INTENT(IN) :: flag
            INTEGER, INTENT(OUT)       :: hdferr    

            INTERFACE
              INTEGER FUNCTION h5pset_alloc_time_c(plist_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_ALLOC_TIME_C'::h5pset_alloc_time_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER, INTENT(IN) :: flag
              END FUNCTION h5pset_alloc_time_c
            END INTERFACE

            hdferr = h5pset_alloc_time_c(plist_id, flag) 
          END SUBROUTINE h5pset_alloc_time_f

!----------------------------------------------------------------------
! Name:		h5pget_alloc_time_f
!
! Purpose: 	Get space allocation time for dataset during creation.
!
! Inputs:  
!		plist_id	- dataset creation property list identifier
! Outputs:  
!               flag            - allocation time flag
!                                 Possible values are:
!				    H5D_ALLOC_TIME_ERROR_F
!				    H5D_ALLOC_TIME_DEFAULT_F
!				    H5D_ALLOC_TIME_EARLY_F
!				    H5D_ALLOC_TIME_LATE_F
!				    H5D_ALLOC_TIME_INCR_F
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 4, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_alloc_time_f(plist_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_alloc_time_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id 
            INTEGER, INTENT(OUT) :: flag
            INTEGER, INTENT(OUT)       :: hdferr    

            INTERFACE
              INTEGER FUNCTION h5pget_alloc_time_c(plist_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_ALLOC_TIME_C'::h5pget_alloc_time_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER, INTENT(OUT) :: flag
              END FUNCTION h5pget_alloc_time_c
            END INTERFACE

            hdferr = h5pget_alloc_time_c(plist_id, flag) 
          END SUBROUTINE h5pget_alloc_time_f

!----------------------------------------------------------------------
! Name:		h5pset_fill_time_f
!
! Purpose: 	Set fill value writing time for dataset
!
! Inputs:  
!		plist_id	- dataset creation property list identifier
!               flag            - fill time flag
!                                 Possible values are:
!				    H5D_FILL_TIME_ERROR_F
!				    H5D_FILL_TIME_ALLOC_F
!				    H5D_FILL_TIME_NEVER_F
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 4, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_fill_time_f(plist_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_time_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id 
            INTEGER, INTENT(IN) :: flag
            INTEGER, INTENT(OUT)       :: hdferr    

            INTERFACE
              INTEGER FUNCTION h5pset_fill_time_c(plist_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FILL_TIME_C'::h5pset_fill_time_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER, INTENT(IN) :: flag
              END FUNCTION h5pset_fill_time_c
            END INTERFACE

            hdferr = h5pset_fill_time_c(plist_id, flag) 
          END SUBROUTINE h5pset_fill_time_f

!----------------------------------------------------------------------
! Name:		h5pget_fill_time_f
!
! Purpose: 	Get fill value writing time for dataset
!
! Inputs:  
!		plist_id	- dataset creation property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!               flag            - fill time flag
!                                 Possible values are:
!				    H5D_FILL_TIME_ERROR_F
!				    H5D_FILL_TIME_ALLOC_F
!				    H5D_FILL_TIME_NEVER_F
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 4, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_fill_time_f(plist_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fill_time_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id 
            INTEGER, INTENT(OUT) :: flag
            INTEGER, INTENT(OUT)       :: hdferr    

            INTERFACE
              INTEGER FUNCTION h5pget_fill_time_c(plist_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILL_TIME_C'::h5pget_fill_time_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER, INTENT(OUT) :: flag
              END FUNCTION h5pget_fill_time_c
            END INTERFACE

            hdferr = h5pget_fill_time_c(plist_id, flag) 
          END SUBROUTINE h5pget_fill_time_f

!----------------------------------------------------------------------
! Name:		h5pset_meta_block_size_f 
!
! Purpose: 	Sets the minimum size of metadata block allocations 
!
! Inputs:  
!		plist_id	- file access property list identifier
!               size		- metatdata block size
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_meta_block_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_meta_block_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: size  ! Block size in bytes; 
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_meta_block_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_META_BLOCK_SIZE_C'::h5pset_meta_block_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_meta_block_size_c
            END INTERFACE

            hdferr = h5pset_meta_block_size_c(plist_id, size) 
          END SUBROUTINE h5pset_meta_block_size_f

!----------------------------------------------------------------------
! Name:		h5pget_meta_block_size_f 
!
! Purpose: 	Gets the minimum size of metadata block allocations 
!
! Inputs:  
!		plist_id	- file access property list identifier
! Outputs:  
!               size		- metatdata block size
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_meta_block_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_meta_block_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) :: size  ! Block size in bytes; 
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_meta_block_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_META_BLOCK_SIZE_C'::h5pget_meta_block_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_meta_block_size_c
            END INTERFACE

            hdferr = h5pget_meta_block_size_c(plist_id, size) 
          END SUBROUTINE h5pget_meta_block_size_f

!----------------------------------------------------------------------
! Name:		h5pset_sieve_buf_size_f 
!
! Purpose: 	Sets the maximum size of the data sieve buffer
!
! Inputs:  
!		plist_id	- file access property list identifier
!               size		- sieve buffer size
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_sieve_buf_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_sieve_buf_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(SIZE_T), INTENT(IN) :: size  ! Buffer size in bytes; 
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_sieve_buf_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SIEVE_BUF_SIZE_C'::h5pset_sieve_buf_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(SIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_sieve_buf_size_c
            END INTERFACE

            hdferr = h5pset_sieve_buf_size_c(plist_id, size) 
          END SUBROUTINE h5pset_sieve_buf_size_f

!----------------------------------------------------------------------
! Name:		h5pget_sieve_buf_size_f 
!
! Purpose: 	Gets the maximum size of the data sieve buffer
!
! Inputs:  
!		plist_id	- file access property list identifier
! Outputs:  
!               size		- sieve buffer size
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_sieve_buf_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_sieve_buf_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(SIZE_T), INTENT(OUT) :: size   ! Buffer size in bytes 
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_sieve_buf_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_SIEVE_BUF_SIZE_C'::h5pget_sieve_buf_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(SIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_sieve_buf_size_c
            END INTERFACE

            hdferr = h5pget_sieve_buf_size_c(plist_id, size) 
          END SUBROUTINE h5pget_sieve_buf_size_f

!----------------------------------------------------------------------
! Name:		h5pset_small_data_block_size_f 
!
! Purpose: 	Sets the minimum size of "small" raw data block
!
! Inputs:  
!		plist_id	- file access property list identifier
!               size		- small raw data block size
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_small_data_block_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_small_data_block_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: size    ! Small raw data block size
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_small_data_block_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SMALL_DATA_BLOCK_SIZE_C'::h5pset_small_data_block_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_small_data_block_size_c
            END INTERFACE

            hdferr = h5pset_small_data_block_size_c(plist_id, size) 
          END SUBROUTINE h5pset_small_data_block_size_f

!----------------------------------------------------------------------
! Name:		h5pget_small_data_block_size_f 
!
! Purpose: 	Gets the minimum size of "small" raw data block
!
! Inputs:  
!		plist_id	- file access property list identifier
! Outputs:  
!               size		- small raw data block size
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_small_data_block_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_small_data_block_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! File access property list identifier 
            INTEGER(HSIZE_T), INTENT(OUT) :: size    ! Small raw data block size
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_small_data_block_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_SMALL_DATA_BLOCK_SIZE_C'::h5pget_small_data_block_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(HSIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_small_data_block_size_c
            END INTERFACE

            hdferr = h5pget_small_data_block_size_c(plist_id, size) 
          END SUBROUTINE h5pget_small_data_block_size_f

!----------------------------------------------------------------------
! Name:		h5pset_hyper_vector_size_f 
!
! Purpose: 	Set the number of "I/O" vectors (vector size)
!
! Inputs:  
!		plist_id	- dataset transfer property list identifier
!               size		- vector size
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_hyper_vector_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_hyper_vector_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! Dataset transfer property list identifier 
            INTEGER(SIZE_T), INTENT(IN) :: size     ! Vector size
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pset_hyper_vector_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PSET_HYPER_VECTOR_SIZE_C'::h5pset_hyper_vector_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(SIZE_T), INTENT(IN) :: size
              END FUNCTION h5pset_hyper_vector_size_c
            END INTERFACE

            hdferr = h5pset_hyper_vector_size_c(plist_id, size) 
          END SUBROUTINE h5pset_hyper_vector_size_f

!----------------------------------------------------------------------
! Name:		h5pget_hyper_vector_size_f 
!
! Purpose: 	Get the number of "I/O" vectors (vector size)
!
! Inputs:  
!		plist_id	- dataset transfer property list identifier
! Outputs:  
!               size		- vector size
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 7, 2002	
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_hyper_vector_size_f(plist_id, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_hyper_vector_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist_id ! Dataset transfer property list identifier 
            INTEGER(SIZE_T), INTENT(OUT) :: size     ! Vector size
            INTEGER, INTENT(OUT)       :: hdferr    ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_hyper_vector_size_c(plist_id, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
        !MS$ATTRIBUTES C,reference,alias:'_H5PGET_HYPER_VECTOR_SIZE_C'::h5pget_hyper_vector_size_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist_id
              INTEGER(SIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_hyper_vector_size_c
            END INTERFACE

            hdferr = h5pget_hyper_vector_size_c(plist_id, size) 
          END SUBROUTINE h5pget_hyper_vector_size_f

!----------------------------------------------------------------------
! Name:		h5pset_integer 
!
! Purpose: 	Sets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
!		value		- value to set property to
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_integer(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            INTEGER,   INTENT(IN) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pset_integer_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_INTEGER_C'::h5pset_integer_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER, INTENT(IN) :: value
              END FUNCTION h5pset_integer_c
            END INTERFACE

            name_len = LEN(name)
        hdferr = h5pset_integer_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pset_integer

!----------------------------------------------------------------------
! Name:		h5pset_real
!
! Purpose: 	Sets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
!		value		- value to set property to
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_real(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            REAL,   INTENT(IN) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pset_real_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_REAL_C'::h5pset_real_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              REAL, INTENT(IN) :: value
              END FUNCTION h5pset_real_c
            END INTERFACE

            name_len = LEN(name)
        hdferr = h5pset_real_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pset_real

!----------------------------------------------------------------------
! Name:		h5pset_double
!
! Purpose: 	Sets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
!		value		- value to set property to
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_double(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            DOUBLE PRECISION,   INTENT(IN) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pset_double_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_DOUBLE_C'::h5pset_double_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              DOUBLE PRECISION, INTENT(IN) :: value
              END FUNCTION h5pset_double_c
            END INTERFACE

            name_len = LEN(name)
        hdferr = h5pset_double_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pset_double

!----------------------------------------------------------------------
! Name:		h5pset_char
!
! Purpose: 	Sets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
!		value		- value to set property to
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pset_char(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            CHARACTER(LEN=*),   INTENT(IN) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len
            INTEGER :: value_len

            INTERFACE
              INTEGER FUNCTION h5psetc_c(prp_id, name, name_len, value, value_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSETC_C'::h5psetc_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              !DEC$ATTRIBUTES reference :: value
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              CHARACTER(LEN=*), INTENT(IN) :: value
              INTEGER, INTENT(IN)         :: value_len
              END FUNCTION h5psetc_c
            END INTERFACE

            name_len = LEN(name)
            value_len = LEN(value)
            hdferr = h5psetc_c(prp_id, name , name_len, value, value_len)
          END SUBROUTINE h5pset_char

!----------------------------------------------------------------------
! Name:		h5pget_integer 
!
! Purpose: 	Gets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
! Outputs:  
!		value		- value of property
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_integer(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            INTEGER,   INTENT(OUT) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pget_integer_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_INTEGER_C'::h5pget_integer_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER, INTENT(OUT) :: value
              END FUNCTION h5pget_integer_c
            END INTERFACE

            name_len = LEN(name)
        hdferr = h5pget_integer_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pget_integer

!----------------------------------------------------------------------
! Name:		h5pget_real
!
! Purpose: 	Gets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
! Outputs:  
!		value		- value of property
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_real(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            REAL,   INTENT(OUT) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pget_real_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_REAL_C'::h5pget_real_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              REAL, INTENT(OUT) :: value
              END FUNCTION h5pget_real_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pget_real_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pget_real

!----------------------------------------------------------------------
! Name:		h5pget_double
!
! Purpose: 	Gets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
! Outputs:  
!		value		- value of property
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_double(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            DOUBLE PRECISION,   INTENT(OUT) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pget_double_c(prp_id, name, name_len, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_DOUBLE_C'::h5pget_double_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              DOUBLE PRECISION, INTENT(OUT) :: value
              END FUNCTION h5pget_double_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pget_double_c(prp_id, name , name_len, value)
          END SUBROUTINE h5pget_double

!----------------------------------------------------------------------
! Name:		h5pget_char
!
! Purpose: 	Gets a property list value
!
! Inputs:  
!		prp_id		- iproperty list identifier to modify
!		name 		- name of property to modify
! Outputs:  
!		value		- value of property
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_char(prp_id, name, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            CHARACTER(LEN=*),   INTENT(OUT) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len
            INTEGER :: value_len

            INTERFACE
              INTEGER FUNCTION h5pgetc_c(prp_id, name, name_len, value, value_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGETC_C'::h5pgetc_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              !DEC$ATTRIBUTES reference :: value
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              CHARACTER(LEN=*), INTENT(OUT) :: value
              INTEGER, INTENT(IN)         :: value_len
              END FUNCTION h5pgetc_c
            END INTERFACE

            name_len = LEN(name)
            value_len = LEN(value)
            hdferr = h5pgetc_c(prp_id, name , name_len, value, value_len)
          END SUBROUTINE h5pget_char

!----------------------------------------------------------------------
! Name:		h5pexist_f 
!
! Purpose: 	Queries whether a property name exists in a property list or class. 
!
! Inputs:  
!		prp_id		- iproperty list identifier to query
!		name 		- name of property to check for
! Outputs:  
!               flag            - logical flag
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pexist_f(prp_id, name, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pexist_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to modify
            LOGICAL, INTENT(OUT) :: flag          ! .TRUE. if exists, .FALSE.
                                                  ! otherwise
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pexist_c(prp_id, name, name_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PEXIST_C'::h5pexist_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              END FUNCTION h5pexist_c
            END INTERFACE
            flag = .FALSE.
            name_len = LEN(name)
            hdferr = h5pexist_c(prp_id, name , name_len)
            if (hdferr > 0) then
                flag = .TRUE.
                hdferr = 0
            endif
          END SUBROUTINE h5pexist_f

!----------------------------------------------------------------------
! Name:		h5pget_size_f 
!
! Purpose: 	Queries the size of a property value in bytes.
!
! Inputs:  
!		prp_id		- property list identifier to query
!		name 		- name of property to query
! Outputs:  
!               size            - size of property in bytes
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_size_f(prp_id, name, size, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_size_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to query
            INTEGER(SIZE_T), INTENT(OUT) :: size  ! Size in bytes
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pget_size_c(prp_id, name, name_len, size)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_SIZE_C'::h5pget_size_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(OUT) :: size
              END FUNCTION h5pget_size_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5pget_size_c(prp_id, name , name_len, size)
          END SUBROUTINE h5pget_size_f

!----------------------------------------------------------------------
! Name:		h5pget_npros_f 
!
! Purpose: 	Queries number of properties in property list or class
!
! Inputs:  
!		prp_id		- iproperty list identifier to query
! Outputs:  
!               nprops          - number of properties in property object
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_nprops_f(prp_id, nprops, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_nprops_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id    ! Property list identifier 
            INTEGER(SIZE_T), INTENT(OUT) :: nprops  ! iNumber of properties
            INTEGER, INTENT(OUT) :: hdferr          ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_nprops_c(prp_id, nprops)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_NPROPS_C'::h5pget_nprops_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(SIZE_T), INTENT(OUT) :: nprops
              END FUNCTION h5pget_nprops_c
            END INTERFACE
            hdferr = h5pget_nprops_c(prp_id, nprops)
          END SUBROUTINE h5pget_nprops_f

!----------------------------------------------------------------------
! Name:		h5pget_class_name_f 
!
! Purpose: 	Queries the ithe name of a class.
!
! Inputs:  
!		prp_id		- property list identifier to query
! Outputs:  
!		name 		- name of a class
!		hdferr:		- error code		
!                                       
!				 	Success:  Actual lenght of the class name
!                                                 If provided buffer "name" is 
!                                                 smaller, than name will be 
!                                                 truncated to fit into
!                                                 provided user buffer
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_class_name_f(prp_id, name, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_class_name_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id  ! Property list identifier 
            CHARACTER(LEN=*), INTENT(INOUT) :: name  ! Buffer to retireve class name
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pget_class_name_c(prp_id, name, name_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_CLASS_NAME_C'::h5pget_class_name_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              CHARACTER(LEN=*), INTENT(INOUT) :: name
              INTEGER, INTENT(IN)         :: name_len
              END FUNCTION h5pget_class_name_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5pget_class_name_c(prp_id, name , name_len)
          END SUBROUTINE h5pget_class_name_f

!----------------------------------------------------------------------
! Name:		h5pget_class_parent_f 
!
! Purpose: 	Retrieves the parent class of a genric property class. 
!
! Inputs:  
!		prp_id		- property list identifier to query
! Outputs:  
!		parent_id 	- identifier of the parent class
!		hdferr:		- error code		
!                                       
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_class_parent_f(prp_id, parent_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_class_parent_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id     ! Property list identifier 
            INTEGER(HID_T), INTENT(OUT) :: parent_id ! Parent class property list 
                                                     ! identifier 
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

            INTERFACE
              INTEGER FUNCTION h5pget_class_parent_c(prp_id, parent_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_CLASS_PARENT_C'::h5pget_class_parent_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HID_T), INTENT(OUT) :: parent_id
              END FUNCTION h5pget_class_parent_c
            END INTERFACE
            hdferr = h5pget_class_parent_c(prp_id, parent_id)
          END SUBROUTINE h5pget_class_parent_f

!----------------------------------------------------------------------
! Name:		h5pisa_class_f 
!
! Purpose: 	Determines whether a property list is a member of a class. 
!
! Inputs:  
!		plist		- property list identifier 
!		pclass		- identifier of the property class
! Outputs:  
!               flag            - .TRUE. if a member, .FALSE. otherwise
!		hdferr:		- error code		
!                                       
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pisa_class_f(plist, pclass, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pisa_class_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist     ! Property list identifier 
            INTEGER(HID_T), INTENT(IN) :: pclass    ! Class identifier
            LOGICAL, INTENT(OUT) :: flag            ! logical flag
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

            INTERFACE
              INTEGER FUNCTION h5pisa_class_c(plist, pclass)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PISA_CLASS_C'::h5pisa_class_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: plist
              INTEGER(HID_T), INTENT(IN) :: pclass
              END FUNCTION h5pisa_class_c
            END INTERFACE
            flag = .FALSE.
            hdferr = h5pisa_class_c(plist, pclass)
            if (hdferr .gt. 0) then
                flag = .TRUE.
                hdferr = 0
            endif
          END SUBROUTINE h5pisa_class_f

!----------------------------------------------------------------------
! Name:		h5pcopy_prop_f 
!
! Purpose: 	Copies a property from one list or class to another.
!
! Inputs:  
!		dst_id		- Identifier of the destination property list
!		src_id		- Identifier of the source property list 
!		name 		- name of the property to copy
! Outputs:  
!		hdferr:		- error code		
!                                       
!				 	Success: 0 
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pcopy_prop_f(dst_id, src_id, name, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pcopy_prop_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: dst_id  ! Destination property list 
                                                  ! identifier 
            INTEGER(HID_T), INTENT(IN) :: src_id  ! Source property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name ! Property name
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pcopy_prop_c(dst_id, src_id, name, name_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCOPY_PROP_C'::h5pcopy_prop_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: dst_id
              INTEGER(HID_T), INTENT(IN) :: src_id
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              END FUNCTION h5pcopy_prop_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5pcopy_prop_c(dst_id, src_id, name , name_len)
          END SUBROUTINE h5pcopy_prop_f

!----------------------------------------------------------------------
! Name:		h5premove_f 
!
! Purpose: 	Removes a property from a property list. 

!
! Inputs:  
!		plid		- Property list identofoer
!		name 		- name of the property to remove
! Outputs:  
!		hdferr:		- error code		
!                                       
!				 	Success: 0 
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5premove_f(plid, name, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5premove_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plid   ! property list identifier
            CHARACTER(LEN=*), INTENT(IN) :: name ! name of property to remove
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5premove_c(plid, name, name_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREMOVE_C'::h5premove_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: plid 
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              END FUNCTION h5premove_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5premove_c(plid, name , name_len)
          END SUBROUTINE h5premove_f

!----------------------------------------------------------------------
! Name:		h5punregister_f 
!
! Purpose: 	Removes a property from a property list class. 

!
! Inputs:  
!		class		- Property list class identifier
!		name 		- name of the property to remove
! Outputs:  
!		hdferr:		- error code		
!                                       
!				 	Success: 0 
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5punregister_f(class, name, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5punregister_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class  ! property list class identifier
            CHARACTER(LEN=*), INTENT(IN) :: name ! name of property to remove
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5punregister_c(class, name, name_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PUNREGISTER_C'::h5punregister_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: class
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              END FUNCTION h5punregister_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5punregister_c(class, name , name_len)
          END SUBROUTINE h5punregister_f

!----------------------------------------------------------------------
! Name:		h5pclose_class_f 
!
! Purpose: 	Closes an existing property list class. 

!
! Inputs:  
!		class		- Property list class identifier
! Outputs:  
!		hdferr:		- error code		
!                                       
!				 	Success: 0 
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pclose_class_f(class, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pclose_class_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class  ! property list class identifier
            INTEGER, INTENT(OUT) :: hdferr  ! Error code

            INTERFACE
              INTEGER FUNCTION h5pclose_class_c(class)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCLOSE_CLASS_C'::h5pclose_class_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: class
              END FUNCTION h5pclose_class_c
            END INTERFACE
            hdferr = h5pclose_class_c(class)
          END SUBROUTINE h5pclose_class_f

!----------------------------------------------------------------------
! Name:		h5pcreate_class_f 
!
! Purpose: 	Create a new property list class

!
! Inputs:  
!		parent		- Property list identifier of the parent class
!                                 Possible values include:
!                                 H5P_ROOT_F
!                                 H5P_FILE_CREATE_F
!                                 H5P_FILE_ACCESS_F
!                                 H5P_DATASET_CREATE_F
!                                 H5P_DATASET_XFER_F
!                                 H5P_FILE_MOUNT_F
!		name 		- name of the class we are creating
! Outputs:  
!               class           - porperty list class identifier
!		hdferr:		- error code		
!                                       
!				 	Success: 0 
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 9, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pcreate_class_f(parent, name, class, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pcreate_class_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: parent  ! parent property list class 
                                                  ! identifier
            CHARACTER(LEN=*), INTENT(IN) :: name  ! name of property tocreate 
            INTEGER(HID_T), INTENT(OUT) :: class  ! property list class identifier
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pcreate_class_c(parent, name, name_len,&
                                                 class) 
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PCREATE_CLASS_C'::h5pcreate_class_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: parent
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(HID_T), INTENT(OUT) :: class
              END FUNCTION h5pcreate_class_c
            END INTERFACE
            name_len = LEN(name)
            hdferr = h5pcreate_class_c(parent, name , name_len, &
                                       class)
          END SUBROUTINE h5pcreate_class_f

!----------------------------------------------------------------------
! Name:		h5pregister_integer
!
! Purpose: 	Registers a permanent property with a property list class.
!
! Inputs:  
!		class		- property list class to register 
!                                 permanent property within
!		name 		- name of property to register
!               size            - size of property in bytes
!		value		- default value for property in newly 
!                                 created property lists
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pregister_integer(class, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pregister_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class   ! Property list class identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to register
            INTEGER(SIZE_T), INTENT(IN) :: size   ! Size of the property value	
            INTEGER,   INTENT(IN) :: value        ! Property value
            INTEGER, INTENT(OUT) :: hdferr        ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pregister_integer_c(class, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREGISTER_INTEGER_C'::h5pregister_integer_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: class
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              INTEGER, INTENT(IN) :: value
              END FUNCTION h5pregister_integer_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pregister_integer_c(class, name , name_len, size, value)
          END SUBROUTINE h5pregister_integer

!----------------------------------------------------------------------
! Name:		h5pregister_real
!
! Purpose: 	Registers a permanent property with a property list class.
!
! Inputs:  
!		class		- property list class to register 
!                                 permanent property within
!		name 		- name of property to register
!               size            - size of property in bytes
!		value		- default value for property in newly 
!                                 created property lists
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pregister_real(class, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pregister_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class   ! Property list class identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to register
            INTEGER(SIZE_T), INTENT(IN) :: size   ! size of the property value	
            REAL,   INTENT(IN) :: value           ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pregister_real_c(class, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREGISTER_REAL_C'::h5pregister_real_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: class
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              REAL, INTENT(IN) :: value
              END FUNCTION h5pregister_real_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pregister_real_c(class, name , name_len, size, value)
          END SUBROUTINE h5pregister_real

!----------------------------------------------------------------------
! Name:		h5pregister_double
!
! Purpose: 	Registers a permanent property with a property list class.
!
! Inputs:  
!		class		- property list class to register 
!                                 permanent property within
!		name 		- name of property to register
!               size            - size of property in bytes
!		value		- default value for property in newly 
!                                 created property lists
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pregister_double(class, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pregister_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class   ! Property list class identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to register
            INTEGER(SIZE_T), INTENT(IN) :: size  ! size of the property value	
            DOUBLE PRECISION,   INTENT(IN) :: value        ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pregister_double_c(class, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREGISTER_DOUBLE_C'::h5pregister_double_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: class
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              DOUBLE PRECISION, INTENT(IN) :: value
              END FUNCTION h5pregister_double_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pregister_double_c(class, name , name_len, size, value)
          END SUBROUTINE h5pregister_double

!----------------------------------------------------------------------
! Name:		h5pregister_char
!
! Purpose: 	Registers a permanent property with a property list class.
!
! Inputs:  
!		class		- property list class to register 
!                                 permanent property within
!		name 		- name of property to register
!               size            - size of property in bytes
!		value		- default value for property in newly 
!                                 created property lists
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pregister_char(class, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pregister_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: class   ! Property list class identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to register
            INTEGER(SIZE_T), INTENT(IN) :: size  ! size of the property value	
            CHARACTER(LEN=*),   INTENT(IN) :: value        ! Property value
            INTEGER, INTENT(OUT) :: hdferr  ! Error code
            INTEGER :: name_len
            INTEGER :: value_len

            INTERFACE
              INTEGER FUNCTION h5pregisterc_c(class, name, name_len, size, value, &
                                              value_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREGISTERC_C'::h5pregisterc_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              !DEC$ATTRIBUTES reference :: value
              INTEGER(HID_T), INTENT(IN) :: class
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              CHARACTER(LEN=*), INTENT(IN) :: value
              INTEGER, INTENT(IN)          :: value_len
              END FUNCTION h5pregisterc_c
            END INTERFACE

            name_len = LEN(name)
            value_len = LEN(value)
            hdferr = h5pregisterc_c(class, name , name_len, size, value, value_len)
          END SUBROUTINE h5pregister_char

!----------------------------------------------------------------------
! Name:		h5pinsert_integer
!
! Purpose: 	Registers a temporary property with a property list class.
!
! Inputs:  
!		plist		- property list identifier
!		name 		- name of property to insert
!               size            - size of property in bytes
!		value		- initial value for the property 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pinsert_integer(plist, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pinsert_integer
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist   ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to insert 
            INTEGER(SIZE_T), INTENT(IN) :: size   ! Size of the property value	
            INTEGER,   INTENT(IN) :: value        ! Property value
            INTEGER, INTENT(OUT) :: hdferr        ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pinsert_integer_c(plist, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PINSERT_INTEGER_C'::h5pinsert_integer_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: plist
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              INTEGER, INTENT(IN) :: value
              END FUNCTION h5pinsert_integer_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pinsert_integer_c(plist, name , name_len, size, value)
          END SUBROUTINE h5pinsert_integer

!----------------------------------------------------------------------
! Name:		h5pinsert_real
!
! Purpose: 	Registers a temporary property with a property list class.
!
! Inputs:  
!		plist		- property list identifier
!                                 permanent property within
!		name 		- name of property to insert
!               size            - size of property in bytes
!		value		- initial value for the property 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pinsert_real(plist, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pinsert_real
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist   ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to insert 
            INTEGER(SIZE_T), INTENT(IN) :: size   ! Size of the property value	
            REAL,   INTENT(IN) :: value           ! Property value
            INTEGER, INTENT(OUT) :: hdferr        ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pinsert_real_c(plist, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PINSERT_REAL_C'::h5pinsert_real_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: plist
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              REAL, INTENT(IN) :: value
              END FUNCTION h5pinsert_real_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pinsert_real_c(plist, name , name_len, size, value)
          END SUBROUTINE h5pinsert_real

!----------------------------------------------------------------------
! Name:		h5pinsert_double
!
! Purpose: 	Registers a temporary property with a property list class.
!
! Inputs:  
!		plist		- property list identifier
!                                 permanent property within
!		name 		- name of property to insert
!               size            - size of property in bytes
!		value		- initial value for the property 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pinsert_double(plist, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pinsert_double
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist   ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name  ! Name of property to insert 
            INTEGER(SIZE_T), INTENT(IN) :: size   ! Size of the property value	
            DOUBLE PRECISION, INTENT(IN) :: value ! Property value
            INTEGER, INTENT(OUT) :: hdferr        ! Error code
            INTEGER :: name_len

            INTERFACE
              INTEGER FUNCTION h5pinsert_double_c(plist, name, name_len, size, value)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PINSERT_DOUBLE_C'::h5pinsert_double_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: plist
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              DOUBLE PRECISION, INTENT(IN) :: value
              END FUNCTION h5pinsert_double_c
            END INTERFACE

            name_len = LEN(name)
            hdferr = h5pinsert_double_c(plist, name , name_len, size, value)
          END SUBROUTINE h5pinsert_double

!----------------------------------------------------------------------
! Name:		h5pinsert_char
!
! Purpose: 	Registers a temporary property with a property list class.
!
! Inputs:  
!		plist		- property list identifier
!                                 permanent property within
!		name 		- name of property to insert
!               size            - size of property in bytes
!		value		- initial value for the property 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!	        October 10, 2002	
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pinsert_char(plist, name, size, value, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pinsert_char
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: plist      ! Property list identifier 
            CHARACTER(LEN=*), INTENT(IN) :: name     ! Name of property to insert 
            INTEGER(SIZE_T), INTENT(IN) :: size      ! Size of property value	
            CHARACTER(LEN=*),   INTENT(IN) :: value  ! Property value
            INTEGER, INTENT(OUT) :: hdferr           ! Error code
            INTEGER :: name_len
            INTEGER :: value_len

            INTERFACE
              INTEGER FUNCTION h5pinsertc_c(plist, name, name_len, size, value, value_len)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PINSERTC_C'::h5pinsertc_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              !DEC$ATTRIBUTES reference :: value
              INTEGER(HID_T), INTENT(IN) :: plist
              CHARACTER(LEN=*), INTENT(IN) :: name
              INTEGER, INTENT(IN)         :: name_len
              INTEGER(SIZE_T), INTENT(IN) :: size 
              CHARACTER(LEN=*), INTENT(IN) :: value
              INTEGER, INTENT(IN)         :: value_len
              END FUNCTION h5pinsertc_c
            END INTERFACE

            name_len = LEN(name)
            value_len = LEN(value)
            hdferr = h5pinsertc_c(plist, name , name_len, size, value, value_len)
          END SUBROUTINE h5pinsert_char

!----------------------------------------------------------------------
! Name:		h5pset_shuffle_f 
!
! Purpose: 	Sets shuffling filter
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		March 12, 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_shuffle_f(prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_shuffle_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_shuffle_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_shuffle_c(prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SHUFFLE_C'::h5pset_shuffle_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              END FUNCTION h5pset_shuffle_c
            END INTERFACE
            hdferr = h5pset_shuffle_c(prp_id)

          END SUBROUTINE h5pset_shuffle_f

!----------------------------------------------------------------------
! Name:		h5pset_edc_check_f 
!
! Purpose: 	Enables/disables error detecting  
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
!               flag            - EDC flag; possible values:
!                                   H5Z_DISABLE_EDC_F
!                                   H5Z_ENABLE_EDC_F
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		March 13, 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_edc_check_f(prp_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_edc_check_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: flag          ! Checksum filter flag
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_edc_check_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_edc_check_c(prp_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_EDC_CHECK_C'::h5pset_edc_check_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: flag 
              END FUNCTION h5pset_edc_check_c
            END INTERFACE
            hdferr = h5pset_edc_check_c(prp_id, flag)

          END SUBROUTINE h5pset_edc_check_f

!----------------------------------------------------------------------
! Name:		h5pget_edc_check_f 
!
! Purpose: 	Queries error detecting  
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		March 13, 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_edc_check_f(prp_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_edc_check_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Dataset transfer property list identifier 
            INTEGER, INTENT(OUT) :: flag        ! Checksum filter flag
                                                 ! May have one of the following values:
                                                 ! H5Z_ERROR_EDC_F
                                                 ! H5Z_DISABLE_EDC_F
                                                 ! H5Z_ENABLE_EDC_F
                                                 ! H5Z_NO_EDC_F
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pget_edc_check_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_edc_check_c(prp_id, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_EDC_CHECK_C'::h5pget_edc_check_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(OUT) :: flag
              END FUNCTION h5pget_edc_check_c
            END INTERFACE
            hdferr = h5pget_edc_check_c(prp_id, flag)

          END SUBROUTINE h5pget_edc_check_f
!----------------------------------------------------------------------
! Name:		h5pset_fletcher32_f 
!
! Purpose: 	Sets Fletcher32 checksum of EDC for a dataset creation 
!               property list.
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		March 13, 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_fletcher32_f(prp_id, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fletcher32_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_fletcher32_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fletcher32_c(prp_id)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FLETCHER32_C'::h5pset_fletcher32_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              END FUNCTION h5pset_fletcher32_c
            END INTERFACE
            hdferr = h5pset_fletcher32_c(prp_id)

          END SUBROUTINE h5pset_fletcher32_f

!----------------------------------------------------------------------
! Name:		h5pset_family_offset_f 
!
! Purpose: 	Sets offset for family file driver.
!
! Inputs:  
!		prp_id		- file creation property list identifier
!               offset		- file offset
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		19 March 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_family_offset_f(prp_id, offset, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_family_offset_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER(HSIZE_T), INTENT(IN) :: offset ! Offset in bytes
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_family_offset_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_family_offset_c(prp_id, offset)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAMILY_OFFSET_C'::h5pset_family_offset_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER(HSIZE_T), INTENT(IN) :: offset 
              END FUNCTION h5pset_family_offset_c
            END INTERFACE
            hdferr = h5pset_family_offset_c(prp_id, offset)

          END SUBROUTINE h5pset_family_offset_f

!----------------------------------------------------------------------
! Name:		h5pset_fapl_multi_l 
!
! Purpose: 	Sets up use of the multi-file driver. 
!
! Inputs:  
!		prp_id		- file creation property list identifier
!               mem_map         - mapping array
!               memb_fapl       - property list for each memory usage type
!               memb_name       - names of member file
!               relax           - flag 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		20 March 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_fapl_multi_l(prp_id, memb_map, memb_fapl, memb_name, memb_addr, relax, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_multi_l
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
            INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_map
            INTEGER(HID_T), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_fapl
            CHARACTER(LEN=*), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_name
            !INTEGER(HADDR_T), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_addr
            REAL, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_addr
            LOGICAL, INTENT(IN) :: relax
            INTEGER, INTENT(OUT) :: hdferr       ! Error code
            INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1) :: lenm
            INTEGER :: maxlen
            INTEGER :: flag
            INTEGER :: i

!            INTEGER, EXTERNAL :: h5pset_fapl_multi_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_multi_c(prp_id, memb_map, memb_fapl, memb_name, lenm, &
                                                   maxlen, memb_addr, flag)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_MULTI_C'::h5pset_fapl_multi_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: memb_name
              INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
              INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_map
              INTEGER(HID_T), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_fapl
              CHARACTER(LEN=*), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_name
              REAL, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(IN) :: memb_addr
              !INTEGER(HADDR_T), DIMENSION(H5FD_MEM_NTYPES_F), INTENT(IN) :: memb_addr
              INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1) :: lenm
              INTEGER :: maxlen
              INTEGER, INTENT(IN) :: flag
              END FUNCTION h5pset_fapl_multi_c
            END INTERFACE
            maxlen = LEN(memb_name(1))
            do i=0, H5FD_MEM_NTYPES_F-1
             lenm(i) = LEN_TRIM(memb_name(i))
            enddo
            flag = 0
            if (relax) flag = 1
            hdferr = h5pset_fapl_multi_c(prp_id, memb_map, memb_fapl, memb_name, lenm, maxlen, memb_addr, flag) 

          END SUBROUTINE h5pset_fapl_multi_l
!----------------------------------------------------------------------
! Name:		h5pset_fapl_multi_s 
!
! Purpose: 	Sets up use of the multi-file driver. 
!
! Inputs:  
!		prp_id		- file creation property list identifier
!               relax           - flag 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		31 March 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_fapl_multi_s(prp_id, relax, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_fapl_multi_s
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
            LOGICAL, INTENT(IN) :: relax
            INTEGER, INTENT(OUT) :: hdferr       ! Error code
            INTEGER :: flag

!            INTEGER, EXTERNAL :: h5pset_fapl_multi_sc
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_fapl_multi_sc(prp_id,flag) 
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_FAPL_MULTI_SC'::h5pset_fapl_multi_sc
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
              INTEGER, INTENT(IN) :: flag
              END FUNCTION h5pset_fapl_multi_sc
            END INTERFACE
            flag = 0
            if (relax) flag = 1
            hdferr = h5pset_fapl_multi_sc(prp_id, flag) 

          END SUBROUTINE h5pset_fapl_multi_s
!----------------------------------------------------------------------
! Name:		h5pget_fapl_multi_f 
!
! Purpose: 	Sets up use of the multi-file driver. 
!
! Inputs:  
!		prp_id		- file creation property list identifier
! Outputs:  
!               mem_map         - mapping array
!               memb_fapl       - property list for each memory usage type
!               memb_name       - names of member file
!               relax           - flag 
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				maxlen_out - maximum length for memb_name array element 
!
! Programmer:	Elena Pourmal
!		24 March 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pget_fapl_multi_f(prp_id, memb_map, memb_fapl, memb_name, memb_addr, relax, hdferr, maxlen_out) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_fapl_multi_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
            INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(OUT) :: memb_map
            INTEGER(HID_T), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(OUT) :: memb_fapl
            CHARACTER(LEN=*), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(OUT) :: memb_name
            !INTEGER(HADDR_T), DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(OUT) :: memb_addr
            REAL, DIMENSION(0:H5FD_MEM_NTYPES_F-1), INTENT(OUT) :: memb_addr
            INTEGER, OPTIONAL, INTENT(OUT) :: maxlen_out 
            LOGICAL, INTENT(OUT) :: relax
            INTEGER, INTENT(OUT) :: hdferr       ! Error code
            INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1) :: lenm
            INTEGER :: maxlen
            INTEGER :: c_maxlen_out 
            INTEGER :: flag
            INTEGER :: i

!            INTEGER, EXTERNAL :: h5pget_fapl_multi_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_fapl_multi_c(prp_id, memb_map, memb_fapl, memb_name, lenm, &
                                                   maxlen, memb_addr, flag, c_maxlen_out)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FAPL_MULTI_C'::h5pget_fapl_multi_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: memb_name
              INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
              INTEGER, DIMENSION(H5FD_MEM_NTYPES_F), INTENT(OUT) :: memb_map
              INTEGER(HID_T), DIMENSION(H5FD_MEM_NTYPES_F), INTENT(OUT) :: memb_fapl
              CHARACTER(LEN=*), DIMENSION(H5FD_MEM_NTYPES_F), INTENT(OUT) :: memb_name
              REAL, DIMENSION(H5FD_MEM_NTYPES_F), INTENT(OUT) :: memb_addr
              INTEGER, DIMENSION(0:H5FD_MEM_NTYPES_F-1) :: lenm
              INTEGER :: maxlen
              INTEGER :: c_maxlen_out 
              INTEGER, INTENT(OUT) :: flag
              END FUNCTION h5pget_fapl_multi_c
            END INTERFACE
            maxlen = LEN(memb_name(0))
            do i=0, H5FD_MEM_NTYPES_F-1
             lenm(i) = LEN_TRIM(memb_name(i))
            enddo
            hdferr = h5pget_fapl_multi_c(prp_id, memb_map, memb_fapl, memb_name, lenm, maxlen, memb_addr, flag, c_maxlen_out) 
            relax = .TRUE.
            if(flag .eq. 0) relax = .FALSE. 
            if(present(maxlen_out)) maxlen_out = c_maxlen_out
          END SUBROUTINE h5pget_fapl_multi_f
!----------------------------------------------------------------------
! Name:		h5pset_szip_f 
!
! Purpose: 	Sets up use of szip compression
!
! Inputs:  
!		prp_id		- dataset creation property list identifier
!               options_mask
!               pixels_per_block - szip parameters
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		April 10 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pset_szip_f(prp_id, options_mask, pixels_per_block, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pset_szip_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Dataset creation property 
                                                 ! list identifier 
            INTEGER, INTENT(IN) :: options_mask
            INTEGER, INTENT(IN) :: pixels_per_block
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5pset_szip_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pset_szip_c(prp_id, options_mask, pixels_per_block) 
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PSET_SZIP_C'::h5pset_szip_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
              INTEGER, INTENT(IN) :: options_mask
              INTEGER, INTENT(IN) :: pixels_per_block
              END FUNCTION h5pset_szip_c
            END INTERFACE
            hdferr = h5pset_szip_c(prp_id, options_mask, pixels_per_block) 

          END SUBROUTINE h5pset_szip_f

!----------------------------------------------------------------------
! Name:		h5pall_filters_avail_f 
!
! Purpose: 	Checks if all filters set in the dataset creation
!               property list are available
!
! Inputs:  
!		prp_id		- data creation property list identifier
! Outputs:  
!               flag            - .TRUE. if all filters are available
!                                 .FALSE. otherwise
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		April 10 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------


          SUBROUTINE h5pall_filters_avail_f(prp_id, flag, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pall_filters_avail_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Dataset creation property 
                                                 ! list identifier 
            LOGICAL, INTENT(OUT) :: flag
            INTEGER, INTENT(OUT) :: hdferr       ! Error code
            INTEGER :: status

!            INTEGER, EXTERNAL :: h5pall_filters_avail_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pall_filters_avail_c(prp_id, status) 
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PALL_FILTERS_AVAIL_C'::h5pall_filters_avail_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id ! File creation property list identifier 
              INTEGER, INTENT(OUT) :: status
              END FUNCTION h5pall_filters_avail_c
            END INTERFACE
            flag = .TRUE.
            hdferr = h5pall_filters_avail_c(prp_id, status) 
            if (status .eq. 0 ) flag = .FALSE.

          END SUBROUTINE h5pall_filters_avail_f

!----------------------------------------------------------------------
! Name:		h5pget_filter_by_id_f 
!
! Purpose: 	Returns information about a filter in a pipeline
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
! Outputs:  
!		filter_id	- filter identifier
!		flags		- bit vector specifying certain general
!				  properties of the filter
!		cd_nelmts	- number of elements in cd_values
!		cd_values	- auxiliary data for the filter
!		namelen		- number of characters in the name buffer
!		name		- buffer to retrieve filter name
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		April 10 2003
!
! Modifications: 
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pget_filter_by_id_f(prp_id, filter_id, flags, cd_nelmts, cd_values, namelen, name, hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pget_filter_by_id_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 

            INTEGER, INTENT(IN) :: filter_id  ! Filter identifier
            INTEGER(SIZE_T), INTENT(INOUT) :: cd_nelmts  !Number of elements in cd_values.
            INTEGER, DIMENSION(*), INTENT(OUT) :: cd_values  !Auxiliary data for the filter.
            INTEGER, INTENT(OUT) :: flags  !Bit vector specifying certain general
                                          !properties of the filter.
            INTEGER(SIZE_T), INTENT(IN) :: namelen !Anticipated number of characters in name.
            CHARACTER(LEN=*), INTENT(OUT) :: name !Name of the filter

            INTEGER, INTENT(OUT) :: hdferr  ! Error code


!            INTEGER, EXTERNAL :: h5pget_filter_by_id_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pget_filter_by_id_c(prp_id, filter_id, flags, cd_nelmts,  &
                                              cd_values, namelen, name)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PGET_FILTER_BY_ID_C'::h5pget_filter_by_id_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: name
              INTEGER(HID_T), INTENT(IN) :: prp_id
              INTEGER, INTENT(IN) :: filter_id 
              INTEGER, DIMENSION(*), INTENT(OUT) :: cd_values
              INTEGER, INTENT(OUT) :: flags 
              INTEGER(SIZE_T), INTENT(INOUT) :: cd_nelmts
              INTEGER(SIZE_T), INTENT(IN) :: namelen
              CHARACTER(LEN=*), INTENT(OUT) :: name
              END FUNCTION h5pget_filter_by_id_c
            END INTERFACE

            hdferr = h5pget_filter_by_id_c(prp_id, filter_id, flags, cd_nelmts,  & 
                                     cd_values, namelen, name)
          END SUBROUTINE h5pget_filter_by_id_f

!----------------------------------------------------------------------
! Name:		h5pmodify_filter_f 
!
! Purpose: 	Adds a filter to the filter pipeline. 
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
!		filter		- filter to be modified
!		flags		- bit vector specifying certain general
!				  properties of the filter
!		cd_nelmts	- number of elements in cd_values
!		cd_values	- auxiliary data for the filter
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Elena Pourmal
!		April 10 2003
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5pmodify_filter_f(prp_id, filter, flags, cd_nelmts, cd_values,  hdferr) 
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5pmodify_filter_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Property list identifier 
            INTEGER, INTENT(IN) :: filter  !Filter to be modified
            INTEGER, INTENT(IN) :: flags  !Bit vector specifying certain general
                                          !properties of the filter.
            INTEGER(SIZE_T), INTENT(IN) :: cd_nelmts  !Number of elements in cd_values.
            INTEGER, DIMENSION(*), INTENT(IN) :: cd_values  !Auxiliary data for the filter.

            INTEGER, INTENT(OUT) :: hdferr  ! Error code

!            INTEGER, EXTERNAL :: h5pmodify_filter_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5pmodify_filter_c(prp_id, filter, flags, cd_nelmts, cd_values)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PMODIFY_FILTER_C'::h5pmodify_filter_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(IN) :: filter 
              INTEGER, INTENT(IN) :: flags 
              INTEGER(SIZE_T), INTENT(IN) :: cd_nelmts 
              INTEGER, DIMENSION(*), INTENT(IN) :: cd_values 
              END FUNCTION h5pmodify_filter_c
            END INTERFACE

            hdferr = h5pmodify_filter_c(prp_id, filter, flags, cd_nelmts, cd_values )
          END SUBROUTINE h5pmodify_filter_f

!----------------------------------------------------------------------
! Name:		h5premove_filter_f 
!
! Purpose: 	Delete one or more filters from the filter pipeline. 
!
! Inputs:  
!		prp_id		- data creation or transfer property list 
!				  identifier
!		filter		- filter to be removed
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
! Optional parameters:
!				NONE
!
! Programmer:	Quincey Koziol
!		January 27 2004
!
! Modifications: 	
!
! Comment:		
!----------------------------------------------------------------------

          SUBROUTINE h5premove_filter_f(prp_id, filter, hdferr)
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5premove_filter_f
!DEC$endif
!
            IMPLICIT NONE
            INTEGER(HID_T), INTENT(IN) :: prp_id ! Dataset creation property list
                                                 ! identifier
            INTEGER, INTENT(IN) :: filter        ! Filter to be removed
            INTEGER, INTENT(OUT) :: hdferr       ! Error code

!            INTEGER, EXTERNAL :: h5premove_filter_c
!  MS FORTRAN needs explicit interface for C functions called here.
!
            INTERFACE
              INTEGER FUNCTION h5premove_filter_c(prp_id, filter)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !MS$ATTRIBUTES C,reference,alias:'_H5PREMOVE_FILTER_C'::h5premove_filter_c
              !DEC$ ENDIF
              INTEGER(HID_T), INTENT(IN) :: prp_id 
              INTEGER, INTENT(IN) :: filter 
              END FUNCTION h5premove_filter_c
            END INTERFACE

            hdferr = h5premove_filter_c(prp_id, filter)
          END SUBROUTINE h5premove_filter_f

     END MODULE H5P
