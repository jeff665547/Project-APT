/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "H5TA.h"
#include <stdlib.h>

/*-------------------------------------------------------------------------
 * Table API example
 *
 * H5TBappend_records
 *
 *-------------------------------------------------------------------------
 */

#define NFIELDS  (hsize_t)     5
#define NRECORDS (hsize_t)     8
#define NRECORDS_ADD (hsize_t) 2
#define TABLE_NAME             "table"

int main( void )
{
 typedef struct Particle 
 {
  char   name[16];
  int    lati;
  int    longi;
  float  pressure;
  double temperature; 
 } Particle;

 Particle  dst_buf[NRECORDS+NRECORDS_ADD];

/* Define an array of Particles */
 Particle  p_data[NRECORDS] = { 
 {"zero",0,0, 0.0f, 0.0},
 {"one",10,10, 1.0f, 10.0},
 {"two",  20,20, 2.0f, 20.0},
 {"three",30,30, 3.0f, 30.0},
 {"four", 40,40, 4.0f, 40.0},
 {"five", 50,50, 5.0f, 50.0},
 {"six",  60,60, 6.0f, 60.0},
 {"seven",70,70, 7.0f, 70.0}
  };

 /* Calculate the size and the offsets of our struct members in memory */
 size_t dst_size =  sizeof( Particle );
 size_t dst_offset[NFIELDS] = { HOFFSET( Particle, name ),
                                HOFFSET( Particle, lati ),
                                HOFFSET( Particle, longi ),
                                HOFFSET( Particle, pressure ),
                                HOFFSET( Particle, temperature )};

 size_t dst_sizes[NFIELDS] = { sizeof( p_data[0].name),
                               sizeof( p_data[0].lati),
                               sizeof( p_data[0].longi),
                               sizeof( p_data[0].pressure),
                               sizeof( p_data[0].temperature)};
 
 /* Define field information */
 const char *field_names[NFIELDS] = 
 { "Name","Latitude", "Longitude", "Pressure", "Temperature" };
 hid_t      field_type[NFIELDS];
 hid_t      string_type;
 hid_t      file_id;
 hsize_t    chunk_size = 10;
 int        *fill_data = NULL;
 int        compress  = 0;
 herr_t     status; 
 int        i;

  /* Append particles */ 
 Particle particle_in[ NRECORDS_ADD ] = 
 {{ "eight",80,80, 8.0f, 80.0},
 {"nine",90,90, 9.0f, 90.0} };

 /* Initialize the field field_type */
 string_type = H5Tcopy( H5T_C_S1 );
 H5Tset_size( string_type, 16 );
 field_type[0] = string_type;
 field_type[1] = H5T_NATIVE_INT;
 field_type[2] = H5T_NATIVE_INT;
 field_type[3] = H5T_NATIVE_FLOAT;
 field_type[4] = H5T_NATIVE_DOUBLE;
   
 /* Create a new file using default properties. */
 file_id = H5Fcreate( "ex_table_02.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

 /* make a table */
 status=H5TBmake_table( "Table Title",file_id,TABLE_NAME,NFIELDS,NRECORDS, 
                        dst_size, field_names, dst_offset, field_type, 
                        chunk_size, fill_data, compress, p_data  );

 /* append two records */
 status=H5TBappend_records(file_id, TABLE_NAME,NRECORDS_ADD, dst_size, dst_offset, dst_sizes, 
  &particle_in );

 /* read the table */
 status=H5TBread_table( file_id, TABLE_NAME, dst_size, dst_offset, dst_sizes, dst_buf );

 /* print it by rows */
 for (i=0; i<NRECORDS+NRECORDS_ADD; i++) {
  printf ("%-5s %-5d %-5d %-5f %-5f", 
   dst_buf[i].name,
   dst_buf[i].lati,
   dst_buf[i].longi,
   dst_buf[i].pressure,
   dst_buf[i].temperature);
  printf ("\n");
 }
 
 /* Close the file. */
 H5Fclose( file_id );

 return 0;
}

