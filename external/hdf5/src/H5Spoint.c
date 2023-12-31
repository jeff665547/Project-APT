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

/*
 * Programmer:  Quincey Koziol <koziol@ncsa.uiuc.edu>
 *              Tuesday, June 16, 1998
 *
 * Purpose:	Point selection data space I/O functions.
 */

#define H5S_PACKAGE		/*suppress error about including H5Spkg	  */


#include "H5private.h"		/* Generic Functions			  */
#include "H5Eprivate.h"		/* Error handling		  */
#include "H5FLprivate.h"	/* Free Lists	  */
#include "H5Iprivate.h"		/* ID Functions		  */
#include "H5MMprivate.h"	/* Memory Management functions		  */
#include "H5Spkg.h"		/* Dataspace functions			  */
#include "H5Vprivate.h"         /* Vector functions */

/* Static function prototypes */

/* Selection callbacks */
static herr_t H5S_point_copy(H5S_t *dst, const H5S_t *src, hbool_t share_selection);
static herr_t H5S_point_get_seq_list(const H5S_t *space, unsigned flags,
    H5S_sel_iter_t *iter, size_t maxseq, size_t maxbytes,
    size_t *nseq, size_t *nbytes, hsize_t *off, size_t *len);
static herr_t H5S_point_release(H5S_t *space);
static htri_t H5S_point_is_valid(const H5S_t *space);
static hssize_t H5S_point_serial_size(const H5S_t *space);
static herr_t H5S_point_serialize(const H5S_t *space, uint8_t *buf);
static herr_t H5S_point_deserialize(H5S_t *space, const uint8_t *buf);
static herr_t H5S_point_bounds(const H5S_t *space, hsize_t *start, hsize_t *end);
static htri_t H5S_point_is_contiguous(const H5S_t *space);
static htri_t H5S_point_is_single(const H5S_t *space);
static htri_t H5S_point_is_regular(const H5S_t *space);
static herr_t H5S_point_iter_init(H5S_sel_iter_t *iter, const H5S_t *space);

/* Selection iteration callbacks */
static herr_t H5S_point_iter_coords(const H5S_sel_iter_t *iter, hsize_t *coords);
static herr_t H5S_point_iter_block(const H5S_sel_iter_t *iter, hsize_t *start, hsize_t *end);
static hsize_t H5S_point_iter_nelmts(const H5S_sel_iter_t *iter);
static htri_t H5S_point_iter_has_next_block(const H5S_sel_iter_t *iter);
static herr_t H5S_point_iter_next(H5S_sel_iter_t *sel_iter, size_t nelem);
static herr_t H5S_point_iter_next_block(H5S_sel_iter_t *sel_iter);
static herr_t H5S_point_iter_release(H5S_sel_iter_t *sel_iter);

/* Selection properties for point selections */
const H5S_select_class_t H5S_sel_point[1] = {{
    H5S_SEL_POINTS,

    /* Methods on selection */
    H5S_point_copy,
    H5S_point_get_seq_list,
    H5S_point_release,
    H5S_point_is_valid,
    H5S_point_serial_size,
    H5S_point_serialize,
    H5S_point_deserialize,
    H5S_point_bounds,
    H5S_point_is_contiguous,
    H5S_point_is_single,
    H5S_point_is_regular,
    H5S_point_iter_init,
}};

/* Iteration properties for point selections */
static const H5S_sel_iter_class_t H5S_sel_iter_point[1] = {{
    H5S_SEL_POINTS,

    /* Methods on selection iterator */
    H5S_point_iter_coords,
    H5S_point_iter_block,
    H5S_point_iter_nelmts,
    H5S_point_iter_has_next_block,
    H5S_point_iter_next,
    H5S_point_iter_next_block,
    H5S_point_iter_release,
}};

/* Declare a free list to manage the H5S_pnt_node_t struct */
H5FL_DEFINE_STATIC(H5S_pnt_node_t);

/* Declare a free list to manage the H5S_pnt_list_t struct */
H5FL_DEFINE_STATIC(H5S_pnt_list_t);


/*-------------------------------------------------------------------------
 * Function:	H5S_point_iter_init
 *
 * Purpose:	Initializes iteration information for point selection.
 *
 * Return:	non-negative on success, negative on failure.
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, June 16, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5S_point_iter_init(H5S_sel_iter_t *iter, const H5S_t *space)
{
    FUNC_ENTER_NOAPI_NOFUNC(H5S_point_iter_init);

    /* Check args */
    assert (space && H5S_SEL_POINTS==H5S_GET_SELECT_TYPE(space));
    assert (iter);

    /* Initialize the number of points to iterate over */
    iter->elmt_left=space->select.num_elem;

    /* Start at the head of the list of points */
    iter->u.pnt.curr=space->select.sel_info.pnt_lst->head;

    /* Initialize type of selection iterator */
    iter->type=H5S_sel_iter_point;

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_init() */


/*-------------------------------------------------------------------------
 * Function:	H5S_point_iter_coords
 *
 * Purpose:	Retrieve the current coordinates of iterator for current
 *              selection
 *
 * Return:	non-negative on success, negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, April 22, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5S_point_iter_coords (const H5S_sel_iter_t *iter, hsize_t *coords)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_coords);

    /* Check args */
    assert (iter);
    assert (coords);

    /* Copy the offset of the current point */
    HDmemcpy(coords,iter->u.pnt.curr->pnt,sizeof(hsize_t)*iter->rank);

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_coords() */


/*-------------------------------------------------------------------------
 * Function:	H5S_point_iter_block
 *
 * Purpose:	Retrieve the current block of iterator for current
 *              selection
 *
 * Return:	non-negative on success, negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, June 2, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5S_point_iter_block (const H5S_sel_iter_t *iter, hsize_t *start, hsize_t *end)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_block);

    /* Check args */
    assert (iter);
    assert (start);
    assert (end);

    /* Copy the current point as a block */
    HDmemcpy(start,iter->u.pnt.curr->pnt,sizeof(hsize_t)*iter->rank);
    HDmemcpy(end,iter->u.pnt.curr->pnt,sizeof(hsize_t)*iter->rank);

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_block() */


/*-------------------------------------------------------------------------
 * Function:	H5S_point_iter_nelmts
 *
 * Purpose:	Return number of elements left to process in iterator
 *
 * Return:	non-negative number of elements on success, zero on failure
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, June 16, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static hsize_t
H5S_point_iter_nelmts (const H5S_sel_iter_t *iter)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_nelmts);

    /* Check args */
    assert (iter);

    FUNC_LEAVE_NOAPI(iter->elmt_left);
}   /* H5S_point_iter_nelmts() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_iter_has_next_block
 PURPOSE
    Check if there is another block left in the current iterator
 USAGE
    htri_t H5S_point_iter_has_next_block(iter)
        const H5S_sel_iter_t *iter;       IN: Pointer to selection iterator
 RETURNS
    Non-negative (TRUE/FALSE) on success/Negative on failure
 DESCRIPTION
    Check if there is another block available in the selection iterator.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t
H5S_point_iter_has_next_block(const H5S_sel_iter_t *iter)
{
    htri_t ret_value=TRUE;   /* Return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_has_next_block);

    /* Check args */
    assert (iter);

    /* Check if there is another point in the list */
    if(iter->u.pnt.curr->next==NULL)
        HGOTO_DONE(FALSE);

done:
    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_iter_has_next_block() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_iter_next
 PURPOSE
    Increment selection iterator
 USAGE
    herr_t H5S_point_iter_next(iter, nelem)
        H5S_sel_iter_t *iter;       IN: Pointer to selection iterator
        size_t nelem;               IN: Number of elements to advance by
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Advance selection iterator to the NELEM'th next element in the selection.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_iter_next(H5S_sel_iter_t *iter, size_t nelem)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_next);

    /* Check args */
    assert (iter);
    assert (nelem>0);

    /* Increment the iterator */
    while(nelem>0) {
        iter->u.pnt.curr=iter->u.pnt.curr->next;
        nelem--;
    } /* end while */

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_next() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_iter_next_block
 PURPOSE
    Increment selection iterator to next block
 USAGE
    herr_t H5S_point_iter_next_block(iter)
        H5S_sel_iter_t *iter;       IN: Pointer to selection iterator
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Advance selection iterator to the next block in the selection.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_iter_next_block(H5S_sel_iter_t *iter)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_next_block);

    /* Check args */
    assert (iter);

    /* Increment the iterator */
    iter->u.pnt.curr=iter->u.pnt.curr->next;

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_next_block() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_iter_release
 PURPOSE
    Release point selection iterator information for a dataspace
 USAGE
    herr_t H5S_point_iter_release(iter)
        H5S_sel_iter_t *iter;       IN: Pointer to selection iterator
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Releases all information for a dataspace point selection iterator
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_iter_release (H5S_sel_iter_t UNUSED * iter)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_iter_release);

    /* Check args */
    assert (iter);

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_iter_release() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_add
 PURPOSE
    Add a series of elements to a point selection
 USAGE
    herr_t H5S_point_add(space, num_elem, coord)
        H5S_t *space;           IN: Dataspace of selection to modify
        size_t num_elem;        IN: Number of elements in COORD array.
        const hsize_t *coord[];    IN: The location of each element selected
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    This function adds elements to the current point selection for a dataspace
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_add (H5S_t *space, H5S_seloper_t op, size_t num_elem, const hsize_t **_coord)
{
    H5S_pnt_node_t *top, *curr, *new_node; /* Point selection nodes */
    const hsize_t *coord=(const hsize_t *)_coord;     /* Pointer to the actual coordinates */
    unsigned i;                 /* Counter */
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_point_add);

    assert(space);
    assert(num_elem>0);
    assert(coord);
    assert(op==H5S_SELECT_SET || op==H5S_SELECT_APPEND || op==H5S_SELECT_PREPEND);

    top=curr=NULL;
    for(i=0; i<num_elem; i++) {
        /* Allocate space for the new node */
        if((new_node = H5FL_MALLOC(H5S_pnt_node_t))==NULL)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate point node");

        if((new_node->pnt = H5MM_malloc(space->extent.rank*sizeof(hsize_t)))==NULL)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate coordinate information");

        /* Copy over the coordinates */
        HDmemcpy(new_node->pnt,coord+(i*space->extent.rank),(space->extent.rank*sizeof(hsize_t)));

        /* Link into list */
        new_node->next=NULL;
        if(top==NULL)
            top=new_node;
        else
            curr->next=new_node;
        curr=new_node;
    } /* end for */

    /* Insert the list of points selected in the proper place */
    if(op==H5S_SELECT_SET || op==H5S_SELECT_PREPEND) {
        /* Append current list, if there is one */
        if(space->select.sel_info.pnt_lst->head!=NULL)
            curr->next=space->select.sel_info.pnt_lst->head;

        /* Put new list in point selection */
        space->select.sel_info.pnt_lst->head=top;
    } /* end if */
    else {  /* op==H5S_SELECT_APPEND */
        new_node=space->select.sel_info.pnt_lst->head;
        if(new_node!=NULL) {
            while(new_node->next!=NULL)
                new_node=new_node->next;

            /* Append new list to point selection */
            new_node->next=top;
        } /* end if */
        else
            space->select.sel_info.pnt_lst->head=top;
    } /* end else */

    /* Set the number of elements in the new selection */
    if(op==H5S_SELECT_SET)
        space->select.num_elem=num_elem;
    else
        space->select.num_elem+=num_elem;

done:
    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_add() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_release
 PURPOSE
    Release point selection information for a dataspace
 USAGE
    herr_t H5S_point_release(space)
        H5S_t *space;       IN: Pointer to dataspace
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Releases all point selection information for a dataspace
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_release (H5S_t *space)
{
    H5S_pnt_node_t *curr, *next;        /* Point selection nodes */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_release);

    /* Check args */
    assert (space);

    /* Delete all the nodes from the list */
    curr=space->select.sel_info.pnt_lst->head;
    while(curr!=NULL) {
        next=curr->next;
        H5MM_xfree(curr->pnt);
        H5FL_FREE(H5S_pnt_node_t,curr);
        curr=next;
    } /* end while */

    /* Free & reset the point list header */
    H5FL_FREE(H5S_pnt_list_t,space->select.sel_info.pnt_lst);
    space->select.sel_info.pnt_lst=NULL;

    /* Reset the number of elements in the selection */
    space->select.num_elem=0;

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_release() */


/*--------------------------------------------------------------------------
 NAME
    H5S_select_elements
 PURPOSE
    Specify a series of elements in the dataspace to select
 USAGE
    herr_t H5S_select_elements(dsid, op, num_elem, coord)
        hid_t dsid;             IN: Dataspace ID of selection to modify
        H5S_seloper_t op;       IN: Operation to perform on current selection
        size_t num_elem;        IN: Number of elements in COORD array.
        const hsize_t **coord; IN: The location of each element selected
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    This function selects array elements to be included in the selection for
    the dataspace.  The COORD array is a 2-D array of size <dataspace rank>
    by NUM_ELEM (ie. a list of coordinates in the dataspace).  The order of
    the element coordinates in the COORD array specifies the order that the
    array elements are iterated through when I/O is performed.  Duplicate
    coordinates are not checked for.  The selection operator, OP, determines
    how the new selection is to be combined with the existing selection for
    the dataspace.  Currently, only H5S_SELECT_SET is supported, which replaces
    the existing selection with the one defined in this call.  When operators
    other than H5S_SELECT_SET are used to combine a new selection with an
    existing selection, the selection ordering is reset to 'C' array ordering.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5S_select_elements (H5S_t *space, H5S_seloper_t op, size_t num_elem,
    const hsize_t **coord)
{
    herr_t ret_value=SUCCEED;  /* return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_select_elements);

    /* Check args */
    assert(space);
    assert(num_elem);
    assert(coord);
    assert(op==H5S_SELECT_SET || op==H5S_SELECT_APPEND || op==H5S_SELECT_PREPEND);

    /* If we are setting a new selection, remove current selection first */
    if(op==H5S_SELECT_SET || H5S_GET_SELECT_TYPE(space)!=H5S_SEL_POINTS) {
        if(H5S_SELECT_RELEASE(space)<0)
            HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL, "can't release point selection");
    } /* end if */

    /* Allocate space for the point selection information if necessary */
    if(H5S_GET_SELECT_TYPE(space)!=H5S_SEL_POINTS || space->select.sel_info.pnt_lst==NULL) {
        if((space->select.sel_info.pnt_lst = H5FL_CALLOC(H5S_pnt_list_t))==NULL)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate element information");
    } /* end if */

    /* Add points to selection */
    if(H5S_point_add(space,op,num_elem,coord)<0)
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTINSERT, FAIL, "can't insert elements");

    /* Set selection type */
    space->select.type=H5S_sel_point;

done:
    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_select_elements() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_copy
 PURPOSE
    Copy a selection from one dataspace to another
 USAGE
    herr_t H5S_point_copy(dst, src)
        H5S_t *dst;  OUT: Pointer to the destination dataspace
        H5S_t *src;  IN: Pointer to the source dataspace
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Copies all the point selection information from the source
    dataspace to the destination dataspace.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_copy(H5S_t *dst, const H5S_t *src, hbool_t UNUSED share_selection)
{
    H5S_pnt_node_t *curr, *new_node, *new_head;    /* Point information nodes */
    herr_t ret_value=SUCCEED;  /* return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_point_copy);

    assert(src);
    assert(dst);

    /* Allocate room for the head of the point list */
    if((dst->select.sel_info.pnt_lst=H5FL_MALLOC(H5S_pnt_list_t))==NULL)
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate point node");

    curr=src->select.sel_info.pnt_lst->head;
    new_head=NULL;
    while(curr!=NULL) {
        /* Create each point */
        if((new_node=H5FL_MALLOC(H5S_pnt_node_t))==NULL)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate point node");
        if((new_node->pnt = H5MM_malloc(src->extent.rank*sizeof(hsize_t)))==NULL)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate coordinate information");
        HDmemcpy(new_node->pnt,curr->pnt,(src->extent.rank*sizeof(hsize_t)));
        new_node->next=NULL;

        /* Keep the order the same when copying */
        if(new_head==NULL)
            new_head=dst->select.sel_info.pnt_lst->head=new_node;
        else {
            new_head->next=new_node;
            new_head=new_node;
        } /* end else */

        curr=curr->next;
    } /* end while */

done:
    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5S_point_copy() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_is_valid
 PURPOSE
    Check whether the selection fits within the extent, with the current
    offset defined.
 USAGE
    htri_t H5S_point_is_valid(space);
        H5S_t *space;             IN: Dataspace pointer to query
 RETURNS
    TRUE if the selection fits within the extent, FALSE if it does not and
        Negative on an error.
 DESCRIPTION
    Determines if the current selection at the current offet fits within the
    extent for the dataspace.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t
H5S_point_is_valid (const H5S_t *space)
{
    H5S_pnt_node_t *curr;      /* Point information nodes */
    unsigned u;                   /* Counter */
    htri_t ret_value=TRUE;     /* return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_is_valid);

    assert(space);

    /* Check each point to determine whether selection+offset is within extent */
    curr=space->select.sel_info.pnt_lst->head;
    while(curr!=NULL) {
        /* Check each dimension */
        for(u=0; u<space->extent.rank; u++) {
            /* Check if an offset has been defined */
            /* Bounds check the selected point + offset against the extent */
            if(((curr->pnt[u]+space->select.offset[u])>space->extent.size[u])
                    || (((hssize_t)curr->pnt[u]+space->select.offset[u])<0))
                HGOTO_DONE(FALSE)
        } /* end for */

        curr=curr->next;
    } /* end while */

done:
    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5S_point_is_valid() */


/*--------------------------------------------------------------------------
 NAME
    H5Sget_select_elem_npoints
 PURPOSE
    Get the number of points in current element selection
 USAGE
    hssize_t H5Sget_select_elem_npoints(dsid)
        hid_t dsid;             IN: Dataspace ID of selection to query
 RETURNS
    The number of element points in selection on success, negative on failure
 DESCRIPTION
    Returns the number of element points in current selection for dataspace.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
hssize_t
H5Sget_select_elem_npoints(hid_t spaceid)
{
    H5S_t	*space = NULL;      /* Dataspace to modify selection of */
    hssize_t ret_value;        /* return value */

    FUNC_ENTER_API(H5Sget_select_elem_npoints, FAIL);
    H5TRACE1("Hs","i",spaceid);

    /* Check args */
    if (NULL == (space=H5I_object_verify(spaceid, H5I_DATASPACE)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a data space");
    if(H5S_GET_SELECT_TYPE(space)!=H5S_SEL_POINTS)
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not an element selection");

    ret_value = H5S_GET_SELECT_NPOINTS(space);

done:
    FUNC_LEAVE_API(ret_value);
}   /* H5Sget_select_elem_npoints() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_serial_size
 PURPOSE
    Determine the number of bytes needed to store the serialized point selection
    information.
 USAGE
    hssize_t H5S_point_serial_size(space)
        H5S_t *space;             IN: Dataspace pointer to query
 RETURNS
    The number of bytes required on success, negative on an error.
 DESCRIPTION
    Determines the number of bytes required to serialize the current point
    selection information for storage on disk.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static hssize_t
H5S_point_serial_size (const H5S_t *space)
{
    H5S_pnt_node_t *curr;       /* Point information nodes */
    hssize_t ret_value;         /* return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_serial_size);

    assert(space);

    /* Basic number of bytes required to serialize point selection:
     *  <type (4 bytes)> + <version (4 bytes)> + <padding (4 bytes)> +
     *      <length (4 bytes)> + <rank (4 bytes)> + <# of points (4 bytes)> = 24 bytes
     */
    ret_value=24;

    /* Count points in selection */
    curr=space->select.sel_info.pnt_lst->head;
    while(curr!=NULL) {
        /* Add 4 bytes times the rank for each element selected */
        ret_value+=4*space->extent.rank;
        curr=curr->next;
    } /* end while */

    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5S_point_serial_size() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_serialize
 PURPOSE
    Serialize the current selection into a user-provided buffer.
 USAGE
    herr_t H5S_point_serialize(space, buf)
        H5S_t *space;           IN: Dataspace pointer of selection to serialize
        uint8 *buf;             OUT: Buffer to put serialized selection into
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Serializes the current element selection into a buffer.  (Primarily for
    storing on disk).
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_serialize (const H5S_t *space, uint8_t *buf)
{
    H5S_pnt_node_t *curr;   /* Point information nodes */
    uint8_t *lenp;          /* pointer to length location for later storage */
    uint32_t len=0;         /* number of bytes used */
    unsigned u;                /* local counting variable */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_serialize);

    assert(space);

    /* Store the preamble information */
    UINT32ENCODE(buf, (uint32_t)H5S_GET_SELECT_TYPE(space));  /* Store the type of selection */
    UINT32ENCODE(buf, (uint32_t)1);  /* Store the version number */
    UINT32ENCODE(buf, (uint32_t)0);  /* Store the un-used padding */
    lenp=buf;           /* keep the pointer to the length location for later */
    buf+=4;             /* skip over space for length */

    /* Encode number of dimensions */
    UINT32ENCODE(buf, (uint32_t)space->extent.rank);
    len+=4;

    /* Encode number of elements */
    UINT32ENCODE(buf, (uint32_t)space->select.num_elem);
    len+=4;

    /* Encode each point in selection */
    curr=space->select.sel_info.pnt_lst->head;
    while(curr!=NULL) {
        /* Add 4 bytes times the rank for each element selected */
        len+=4*space->extent.rank;

        /* Encode each point */
        for(u=0; u<space->extent.rank; u++)
            UINT32ENCODE(buf, (uint32_t)curr->pnt[u]);

        curr=curr->next;
    } /* end while */

    /* Encode length */
    UINT32ENCODE(lenp, (uint32_t)len);  /* Store the length of the extra information */

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_point_serialize() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_deserialize
 PURPOSE
    Deserialize the current selection from a user-provided buffer.
 USAGE
    herr_t H5S_point_deserialize(space, buf)
        H5S_t *space;           IN/OUT: Dataspace pointer to place selection into
        uint8 *buf;             IN: Buffer to retrieve serialized selection from
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Deserializes the current selection into a buffer.  (Primarily for retrieving
    from disk).
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_deserialize (H5S_t *space, const uint8_t *buf)
{
    H5S_seloper_t op=H5S_SELECT_SET;    /* Selection operation */
    uint32_t rank;           /* Rank of points */
    size_t num_elem=0;      /* Number of elements in selection */
    hsize_t *coord=NULL, *tcoord;   /* Pointer to array of elements */
    unsigned i,j;              /* local counting variables */
    herr_t ret_value;          /* return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_point_deserialize);

    /* Check args */
    assert(space);
    assert(buf);

    /* Deserialize points to select */
    buf+=16;    /* Skip over selection header */
    UINT32DECODE(buf,rank);  /* decode the rank of the point selection */
    if(rank!=space->extent.rank)
        HGOTO_ERROR(H5E_DATASPACE, H5E_BADRANGE, FAIL, "rank of pointer does not match dataspace");
    UINT32DECODE(buf,num_elem);  /* decode the number of points */

    /* Allocate space for the coordinates */
    if((coord = H5MM_malloc(num_elem*rank*sizeof(hsize_t)))==NULL)
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't allocate coordinate information");

    /* Retrieve the coordinates from the buffer */
    for(tcoord=coord,i=0; i<num_elem; i++)
        for(j=0; j<(unsigned)rank; j++,tcoord++)
            UINT32DECODE(buf, *tcoord);

    /* Select points */
    if((ret_value=H5S_select_elements(space,op,num_elem,(const hsize_t **)coord))<0)
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL, "can't change selection");

done:
    /* Free the coordinate array if necessary */
    if(coord!=NULL)
        H5MM_xfree(coord);

    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_deserialize() */


/*--------------------------------------------------------------------------
 NAME
    H5S_get_select_elem_pointlist
 PURPOSE
    Get the list of element points currently selected
 USAGE
    herr_t H5S_get_select_elem_pointlist(space, hsize_t *buf)
        H5S_t *space;           IN: Dataspace pointer of selection to query
        hsize_t startpoint;     IN: Element point to start with
        hsize_t numpoints;      IN: Number of element points to get
        hsize_t *buf;           OUT: List of element points selected
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
        Puts a list of the element points into the user's buffer.  The points
    start with the 'startpoint'th block in the list of points and put
    'numpoints' number of points into the user's buffer (or until the end of
    the list of points, whichever happen first)
        The point coordinates have the same dimensionality (rank) as the
    dataspace they are located within.  The list of points is formatted as
    follows: <coordinate> followed by the next coordinate, etc. until all the
    point information in the selection have been put into the user's buffer.
        The points are returned in the order they will be interated through
    when a selection is read/written from/to disk.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_get_select_elem_pointlist(H5S_t *space, hsize_t startpoint, hsize_t numpoints, hsize_t *buf)
{
    H5S_pnt_node_t *node;       /* Point node */
    int rank;                  /* Dataspace rank */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_get_select_elem_pointlist);

    assert(space);
    assert(buf);

    /* Get the dataspace extent rank */
    rank=space->extent.rank;

    /* Get the head of the point list */
    node=space->select.sel_info.pnt_lst->head;

    /* Iterate to the first point to return */
    while(node!=NULL && startpoint>0) {
        startpoint--;
        node=node->next;
      } /* end while */

    /* Iterate through the node, copying each hyperslab's information */
    while(node!=NULL && numpoints>0) {
        HDmemcpy(buf,node->pnt,sizeof(hsize_t)*rank);
        buf+=rank;
        numpoints--;
        node=node->next;
      } /* end while */

    FUNC_LEAVE_NOAPI(SUCCEED);
}   /* H5S_get_select_elem_pointlist() */


/*--------------------------------------------------------------------------
 NAME
    H5Sget_select_elem_pointlist
 PURPOSE
    Get the list of element points currently selected
 USAGE
    herr_t H5Sget_select_elem_pointlist(dsid, hsize_t *buf)
        hid_t dsid;             IN: Dataspace ID of selection to query
        hsize_t startpoint;     IN: Element point to start with
        hsize_t numpoints;      IN: Number of element points to get
        hsize_t *buf;           OUT: List of element points selected
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
        Puts a list of the element points into the user's buffer.  The points
    start with the 'startpoint'th block in the list of points and put
    'numpoints' number of points into the user's buffer (or until the end of
    the list of points, whichever happen first)
        The point coordinates have the same dimensionality (rank) as the
    dataspace they are located within.  The list of points is formatted as
    follows: <coordinate> followed by the next coordinate, etc. until all the
    point information in the selection have been put into the user's buffer.
        The points are returned in the order they will be interated through
    when a selection is read/written from/to disk.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5Sget_select_elem_pointlist(hid_t spaceid, hsize_t startpoint, hsize_t numpoints, hsize_t *buf)
{
    H5S_t	*space = NULL;      /* Dataspace to modify selection of */
    herr_t ret_value;        /* return value */

    FUNC_ENTER_API(H5Sget_select_elem_pointlist, FAIL);
    H5TRACE4("e","ihh*h",spaceid,startpoint,numpoints,buf);

    /* Check args */
    if(buf==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid pointer");
    if (NULL == (space=H5I_object_verify(spaceid, H5I_DATASPACE)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a data space");
    if(H5S_GET_SELECT_TYPE(space)!=H5S_SEL_POINTS)
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a point selection");

    ret_value = H5S_get_select_elem_pointlist(space,startpoint,numpoints,buf);

done:
    FUNC_LEAVE_API(ret_value);
}   /* H5Sget_select_elem_pointlist() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_bounds
 PURPOSE
    Gets the bounding box containing the selection.
 USAGE
    herr_t H5S_point_bounds(space, start, end)
        H5S_t *space;           IN: Dataspace pointer of selection to query
        hsize_t *start;         OUT: Starting coordinate of bounding box
        hsize_t *end;           OUT: Opposite coordinate of bounding box
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
    Retrieves the bounding box containing the current selection and places
    it into the user's buffers.  The start and end buffers must be large
    enough to hold the dataspace rank number of coordinates.  The bounding box
    exactly contains the selection, ie. if a 2-D element selection is currently
    defined with the following points: (4,5), (6,8) (10,7), the bounding box
    with be (4, 5), (10, 8).
        The bounding box calculations _does_ include the current offset of the
    selection within the dataspace extent.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_bounds(const H5S_t *space, hsize_t *start, hsize_t *end)
{
    H5S_pnt_node_t *node;       /* Point node */
    int rank;                   /* Dataspace rank */
    int i;                      /* index variable */
    herr_t ret_value=SUCCEED;   /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_point_bounds);

    assert(space);
    assert(start);
    assert(end);

    /* Get the dataspace extent rank */
    rank=space->extent.rank;

    /* Set the start and end arrays up */
    for(i=0; i<rank; i++) {
        start[i]=HSIZET_MAX;
        end[i]=0;
    } /* end for */

    /* Iterate through the node, checking the bounds on each element */
    node=space->select.sel_info.pnt_lst->head;
    while(node!=NULL) {
        for(i=0; i<rank; i++) {
            /* Check for offset moving selection negative */
            if(((hssize_t)node->pnt[i]+space->select.offset[i])<0)
                HGOTO_ERROR(H5E_DATASPACE, H5E_BADRANGE, FAIL, "offset moves selection out of bounds")

            if(start[i]>(node->pnt[i]+space->select.offset[i]))
                start[i]=node->pnt[i]+space->select.offset[i];
            if(end[i]<(node->pnt[i]+space->select.offset[i]))
                end[i]=node->pnt[i]+space->select.offset[i];
        } /* end for */
        node=node->next;
      } /* end while */

done:
    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_bounds() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_is_contiguous
 PURPOSE
    Check if a point selection is contiguous within the dataspace extent.
 USAGE
    htri_t H5S_point_is_contiguous(space)
        H5S_t *space;           IN: Dataspace pointer to check
 RETURNS
    TRUE/FALSE/FAIL
 DESCRIPTION
    Checks to see if the current selection in the dataspace is contiguous.
    This is primarily used for reading the entire selection in one swoop.
    This code currently doesn't properly check for contiguousness when there is
    more than one point, as that would take a lot of extra coding that we
    don't need now.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t
H5S_point_is_contiguous(const H5S_t *space)
{
    htri_t ret_value;  /* return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_is_contiguous);

    assert(space);

    /* One point is definitely contiguous */
    if(space->select.num_elem==1)
    	ret_value=TRUE;
    else	/* More than one point might be contiguous, but it's complex to check and we don't need it right now */
    	ret_value=FALSE;

    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_is_contiguous() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_is_single
 PURPOSE
    Check if a point selection is single within the dataspace extent.
 USAGE
    htri_t H5S_point_is_single(space)
        H5S_t *space;           IN: Dataspace pointer to check
 RETURNS
    TRUE/FALSE/FAIL
 DESCRIPTION
    Checks to see if the current selection in the dataspace is a single block.
    This is primarily used for reading the entire selection in one swoop.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t
H5S_point_is_single(const H5S_t *space)
{
    htri_t ret_value;  /* return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_is_single);

    assert(space);

    /* One point is definitely 'single' :-) */
    if(space->select.num_elem==1)
    	ret_value=TRUE;
    else
    	ret_value=FALSE;

    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_is_single() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_is_regular
 PURPOSE
    Check if a point selection is "regular"
 USAGE
    htri_t H5S_point_is_regular(space)
        const H5S_t *space;     IN: Dataspace pointer to check
 RETURNS
    TRUE/FALSE/FAIL
 DESCRIPTION
    Checks to see if the current selection in a dataspace is the a regular
    pattern.
    This is primarily used for reading the entire selection in one swoop.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    Doesn't check for points selected to be next to one another in a regular
    pattern yet.
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t
H5S_point_is_regular(const H5S_t *space)
{
    htri_t ret_value;  /* return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5S_point_is_regular);

    /* Check args */
    assert(space);

    /* Only simple check for regular points for now... */
    if(space->select.num_elem==1)
    	ret_value=TRUE;
    else
    	ret_value=FALSE;

    FUNC_LEAVE_NOAPI(ret_value);
}   /* H5S_point_is_regular() */


/*--------------------------------------------------------------------------
 NAME
    H5Sselect_elements
 PURPOSE
    Specify a series of elements in the dataspace to select
 USAGE
    herr_t H5Sselect_elements(dsid, op, num_elem, coord)
        hid_t dsid;             IN: Dataspace ID of selection to modify
        H5S_seloper_t op;       IN: Operation to perform on current selection
        size_t num_elem;        IN: Number of elements in COORD array.
        const hsize_t **coord; IN: The location of each element selected
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    This function selects array elements to be included in the selection for
    the dataspace.  The COORD array is a 2-D array of size <dataspace rank>
    by NUM_ELEM (ie. a list of coordinates in the dataspace).  The order of
    the element coordinates in the COORD array specifies the order that the
    array elements are iterated through when I/O is performed.  Duplicate
    coordinates are not checked for.  The selection operator, OP, determines
    how the new selection is to be combined with the existing selection for
    the dataspace.  Currently, only H5S_SELECT_SET is supported, which replaces
    the existing selection with the one defined in this call.  When operators
    other than H5S_SELECT_SET are used to combine a new selection with an
    existing selection, the selection ordering is reset to 'C' array ordering.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5Sselect_elements(hid_t spaceid, H5S_seloper_t op, size_t num_elem,
    const hsize_t **coord)
{
    H5S_t	*space = NULL;  /* Dataspace to modify selection of */
    herr_t ret_value;  /* return value */

    FUNC_ENTER_API(H5Sselect_elements, FAIL);
    H5TRACE4("e","iSsz**h",spaceid,op,num_elem,coord);

    /* Check args */
    if (NULL == (space=H5I_object_verify(spaceid, H5I_DATASPACE)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a dataspace");
    if (H5S_SCALAR==H5S_GET_EXTENT_TYPE(space))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "point doesn't support H5S_SCALAR space");
    if(coord==NULL || num_elem==0)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "elements not specified");
    if(!(op==H5S_SELECT_SET || op==H5S_SELECT_APPEND || op==H5S_SELECT_PREPEND))
        HGOTO_ERROR(H5E_ARGS, H5E_UNSUPPORTED, FAIL, "unsupported operation attempted");

    /* Call the real element selection routine */
    if((ret_value=H5S_select_elements(space,op,num_elem,coord))<0)
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL, "can't select elements");

done:
    FUNC_LEAVE_API(ret_value);
}   /* H5Sselect_elements() */


/*--------------------------------------------------------------------------
 NAME
    H5S_point_get_seq_list
 PURPOSE
    Create a list of offsets & lengths for a selection
 USAGE
    herr_t H5S_point_get_seq_list(space,flags,iter,maxseq,maxelem,nseq,nelem,off,len)
        H5S_t *space;           IN: Dataspace containing selection to use.
        unsigned flags;         IN: Flags for extra information about operation
        H5S_sel_iter_t *iter;   IN/OUT: Selection iterator describing last
                                    position of interest in selection.
        size_t maxseq;          IN: Maximum number of sequences to generate
        size_t maxelem;         IN: Maximum number of elements to include in the
                                    generated sequences
        size_t *nseq;           OUT: Actual number of sequences generated
        size_t *nelem;          OUT: Actual number of elements in sequences generated
        hsize_t *off;           OUT: Array of offsets
        size_t *len;            OUT: Array of lengths
 RETURNS
    Non-negative on success/Negative on failure.
 DESCRIPTION
    Use the selection in the dataspace to generate a list of byte offsets and
    lengths for the region(s) selected.  Start/Restart from the position in the
    ITER parameter.  The number of sequences generated is limited by the MAXSEQ
    parameter and the number of sequences actually generated is stored in the
    NSEQ parameter.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5S_point_get_seq_list(const H5S_t *space, unsigned flags, H5S_sel_iter_t *iter,
    size_t maxseq, size_t maxelem, size_t *nseq, size_t *nelem,
    hsize_t *off, size_t *len)
{
    size_t io_left;             /* The number of bytes left in the selection */
    size_t start_io_left;       /* The initial number of bytes left in the selection */
    H5S_pnt_node_t *node;       /* Point node */
    hsize_t dims[H5O_LAYOUT_NDIMS];     /* Total size of memory buf */
    int	ndims;                  /* Dimensionality of space*/
    hsize_t	acc;            /* Coordinate accumulator */
    hsize_t	loc;            /* Coordinate offset */
    size_t      curr_seq;       /* Current sequence being operated on */
    int         i;              /* Local index variable */
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER_NOAPI_NOINIT(H5S_point_get_seq_list);

    /* Check args */
    assert(space);
    assert(iter);
    assert(maxseq>0);
    assert(maxelem>0);
    assert(nseq);
    assert(nelem);
    assert(off);
    assert(len);

    /* Choose the minimum number of bytes to sequence through */
    H5_CHECK_OVERFLOW(iter->elmt_left,hsize_t,size_t);
    start_io_left=io_left=(size_t)MIN(iter->elmt_left,maxelem);

    /* Get the dataspace dimensions */
    if ((ndims=H5S_get_simple_extent_dims (space, dims, NULL))<0)
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTINIT, FAIL, "unable to retrieve data space dimensions");

    /* Walk through the points in the selection, starting at the current */
    /*  location in the iterator */
    node=iter->u.pnt.curr;
    curr_seq=0;
    while(node!=NULL) {
        /* Compute the offset of each selected point in the buffer */
        for(i=ndims-1,acc=iter->elmt_size,loc=0; i>=0; i--) {
            loc+=(node->pnt[i]+space->select.offset[i])*acc;
            acc*=dims[i];
        } /* end for */

        /* Check if this is a later point in the selection */
        if(curr_seq>0) {
            /* If a sorted sequence is requested, make certain we don't go backwards in the offset */
            if((flags&H5S_GET_SEQ_LIST_SORTED) && loc<off[curr_seq-1])
                break;

            /* Check if this point extends the previous sequence */
            /* (Unlikely, but possible) */
            if(loc==(off[curr_seq-1]+len[curr_seq-1])) {
                /* Extend the previous sequence */
                len[curr_seq-1]+=iter->elmt_size;
            } /* end if */
            else {
                /* Add a new sequence */
                off[curr_seq]=loc;
                len[curr_seq]=iter->elmt_size;

                /* Increment sequence count */
                curr_seq++;
            } /* end else */
        } /* end if */
        else {
            /* Add a new sequence */
            off[curr_seq]=loc;
            len[curr_seq]=iter->elmt_size;

            /* Increment sequence count */
            curr_seq++;
        } /* end else */

        /* Decrement number of elements left to process */
        io_left--;

        /* Move the iterator */
        iter->u.pnt.curr=node->next;
        iter->elmt_left--;

        /* Check if we're finished with all sequences */
        if(curr_seq==maxseq)
            break;

        /* Check if we're finished with all the elements available */
        if(io_left==0)
            break;

        /* Advance to the next point */
        node=node->next;
      } /* end while */

    /* Set the number of sequences generated */
    *nseq=curr_seq;

    /* Set the number of elements used */
    *nelem=start_io_left-io_left;

done:
    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5S_point_get_seq_list() */

