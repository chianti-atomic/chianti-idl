
;
; wmb_h5tb_combine_tables
; 
; Purpose: Combines records from two tables into a third. 
;
; Description: wmb_h5tb_combine_tables combines records from two datasets 
;              named dset_name1 and dset_name2, to a new table named dset_name3.
;              These tables can be located on different files, identified by 
;              loc_id1 and loc_id2 (identifiers obtained with H5Fcreate). 
;              They can also be located on the same file. In this case one 
;              uses the same identifier for both parameters loc_id1 and 
;              loc_id2. If two files are used, the third table is written 
;              into the first file. 
;              
; Parameters:
; 
; loc_id1
;     IN: Identifier of the file or group in which the first table is located. 
; dset_name1
;     IN: The name of the first table to combine. 
; loc_id2
;     IN: Identifier of the file or group in which the second table is located. 
; dset_name2
;     IN: The name of the second table to combine. 
; dset_name3
;     IN: The name of the new table.
; new_table_title
;     IN: Title string for the new table.
; chunk_size
;     IN: The chunk size for the new table.
; compress
;     IN: Flag indicating whether the new table should be compressed.
;

pro wmb_h5tb_combine_tables, loc_id1, $
                             dset_name1, $
                             loc_id2, $
                             dset_name2, $
                             dset_name3, $
                             new_table_title, $
                             chunk_size, $
                             compress

    compile_opt idl2, strictarrsubs

    dims = ulon64arr(1)
    maxdims = ulon64arr(1)
    dims_chunk = ulon64arr(1)
    mem_size = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)

    maxdims[0] = -1
    dims_chunk[0] = chunk_size

    ;-------------------------------------------------------------------------
    ; get information about the two tables
    ;-------------------------------------------------------------------------
    
    wmb_h5tb_get_table_info, loc_id1, dset_name1, nfields1, nrecords1
    wmb_h5tb_get_table_info, loc_id2, dset_name2, nfields2, nrecords2
    
    wmb_h5tb_get_field_info, loc_id1, dset_name1, recdef1
    wmb_h5tb_get_field_info, loc_id2, dset_name2, recdef2
    
    
    ;-------------------------------------------------------------------------
    ; compare the record definitions to ensure that they match
    ;-------------------------------------------------------------------------

    recdef_matched = wmb_compare_struct(recdef1, $
                                        recdef2, $
                                        /compare_tag_names)
                                             
    if recdef_matched eq 0 then message, 'Unmatched record definitions'

    nfields = nfields1
    
    
    ;-------------------------------------------------------------------------
    ; open the first dataset
    ;-------------------------------------------------------------------------

    ; open the dataset
    
    did_1 = h5d_open(loc_id1, dset_name1)
    
    ; get the datatype
    
    tid_1 = h5d_get_type(did_1)
    
    ; get the file data space handle
    
    sid_1 = h5d_get_space(did_1)
    
    
    ;-------------------------------------------------------------------------
    ; make the merged table with no data originally
    ;-------------------------------------------------------------------------

    ; clone the type id
    
    tid_3 = h5t_copy(tid_1)
    
    ;-------------------------------------------------------------------------
    ; here we do not clone the file space from the 1st dataset, because we 
    ; want to create an empty table. Instead we create a new dataspace with 
    ; zero records and expandable.
    ;-------------------------------------------------------------------------

    dims[0] = 0
    
    sid_3 = h5s_create_simple(dims, max_dimensions = maxdims)
    
    if compress then begin
    
        compr_level = 6
        did_3 = h5d_create(loc_id1, dset_name3, tid_3, sid_3, $
                           CHUNK_DIMENSIONS = dims_chunk, GZIP = compr_level)
    
    endif else begin
    
        did_3 = h5d_create(loc_id1, dset_name3, tid_3, sid_3, $
                           CHUNK_DIMENSIONS = dims_chunk)
    
    endelse


    ;-------------------------------------------------------------------------
    ; attach the conforming table attributes
    ;-------------------------------------------------------------------------

    wmb_h5tb_attach_attributes, new_table_title, $
                                loc_id1, $
                                dset_name3, $
                                nfields, $
                                tid_3


    ;-------------------------------------------------------------------------
    ; read data from 1st table
    ;-------------------------------------------------------------------------
    
    ; define a hyperslab in the dataset of the size of the records
    
    offset[0] = 0
    count[0]  = nrecords1
    
    h5s_select_hyperslab, sid_1, offset, count, /RESET
    
    ; create a memory dataspace handle
    
    mem_size[0] = nrecords1
    m_sid = h5s_create_simple(mem_size)
    
    tmp_buf = h5d_read(did_1, MEMORY_SPACE = m_sid, FILE_SPACE = sid_1)
    
    
    ;-------------------------------------------------------------------------
    ; save data from 1st table into new table
    ;-------------------------------------------------------------------------

    wmb_h5tb_append_records, loc_id1, dset_name3, nrecords1, tmp_buf
    
    
    ;-------------------------------------------------------------------------
    ; release resources from 1st table
    ;-------------------------------------------------------------------------
    
    h5s_close, m_sid
    h5s_close, sid_1
    h5t_close, tid_1
    h5d_close, did_1


    ;-------------------------------------------------------------------------
    ; open the second dataset
    ;-------------------------------------------------------------------------

    ; open the dataset
    
    did_2 = h5d_open(loc_id2, dset_name2)
    
    
    ; get the file data space handle
    
    sid_2 = h5d_get_space(did_2)
    
    
    ;-------------------------------------------------------------------------
    ; read data from 2nd table
    ;-------------------------------------------------------------------------
    
    ; define a hyperslab in the dataset of the size of the records
    
    offset[0] = 0
    count[0]  = nrecords2
    
    h5s_select_hyperslab, sid_2, offset, count, /RESET
    
    ; create a memory dataspace handle
    
    mem_size[0] = nrecords2
    m_sid = h5s_create_simple(mem_size)
    
    tmp_buf = h5d_read(did_2, MEMORY_SPACE = m_sid, FILE_SPACE = sid_2)
    
    
    ;-------------------------------------------------------------------------
    ; save data from 2nd table into new table
    ;-------------------------------------------------------------------------

    wmb_h5tb_append_records, loc_id1, dset_name3, nrecords2, tmp_buf


    ;-------------------------------------------------------------------------
    ; release resources from 2st table
    ;-------------------------------------------------------------------------
    
    h5s_close, m_sid
    h5s_close, sid_2
    h5d_close, did_2


    ;-------------------------------------------------------------------------
    ; release resources from 3rd table
    ;-------------------------------------------------------------------------
    
    h5s_close, sid_3
    h5t_close, tid_3
    h5d_close, did_3
    
end

