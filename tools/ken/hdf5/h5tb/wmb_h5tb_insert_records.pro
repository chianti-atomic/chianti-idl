
;
; wmb_h5tb_insert_records
; 
; Purpose: Insert records. 
;
; Description: wmb_h5tb_insert_records inserts records into the middle of the 
;              table ("pushing down" all the records after it).  
;              
; Parameters:
; 
; loc_id
;     IN: Identifier of the file or group in which the table is located. 
; dset_name
;     IN: The name of the dataset. 
; start
;     IN: The position to insert. 
; nrecords
;     IN: The number of records to insert. 
; databuffer
;     IN: Buffer with data. 
;

pro wmb_h5tb_insert_records, loc_id, $
                             dset_name, $
                             start, $
                             nrecords, $
                             databuffer

    compile_opt idl2, strictarrsubs

    dims = ulon64arr(1)
    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)

    ;-------------------------------------------------------------------------
    ; read the records after the inserted one(s)
    ;-------------------------------------------------------------------------

    ; get the dimensions
    
    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, ntotal_records
    
    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    

    ; read the records after the inserted one(s)
    
    read_nrecords = ntotal_records - start
    
    wmb_h5tb_read_records, loc_id, $
                           dset_name, $
                           start, $
                           read_nrecords, $
                           tmp_buf

    ; extend the dataset
    
    dims[0] = ntotal_records + nrecords
    h5d_extend, did, dims
    
    
    ;-------------------------------------------------------------------------
    ; write the inserted records
    ;-------------------------------------------------------------------------

    ; create a simple memory data space
    
    mem_dims[0] = nrecords
    m_sid = h5s_create_simple(mem_dims)
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    
    ; define a hyperslab in the dataset to write the new data */
    
    offset[0] = start
    count[0]  = nrecords
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    h5d_write, did, databuffer, MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
    
    ; terminate access to the dataspace
    
    h5s_close, m_sid
    h5s_close, sid
    
    
    ;-------------------------------------------------------------------------
    ; write the "pushed down" records
    ;-------------------------------------------------------------------------
    
    ; create a simple memory data space
    
    mem_dims[0] = read_nrecords
    m_sid = h5s_create_simple(mem_dims)
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    
    ; define a hyperslab in the dataset to write the new data
    
    offset[0] = start + nrecords
    count[0]  = read_nrecords
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    h5d_write, did, tmp_buf, MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
    
    
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did
    
end

