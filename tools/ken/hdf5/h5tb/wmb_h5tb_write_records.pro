
;
; wmb_h5tb_write_records
; 
; Purpose: Overwrites records. 
;
; Description: wmb_h5tb_write_records overwrites records starting at the 
;              zero index position start of the table named table_name 
;              attached to the object specified by the identifier loc_id. 
;
; Parameters:
;
;    loc_id
;        IN: Identifier of the file or group where the table is located. 
;    dset_name
;        IN: The name of the dataset to overwrite. 
;    start
;        IN: The zero based index record to start writing.  
;    nrecords
;        IN: The number of records to write. 
;    databuffer
;        IN: Buffer with data. 
;

pro wmb_h5tb_write_records, loc_id, $
                            dset_name, $
                            start, $
                            nrecords, $
                            databuffer

    compile_opt idl2, strictarrsubs

    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)
    
    ;-------------------------------------------------------------------------
    ; write the data buffer to disk
    ;-------------------------------------------------------------------------

    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    
    ; create a simple memory data space
    
    mem_dims[0] = nrecords
    m_sid = h5s_create_simple(mem_dims)
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    
    ; check that we are not overwriting the end of the dataset
    
    dims = h5s_get_simple_extent_dims(sid)
    
    if (start + nrecords) gt dims[0] then $
        message, 'Data buffer exceeds table size'
    
    ; define a hyperslab in the dataset to write the new data */
    
    offset[0] = start
    count[0]  = nrecords
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    h5d_write, did, databuffer, MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
    
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did

end
