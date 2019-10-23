
;
; wmb_h5tb_write_records_index
; 
; Purpose: Overwrites records. 
;
; Description: wmb_h5tb_write_records_range overwrites records according 
;              to a zero-based index or index array 
;
; Parameters:
;
;    loc_id
;        IN: Identifier of the file or group where the table is located. 
;    dset_name
;        IN: The name of the dataset to overwrite. 
;    index
;        IN: The zero based index or index array defining the write locations
;    databuffer
;        IN: Buffer with data. 
;

pro wmb_h5tb_write_records_index, loc_id, $
                                  dset_name, $
                                  index, $
                                  databuffer

    compile_opt idl2, strictarrsubs

    n_writes = N_elements(index)

    mem_dims = ulon64arr(1)
    elements = ulon64arr(1,n_writes)
    elements[0] = index
    
    
    ; get the number of records and fields in the table
    
    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, table_size
    
    
    ; make sure the indices are in bounds
    
    chk_min_error = min(index) lt 0
    chk_max_error = max(index) gt (table_size-1)
    
    if chk_min_error or chk_max_error then message, 'Invalid index value'
    
    
    ; check that the size of the databuffer matches the size of the 
    ; index array
    
    if N_elements(databuffer) ne n_writes then $
        message, 'Data buffer size does not match index array'
        

    ;-------------------------------------------------------------------------
    ; write the data buffer to disk
    ;-------------------------------------------------------------------------


    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    
    ; create a simple memory data space
    
    mem_dims[0] = n_writes
    m_sid = h5s_create_simple(mem_dims)
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    

    ; define the elements of the data to retrieve
    
    h5s_select_elements, sid, elements, /RESET 
    
    
    ; write the data
    
    h5d_write, did, databuffer, MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
                   

    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did

end
