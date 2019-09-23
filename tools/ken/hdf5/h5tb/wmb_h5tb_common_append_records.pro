
; wmb_h5tb_common_append_records
; 
; Purpose: Common code for reading records shared between H5PT and H5TB.
;


pro wmb_h5tb_common_append_records, dataset_id, $
                                    nrecords, $
                                    orig_table_size, $
                                    databuffer
    
    compile_opt idl2, strictarrsubs
    
    dims = ulon64arr(1)
    mem_dims = ulon64arr(1)
    offset =  ulon64arr(1)
    count = ulon64arr(1)
    
    ; extend the dataset
    
    dims[0] = nrecords + orig_table_size
    h5d_extend, dataset_id, dims
    
    
    ; create a simple memory data space
    
    mem_dims[0] =  nrecords
    m_sid = h5s_create_simple(mem_dims)
    
    
    ; get a copy of the new file data space for writing
    
    sid = h5d_get_space(dataset_id)
    
    
    ; define a hyperslab in the dataset
    
    offset[0] = orig_table_size
    count[0] = nrecords
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    
    ; write the records
    
    h5d_write, dataset_id, databuffer, MEMORY_SPACE_ID = m_sid, $
                                       FILE_SPACE_ID = sid
                                       

    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    
end


