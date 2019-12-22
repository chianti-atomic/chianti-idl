
;
; wmb_h5tb_common_read_records
; 
; Purpose: Common code for reading records shared between H5PT and H5TB. 
;
;              

pro wmb_h5tb_common_read_records, dataset_id, $
                                  start, $
                                  nrecords, $
                                  table_size, $
                                  databuffer
     
    compile_opt idl2, strictarrsubs
                            
    count = ulon64arr(1)
    offset = ulon64arr(1)
    mem_size = ulon64arr(1)
    
    ; make sure the read request is in bounds
    
    if start + nrecords gt table_size then message, 'Read request out of bounds'
    
    ; get the dataspace handle
    
    sid = h5d_get_space(dataset_id)
    
    ; define a hyperslab in the dataset of the size of the records
    
    offset[0] = start
    count[0] = nrecords
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    ; create a memory dataspace handle
    
    mem_size[0] = count[0]
    m_sid = h5s_create_simple(mem_size)
    
    databuffer = h5d_read(dataset_id, MEMORY_SPACE=m_sid, FILE_SPACE=sid)
                            
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    
end 
