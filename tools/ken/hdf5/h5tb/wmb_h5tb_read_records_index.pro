
;
; wmb_h5tb_read_records_index
; 
;              

pro wmb_h5tb_read_records_index, loc_id, $
                                 dset_name, $
                                 index, $
                                 databuffer
     
    compile_opt idl2, strictarrsubs
    
    
    n_reads = N_elements(index)
    
    mem_size = ulon64arr(1)
    elements = ulon64arr(1,n_reads)
    elements[0] = index
    
    
    ; get the number of records and fields in the table
    
    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, table_size
    
    
    ; make sure the indices are in bounds
    
    chk_min_error = min(index) lt 0
    chk_max_error = max(index) gt (table_size-1)
    
    if chk_min_error or chk_max_error then message, 'Invalid index value'
    
    
    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    

    ; get the dataspace handle
    
    sid = h5d_get_space(did)
    
    
    ; define the elements of the data to retrieve
    
    h5s_select_elements, sid, elements, /RESET
    
    
    ; create a memory dataspace handle
    
    mem_size[0] = n_reads
    m_sid = h5s_create_simple(mem_size)
    
    databuffer = h5d_read(did, MEMORY_SPACE=m_sid, FILE_SPACE=sid)
                                            
                            
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did    
    
end 
