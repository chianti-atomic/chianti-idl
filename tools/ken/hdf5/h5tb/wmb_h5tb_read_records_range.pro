
;
; wmb_h5tb_read_records_range
; 
;              

pro wmb_h5tb_read_records_range,  loc_id, $
                                  dset_name, $
                                  firstrecord, $
                                  lastrecord, $
                                  input_stride, $
                                  databuffer
     
    compile_opt idl2, strictarrsubs
    
    
            
    count = ulon64arr(1)
    offset = ulon64arr(1)
    rangestride = ulon64arr(1)
    mem_size = ulon64arr(1)
        
    
    
    ; check the stride (it may be positive or negative)
    
    if input_stride eq 0 then message, 'Invalid stride value'
                            
    negative_stride_flag = (input_stride lt 0)

    abs_stride = ulong64(abs(input_stride))



    ; get the number of records and fields in the table

    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, table_size
    
    
    
    ; make sure the read request is in bounds
    
    if firstrecord lt 0 or firstrecord ge table_size then $
        message, 'Read request out of bounds'
    
    if lastrecord lt 0 or lastrecord ge table_size then $
        message, 'Read request out of bounds'

    if negative_stride_flag eq 0 then begin
    
        ; positive stride
        if lastrecord lt firstrecord then message, 'Invalid read range'
    
    endif else begin
        
        ; negative stride
        if lastrecord gt firstrecord then message, 'Invalid read range'
        
    endelse
    
    
    
    ; calculate the number of records to return
    
    nrecs_output = ceil((abs(firstrecord-lastrecord)+1)/double(abs_stride),/L64)
    
    
    
    ; if the stride is negative, we have to recalculate the starting point
    ; of the data read
    
    if negative_stride_flag eq 0 then begin
    
        ; positive stride
        disk_start_record = firstrecord
    
    endif else begin
        
        ; negative stride
        disk_start_record = firstrecord - ((nrecs_output-1)*abs_stride)
        
    endelse
    
    
    
    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    

    
    ; get the dataspace handle
    
    sid = h5d_get_space(did)
    
    
    
    ; define a hyperslab in the dataset of the size of the records
    
    offset[0] = disk_start_record
    count[0] = nrecs_output
    rangestride[0] = abs_stride
    
    h5s_select_hyperslab, sid, offset, count, STRIDE=rangestride, /RESET
    
    
    
    ; create a memory dataspace handle
    
    mem_size[0] = count[0]
    m_sid = h5s_create_simple(mem_size)
    
    databuffer = h5d_read(did, MEMORY_SPACE=m_sid, FILE_SPACE=sid)
                            
    
    
    ; if the stride is negative, reverse the output array
    
    if negative_stride_flag eq 1 then begin
        
        reversed_databuffer = databuffer[-1:0:-1]
        databuffer = temporary(reversed_databuffer)
        
    endif
                            
                            
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did    
    
end 
