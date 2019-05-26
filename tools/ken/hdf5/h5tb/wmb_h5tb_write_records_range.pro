
;
; wmb_h5tb_write_records_range
; 
; Purpose: Overwrites records. 
;
; Description: wmb_h5tb_write_records_range overwrites records starting at the 
;              zero index position start of the table named table_name 
;              attached to the object specified by the identifier loc_id. 
;
; Parameters:
;
;    loc_id
;        IN: Identifier of the file or group where the table is located. 
;    dset_name
;        IN: The name of the dataset to overwrite. 
;    firstrecord
;        IN: The zero based index record to start writing.  
;    lastrecord
;        IN: The zero based index record to end writing.  
;    stride
;        IN: The stride of the range to write.
;    databuffer
;        IN: Buffer with data. 
;

pro wmb_h5tb_write_records_range, loc_id, $
                                  dset_name, $
                                  firstrecord, $
                                  lastrecord, $
                                  input_stride, $
                                  databuffer

    compile_opt idl2, strictarrsubs

    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)
    rangestride = ulon64arr(1)    
    
    
    ; check the stride (it may be positive or negative)
    
    if input_stride eq 0 then message, 'Invalid stride value'
                            
    negative_stride_flag = (input_stride lt 0)

    abs_stride = ulong64(abs(input_stride))
    
    
    ; get the number of records and fields in the table

    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, table_size
    
    
    ; make sure the write request is in bounds
    
    if firstrecord lt 0 or firstrecord ge table_size then $
        message, 'Write request out of bounds'
    
    if lastrecord lt 0 or lastrecord ge table_size then $
        message, 'Write request out of bounds'

    if negative_stride_flag eq 0 then begin
    
        ; positive stride
        if lastrecord lt firstrecord then message, 'Invalid write range'
    
    endif else begin
        
        ; negative stride
        if lastrecord gt firstrecord then message, 'Invalid write range'
        
    endelse
    
    
    ; calculate the number of records to write
    
    nrecs_write = ceil((abs(firstrecord-lastrecord)+1)/double(abs_stride),/L64)
    
    
    ; check that the size of the databuffer matches the size of the range
    
    if N_elements(databuffer) ne nrecs_write then $
        message, 'Data buffer size does not match input range'
        
    
    ; if the stride is negative, we have to recalculate the starting point
    ; of the data write
    
    if negative_stride_flag eq 0 then begin
    
        ; positive stride
        disk_start_record = firstrecord
    
    endif else begin
        
        ; negative stride
        disk_start_record = firstrecord - ((nrecs_write-1)*abs_stride)
        
    endelse
    
    
    ; if the stride is negative, reverse the input array
    
    if negative_stride_flag eq 1 then begin
        
        reversed_databuffer = databuffer[-1:0:-1]
        
    endif


    ;-------------------------------------------------------------------------
    ; write the data buffer to disk
    ;-------------------------------------------------------------------------


    offset[0] = disk_start_record
    count[0] = nrecs_write
    rangestride[0] = abs_stride

    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    
    ; create a simple memory data space
    
    mem_dims[0] = nrecs_write
    m_sid = h5s_create_simple(mem_dims)
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    

    ; define a hyperslab in the dataset to write the new data

    h5s_select_hyperslab, sid, offset, count, STRIDE=rangestride, /RESET    
    
    
    ; write the data
    
    if negative_stride_flag eq 0 then begin
    
        h5d_write, did, databuffer, $
                   MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
                   
    endif else begin
        
        h5d_write, did, reversed_databuffer, $
                   MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
        
    endelse
    
    
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did

end
