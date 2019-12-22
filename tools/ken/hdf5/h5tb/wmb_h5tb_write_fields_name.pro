
;
; wmb_h5tb_write_fields_name
; 
; Purpose: Overwrites fields. 
;
; Description: wmb_h5tb_write_fields_name overwrites one or several fields 
;              contained in the buffer field_names from a dataset named 
;              table_name attached to the object specified by the 
;              identifier loc_id.    
;              
; Notes: It is not possible in this implementation of H5TB to selectively
;        overwrite individual fields within the table.  This is because
;        as of IDL 8.2 it is not possible to create a compound datatype 
;        where the field offsets are individually specified.  Hence, 
;        this procedure will load the specified data block from the table
;        (including all of the fields in the data block), modify the 
;        selected fields, and rewrite the entire block to the table.
;
; Parameters:
;
;    loc_id
;        IN: Identifier of the file or group where the table is located. 
;    dset_name
;        IN: The name of the dataset to overwrite. 
;    overwrite_field_names
;        IN: The indexes of the fields to write. 
;    start
;        IN: The zero based index record to start writing.  
;    nrecords
;        IN: The number of records to write. 
;    databuffer
;        IN: Buffer with data. 
;


pro wmb_h5tb_write_fields_name, loc_id, $
                                dset_name, $
                                overwrite_field_names, $
                                start, $
                                nrecords, $
                                databuffer

    compile_opt idl2, strictarrsubs

    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)

    ;-------------------------------------------------------------------------
    ; get the field names from the data table
    ;-------------------------------------------------------------------------
    
    wmb_h5tb_get_field_info, loc_id, dset_name, record_definition
    table_field_names = strupcase(tag_names(record_definition))

    ;-------------------------------------------------------------------------
    ; read the portion of the table to be overwritten
    ;-------------------------------------------------------------------------
    
    wmb_h5tb_read_records, loc_id, $
                           dset_name, $
                           start, $
                           nrecords, $
                           tmp_buf
                           
                           
    ;-------------------------------------------------------------------------
    ; modify the data buffer in memory
    ;-------------------------------------------------------------------------

    n_overwrite_fields = n_elements(overwrite_field_names)

    for i = 0, n_overwrite_fields-1 do begin
    
        tmpstr = strupcase(overwrite_field_names[i])
        fieldindex = where(table_field_names eq tmpstr, tmpcnt)
        if tmpcnt eq 0 then message, 'Field name not found'
        if tmpcnt gt 1 then message, 'Duplicate field name found'
        tmp_buf.(fieldindex)[0] = databuffer.(i)
    
    endfor


    ;-------------------------------------------------------------------------
    ; write the data buffer back to disk
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
    
    h5d_write, did, tmp_buf, MEMORY_SPACE_ID = m_sid, FILE_SPACE_ID = sid
    
    ; close
    
    h5s_close, m_sid
    h5s_close, sid
    h5d_close, did
    
end

