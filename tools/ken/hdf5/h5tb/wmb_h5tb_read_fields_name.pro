
;
; wmb_h5tb_read_fields_name
; 
; Purpose: Reads one or several fields. The fields are identified by name. 
;
; Description: wmb_h5tb_read_fields_name reads the fields identified by
;              read_field_names from a dataset named table_name attached to the 
;              object specified by the identifier loc_id.
;
; Parameters:
;
;    loc_id
;        IN: Identifier of the file or group where the table is located. 
;    dset_name
;        IN: The name of the dataset to overwrite. 
;    read_field_names
;        IN: The indexes of the fields to write. 
;    start
;        IN: The zero based index record to start reading.  
;    nrecords
;        IN: The number of records to read. 
;    databuffer
;       OUT: Buffer with data. 
;


pro wmb_h5tb_read_fields_name, loc_id, $
                               dset_name, $
                               read_field_names, $
                               start, $
                               nrecords, $
                               databuffer

    compile_opt idl2, strictarrsubs

    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)

    ; get the field info

    wmb_h5tb_get_field_info, loc_id, dset_name, record_definition
    table_field_names = tag_names(record_definition)
    uc_field_names = strupcase(table_field_names)

    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)

    
    ; create a read datatype id

    n_read_fields = n_elements(read_field_names)
    
    for i = 0, n_read_fields-1 do begin

        tmpstr = strupcase(read_field_names[i])
        fieldindex = where(uc_field_names eq tmpstr, tmpcnt)
        if tmpcnt eq 0 then message, 'Field name not found'
        if tmpcnt gt 1 then message, 'Duplicate field name found'
        
        tmp_name = table_field_names[fieldindex]
        data_def = record_definition.(fieldindex)
    
        if i eq 0 then tmpstruct = create_struct(tmp_name,data_def) $
                  else tmpstruct = create_struct(tmpstruct, tmp_name, data_def)
    
        tmpstruct.(i) = data_def
    
    endfor
    
    read_type_id = h5t_idl_create(tmpstruct)


    ; get the file data space
    
    sid = h5d_get_space(did)
    
    
    ; define a hyperslab in the dataset from which to read the new data
    
    offset[0] = start
    count[0]  = nrecords
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    
    ; create a simple memory data space
    
    mem_dims[0] = nrecords
    m_sid = h5s_create_simple(mem_dims)
    
    
    ; read
    
    databuffer = h5d_read(did, read_type_id, MEMORY_SPACE=m_sid, FILE_SPACE=sid)

    
    ; close - note that due to a bug in h5d_read, h5d_read closes the 
    ; read_type_id and hence it should not need to be closed again
    
    h5s_close, m_sid
    h5s_close, sid
    if h5i_get_type(read_type_id) eq 'DATATYPE' then h5t_close, read_type_id
    h5d_close, did
    
end

