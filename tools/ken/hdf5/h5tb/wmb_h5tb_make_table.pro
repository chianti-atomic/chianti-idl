
; wmb_h5tb_make_table
;
; Purpose:      Creates and writes a table. 
;
; Description:  wmb_h5tb_make_table creates and writes a dataset named 
;               dset_name attached to the object specified by the 
;               identifier loc_id. 
;              
; Notes: This table creation function does not support fill values.
;
; Parameters:
; 
; table_title
;     IN: The title of the table.
; loc_id
;     IN: Identifier of the file or group to create the table within.
; dset_name
;     IN: The name of the dataset to create.
; record_definition
;     IN: A structure variable which defines the field names and data types 
;         of each record.
; nrecords
;     IN: The number of records in the table.
; chunk_size
;     IN: The chunk size.  
; compress
;     IN: Flag that turns compression on or off.
; databuffer 
;     IN: Data to be written to the table (optional).
;     


pro wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
                       

    compile_opt idl2, strictarrsubs

    dims = ulon64arr(1)
    maxdims = lon64arr(1)
    dims_chunk = ulon64arr(1)


    const_TABLE_CLASS = 'TABLE'
    const_TABLE_VER = 3.0
    const_HLTB_MAX_FIELD_LEN = 255
         
         
    dims[0] = nrecords
    dims_chunk[0] = chunk_size
    
    
    nfields = n_tags(record_definition)
    
    ; create the memory data type
    
    ; record_definition is a structure, complete with field names, which 
    ; fully defines the record data type
    
    mem_type_id = h5t_idl_create(record_definition)

    ; create a simple data space with unlimited size
    
    maxdims[0] = -1
    sid = h5s_create_simple(dims, max_dimensions = maxdims)
    
    ; create the dataset
    
    if compress eq 0 then begin
    
        did = h5d_create(loc_id, dset_name, mem_type_id, sid, $
                         CHUNK_DIMENSIONS=dims_chunk)
        
    endif else begin
    
        compr_level = 6
        did = h5d_create(loc_id, dset_name, mem_type_id, sid, $
                         CHUNK_DIMENSIONS=dims_chunk, GZIP=compr_level)
    
    endelse
    
    ; only write if there is something to write

    if N_elements(databuffer) ne 0 then begin
    
        h5d_write, did, databuffer

    endif


    ; terminate access to the data space
    h5s_close, sid

    ; end access to the dataset
    h5d_close, did

    ; attach the CLASS, VERSION, and TITLE attributes
    
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'CLASS', const_TABLE_CLASS
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'VERSION', const_TABLE_VER
    wmb_h5lt_set_attribute_string, loc_id, dset_name, 'TITLE', table_title
    
    
    ; attach the FIELD_ name attribute
    
    for i = 0, nfields-1 do begin
    
        member_name = h5t_get_member_name(mem_type_id, i)
        tmpstr = 'FIELD_' + strtrim(string(i),2) + '_NAME'
        wmb_h5lt_set_attribute_string, loc_id, dset_name, tmpstr, member_name

    endfor


    ; release the datatype
    
    h5t_close, mem_type_id
    
end
