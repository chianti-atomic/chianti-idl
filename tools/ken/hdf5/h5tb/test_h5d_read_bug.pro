
pro test_h5d_read_bug

    compile_opt idl2, strictarrsubs

    ;-----------------------------------------------------------------
    ; create a new HDF5 file
    ;-----------------------------------------------------------------

    fn = dialog_pickfile(FILTER='*.h5', /WRITE)
    
    fn_info = file_info(fn)
    
    if fn_info.exists then begin

        message, 'Please create a new file'
    
    endif else begin
    
        ; create a new hdf5 file
        fid = h5f_create(fn)
        
    endelse

    ;-----------------------------------------------------------------
    ; define the compound data type structure
    ;-----------------------------------------------------------------
    
    field_names = ['testfield_int', $
                   'testfield_float', $
                   'testfield_string', $
                   'testfield_long', $
                   'testfield_double']
    
    field1_value = 0S
    field2_value = 0.0
    field3_value = ''
    field4_value = 0L
    field5_value = 0.0D
    
    record_definition = create_struct(field_names, $
                                      field1_value, $
                                      field2_value, $
                                      field3_value, $
                                      field4_value, $
                                      field5_value)

    ; set the maximum size of the text_string field to 10 characters, by 
    ; filling the 'testfield_string' field with a 10-character string
                 
    record_definition.testfield_string = '++++++++++'
    
    
    ;-----------------------------------------------------------------
    ; create a buffer of data to write to the table
    ;-----------------------------------------------------------------
    
    nrecords = 100
    
    databuffer = replicate(record_definition,nrecords)
    
    databuffer.testfield_int[0] = indgen(nrecords)
    databuffer.testfield_float[0] = findgen(nrecords) * 2
    databuffer.testfield_string[0] = 'HelloWorld'
    databuffer.testfield_long[0] = lindgen(nrecords) * 3
    databuffer.testfield_double[0] = dindgen(nrecords) - nrecords

    ;-----------------------------------------------------------------
    ; create the table
    ;-----------------------------------------------------------------
    
    loc_id = fid
    dset_name = 'test_hdf5_table'
    table_title = 'first table'
    chunk_size = nrecords / 2
    
    dims = ulon64arr(1)
    maxdims = ulon64arr(1)
    dims_chunk = ulon64arr(1)
    
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
    
    did = h5d_create(loc_id, dset_name, mem_type_id, sid, $
                     CHUNK_DIMENSIONS=dims_chunk)
        
        
    ; write the data
    h5d_write, did, databuffer
    
    ; terminate access to the data space
    h5s_close, sid

    ; end access to the dataset
    h5d_close, did
    
    ; release the datatype
    h5t_close, mem_type_id
    
    
    ;-----------------------------------------------------------------
    ; read back part of the table
    ;-----------------------------------------------------------------

    dset_name = 'test_hdf5_table'
    
    start_read_position = 0
    num_read_records = 10

    mem_dims = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)
    
    loc_id = fid
    did = h5d_open(loc_id, dset_name)
    
    ; choose the fields to read - create a read datatype id
    
    read_field_names = ['testfield_int', $
                        'testfield_string']

    read_field1_value = 0S
    read_field3_value = ''

    tmpstruct = create_struct(read_field_names, $
                              read_field1_value, $
                              read_field3_value)
                 
    tmpstruct.testfield_string = '++++++++++'
    
    read_type_id = h5t_idl_create(tmpstruct)
    
    
    ; get the file data space
    
    sid = h5d_get_space(did)
    
    
    ; define a hyperslab in the dataset from which to read the new data
    
    offset[0] = start_read_position
    count[0]  = num_read_records
    
    h5s_select_hyperslab, sid, offset, count, /RESET
    
    
    ; create a simple memory data space
    
    mem_dims[0] = num_read_records
    m_sid = h5s_create_simple(mem_dims)
    
    
    ; read the data
    
    databuffer = h5d_read(did, read_type_id, MEMORY_SPACE=m_sid, FILE_SPACE=sid)

    ; display the data
    
    print, databuffer
    
    ; close the dataspace, datatype, and dataset references - note that due 
    ; to a bug in h5d_read, h5d_read appears to close the read_type_id ?
    
    h5s_close, m_sid
    h5s_close, sid
    

    ; BUG!  The program crashes here, because read_type_id is no longer a 
    ; valid datatype ID after the call to h5d_read.  Comment out this line
    ; and the program executes without error.  
    
    ; Note that the datatype_id parameter should be added to h5d_write - 
    ; allowing parts of existing data tables to be selectively overwritten.
    h5t_close, read_type_id


    h5d_close, did
    
end