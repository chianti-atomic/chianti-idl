;-----------------------------------------------------------------------------
; Title: wmb_lib HDF5 libraray
;
; Purpose: Reading/writing/modification of HDF5 tables in IDL.
;
; Description:  The wmb_hdf5* library implements the basic functions for 
;               reading and writing HDF5 tables from IDL.  This library was
;               created by porting the C code from the HDF5 distribution
;               (H5TB.c, from HDF5 version 1.8.10) and modifying it to work 
;               with the HDF5 API currently supported by IDL (as of IDL 
;               version 8.2, most of the HDF5 version 1.6 API is supported).
;           
;       The C source code which this library is based on can be found
;       at the following location:
;       
;       http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/hl/src/H5TB.c
;
;       By translating this code into IDL, the following functions were ported:
;           
;       Table creation:
;       
;           wmb_h5tb_make_table.pro
;
;       Storage:
;       
;           wmb_h5tb_append_records.pro
;           wmb_h5tb_write_records.pro
;           wmb_h5tb_write_fields_index.pro
;           wmb_h5tb_write_fields_name.pro
;           
;       Retrieval: 
;       
;           wmb_h5tb_read_table.pro
;           wmb_h5tb_read_records.pro
;           wmb_h5tb_read_fields_index.pro
;           wmb_h5tb_read_fields_name.pro
;           
;       Query:
;       
;           wmb_h5tb_get_table_info.pro
;           wmb_h5tb_get_field_info.pro
;           
;       Modification:
;       
;           wmb_h5tb_insert_records.pro
;           wmb_h5tb_add_records_from.pro
;           wmb_h5tb_combine_tables.pro
;
; Limitations:
;
;       Several of the H5TB* functiions from the could not
;       be ported, due to limitations of the current version of the 
;       IDL HDF5 implementation.  These functions are listed below, 
;       along with the reasons why they could not be ported.
;       
;       H5TB_delete_record:  Currently, IDL does not implement the
;                            H5D_set_extent function, and therefore
;                            there is no way to reduce the size of
;                            an existing dataset.
;
;       H5TB_insert_field:   The version of H5D_write implemented by
;       and                  IDL does not allow the specification of 
;       H5TB_delete_field    a "write datatype id", which is required
;                            to enable the selective writing of a subset
;                            of the fields within a table.  Also, IDL does 
;                            not allow the chunk dimensions and the 
;                            compression state of an existing dataset 
;                            to be queried.
;           
;       Finally, this library does not support fill values, since there
;       is no way of specifying a fill value through the current HDF5
;       API implemented by IDL.
;       
;       
; Dependencies:
; 
;       The wmb_h5tb* library depends on the following procedures from 
;       the wmb_lib:
;       
;           wmb_h5lt_find_attribute.pro
;           wmb_h5lt_get_attribute.pro
;           wmb_h5lt_set_attribute.pro
;           wmb_h5o_open.pro
;           wmb_h5o_close.pro
;           wmb_compare_struct.pro
;           
;           
; Defining the structure of an HDF5 table:
;                                
;       An HDF5 table can be thought of as a one-dimensional array of 
;       structure variables.  The length of the table can be extended 
;       up to an unlimited dimension.  
;        
;       Before creating an HDF5 table, its data structure must be
;       defined.  This is done by creating a structure-type variable 
;       which defines the datatype of each row of the table.
;        
;       The field names of this stucture variable are used to set the 
;       field names of the table.  The data values of the structure 
;       define the datatype of each table field.  Here we refer to this
;       structure as the "record_definition" structure.
;        
;       String datatypes require special consideration.  At present, only
;       fixed-length string types are supported in the wmb_h5tb* library.
;       Therefore, if strings are to be stored in the table, the maximum
;       string length for each field must be set when the table is created.  
;       This is done by writing string data into the record definition 
;       structure.  For each string field, a string value with a number
;       of characters equal to the desired fixed length of the field is 
;       assigned to the record definition structure. 
;
;
; Browsing HDF5 table data:
;       
;       The built in IDL HDF5 browser (H5_BROWSER) is not capable of 
;       displaying HDF5 tables.  Therefore, to view and edit the HDF5 
;       tables users are encouraged to download HDFView, a free 
;       Java-based HDF5 browser which is fully compatible with all 
;       HDF5 files.  HDFView is available at the address below.
;       
;       http://www.hdfgroup.org/hdf-java-html/hdfview/
;       
; Version: 1.0, 14 May 2013, WMB
;-----------------------------------------------------------------------------


;-----------------------------------------------------------------------------
;
;   Make a new HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_make_table_example

    compile_opt idl2, strictarrsubs

    ; choose an output file name

    fn = dialog_pickfile(FILTER='*.h5', /WRITE)
    
    fn_info = file_info(fn)
    
    if fn_info.exists then begin
        chk_file = file_test(fn, /WRITE)
        if ~chk_file then message, 'Error opening file'
    endif
    
    if fn_info.exists then begin

        ; open the hdf5 file
        fid = h5f_open(fn, /WRITE)
    
    endif else begin
    
        ; create a new hdf5 file
        fid = h5f_create(fn)
        
    endelse
    
    ;-----------------------------------------------------------------
    ; Define the table structure
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
    ; Create a buffer of data to write to the table
    ;-----------------------------------------------------------------
    
    nrecords = 100543
    
    databuffer = replicate(record_definition,nrecords)
    
    databuffer.testfield_int[0] = indgen(nrecords)
    databuffer.testfield_float[0] = findgen(nrecords) * 2
    databuffer.testfield_string[0] = 'HelloWorld'
    databuffer.testfield_long[0] = lindgen(nrecords) * 3
    databuffer.testfield_double[0] = dindgen(nrecords) - nrecords
    
    ;-----------------------------------------------------------------
    ; Create the table
    ;-----------------------------------------------------------------
    
    loc_id = fid
    dset_name = 'test_hdf5_table'
    table_title = 'first table'
    chunk_size = nrecords / 2
    compress = 0
    
    wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
    
    ; close the file

    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Read from an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_table_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ)
    if chk_file eq 0 then message, 'Error opening file'

    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will open the dataset 'test_hdf5_table'
    
    dset_name = 'test_hdf5_table'

    ; get some information about the table size
    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, nrecords

    ; get the field names and record data structure
    wmb_h5tb_get_field_info, loc_id, dset_name, record_definition

    ; read the table data
    wmb_h5tb_read_table, loc_id, dset_name, databuffer

    ; print some information about the data we read
    
    print, 'Number of records read: ', nrecords
    print, 'Number of fields in the table: ', nfields
    print, 'Field names: ', tag_names(record_definition)
    print, 'First record: ', databuffer[0]
    
    ; close the file
    
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Read some records from a HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_records_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    

    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    
    start_read_position = 0
    num_read_records = 5
    
    ; overwrite the fields in the table
    
    wmb_h5tb_read_records, loc_id, $
                           dset_name, $
                           start_read_position, $
                           num_read_records, $
                           databuffer

    ; print the output data
    
    print, databuffer

    ; close the file
             
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Read a range of records from a HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_records_range_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    

    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    
    start_read_position = 0
    end_read_position = 10
    stride = 2
    
    ; overwrite the fields in the table
    
    wmb_h5tb_read_records_range, loc_id, $
                                 dset_name, $
                                 start_read_position, $
                                 end_read_position, $
                                 stride, $
                                 databuffer

    ; print the output data
    
    print, databuffer

    ; close the file
             
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Read records from an HDF5 table according to an index array
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_records_index_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    

    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    
    index_array = [5,4,6]
    
    ; overwrite the fields in the table
    
    wmb_h5tb_read_records_index, loc_id, $
                                 dset_name, $
                                 index_array, $
                                 databuffer

    print, size(databuffer)

    ; print the output data
    
    print, databuffer

    ; close the file
             
    h5f_close, fid
    
    bb = databuffer[0]
    
    print, size(bb)
    
end




;-----------------------------------------------------------------------------
;
;   Read a subset of the fields (by name) from a HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_fields_by_name_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    ; choose the fields to read
    
    read_field_names = ['testfield_int', $
                        'testfield_string']
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    
    start_read_position = 0
    num_read_records = 10
    
    ; overwrite the fields in the table
    
    wmb_h5tb_read_fields_name, loc_id, $
                               dset_name, $
                               read_field_names, $
                               start_read_position, $
                               num_read_records, $
                               databuffer

    ; print the output data
    
    print, databuffer

    ; close the file
             
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Read a subset of the fields (by index) from a HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_read_fields_by_index_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    ; choose the fields to read
    
    read_field_index = [0,2]
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    
    start_read_position = 0
    num_read_records = 10
    
    ; overwrite the fields in the table
    
    wmb_h5tb_read_fields_index, loc_id, $
                                dset_name, $
                                read_field_index, $
                                start_read_position, $
                                num_read_records, $
                                databuffer

    ; print the output data
    
    print, databuffer

    ; close the file
             
    h5f_close, fid
    
end




;-----------------------------------------------------------------------------
;
;   Append new records onto an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_append_records_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    
    ; create some new records to append
    
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
    
    n_new_records = 30
    
    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_int[*] = 77
    newdata.testfield_float[*] = 77.0
    newdata.testfield_string[*] = 'AppendData'
    newdata.testfield_long[*] = 77L
    newdata.testfield_double[*] = 77.7D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will add the new records to the dataset created previously
    
    dset_name = 'test_hdf5_table'

    ; append the record

    wmb_h5tb_append_records, loc_id, $
                             dset_name, $
                             n_new_records, $
                             newdata
    
    ; close the file          
     
    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Insert new records into an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_insert_records_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    
    ; create some new records to insert
    
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
    
    n_new_records = 20
    
    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_int[*] = 3
    newdata.testfield_float[*] = 5.0
    newdata.testfield_string[*] = 'NewData'
    newdata.testfield_long[*] = 7L
    newdata.testfield_double[*] = 9.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will add the new records to the dataset created previously
    
    dset_name = 'test_hdf5_table'

    ; insert the new records at the top of the table

    start = 0

    wmb_h5tb_insert_records, loc_id, $
                             dset_name, $
                             start, $
                             n_new_records, $
                             newdata
    
    ; close the file          
     
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Overwrite some records in an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_overwrite_records_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    
    ; create some new records to write
    
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
    
    n_new_records = 5
    
    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_int[*] = 1
    newdata.testfield_float[*] = 1.0
    newdata.testfield_string[*] = 'Overwrite'
    newdata.testfield_long[*] = 1L
    newdata.testfield_double[*] = 1.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will add the new records to the dataset created previously
    
    dset_name = 'test_hdf5_table'

    ; insert the new records at the top of the table

    start = 20

    wmb_h5tb_write_records, loc_id, $
                            dset_name, $
                            start, $
                            n_new_records, $
                            newdata
    
    ; close the file          
     
    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Overwrite some records in an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_overwrite_records_range_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    
    ; create some new records to write
    
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
    
    n_new_records = 5
    
    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_int[0] = indgen(n_new_records)
    newdata.testfield_float[*] = 1.0
    newdata.testfield_string[*] = 'Overwrite'
    newdata.testfield_long[*] = 1L
    newdata.testfield_double[*] = 1.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will add the new records to the dataset created previously
    
    dset_name = 'test_hdf5_table'



    firstrecord = 13
    lastrecord = 5
    stride = -2

    wmb_h5tb_write_records_range, loc_id, $
                                  dset_name, $
                                  firstrecord, $
                                  lastrecord, $
                                  stride, $
                                  newdata
    
    ; close the file          
     
    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Overwrite some records in an existing HDF5 table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_overwrite_records_index_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name
    
    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    
    ; create some new records to write
    
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
    
    n_new_records = 5
    
    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_int[0] = indgen(n_new_records)
    newdata.testfield_float[*] = 1.0
    newdata.testfield_string[*] = 'Overwrite'
    newdata.testfield_long[*] = 1L
    newdata.testfield_double[*] = 1.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    
    ; we will add the new records to the dataset created previously
    
    dset_name = 'test_hdf5_table'



    index = [15,17,26,25,22]

    wmb_h5tb_write_records_index, loc_id, $
                                  dset_name, $
                                  index, $
                                  newdata
    
    ; close the file          
     
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Overwrite some fields (by name) of an existing HDF5 table
;   
;   Note that due to limitations in the current IDL implementation of HDF5,
;   (as of IDL 8.2) this procedure requires the entire record to be loaded 
;   into memory before being modified and re-written to the file.
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_overwrite_fields_by_name_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    ; create some new data

    n_new_records = 10
    
    overwrite_field_names = ['testfield_long', $
                             'testfield_double']
    
    field4_value = 0L
    field5_value = 0.0D
    
    record_definition = create_struct(overwrite_field_names, $
                                      field4_value, $
                                      field5_value)

    newdata = replicate(record_definition, n_new_records)
    
    newdata.testfield_long[*] = 99L
    newdata.testfield_double[*] = 99.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    start_overwrite_position = 0
    
    ; overwrite the fields in the table
    
    wmb_h5tb_write_fields_name, loc_id, $
                                dset_name, $
                                overwrite_field_names, $
                                start_overwrite_position, $
                                n_new_records, $
                                newdata

    ; close the file
             
    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Overwrite some fields (by index) of an existing HDF5 table
;   
;   Note that due to limitations in the current IDL implementation of HDF5,
;   (as of IDL 8.2) this procedure requires the entire record to be loaded 
;   into memory before being modified and re-written to the file.
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_overwrite_fields_by_index_example

    compile_opt idl2, strictarrsubs

    ; choose an input file name

    fn = dialog_pickfile(FILTER='*.h5', /READ)
    chk_file = file_test(fn, /READ, /WRITE)
    if chk_file eq 0 then message, 'Error opening file'
    
    ; create some new data

    n_new_records = 5
    
    overwrite_field_index = [3,4]

    record_definition = {a:0L, b:0.0D}

    newdata = replicate(record_definition, n_new_records)
    
    newdata.a[*] = 99L
    newdata.b[*] = 99.0D
    
    ; open the hdf5 file
    
    fid = h5f_open(fn, /WRITE)
    loc_id = fid
    dset_name = 'test_hdf5_table'
    start_overwrite_position = 20
    
    ; overwrite the fields in the table
    
    wmb_h5tb_write_fields_index, loc_id, $
                                 dset_name, $
                                 overwrite_field_index, $
                                 start_overwrite_position, $
                                 n_new_records, $
                                 newdata

    ; close the file
             
    h5f_close, fid
    
end



;-----------------------------------------------------------------------------
;
;   Add records from one table to another table
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_add_records_from_example

    compile_opt idl2, strictarrsubs

    ; choose an output file name

    fn = dialog_pickfile(FILTER='*.h5', /WRITE)
    
    fn_info = file_info(fn)
    
    if fn_info.exists then begin
        chk_file = file_test(fn, /WRITE)
        if ~chk_file then message, 'Error opening file'
    endif
    
    if fn_info.exists then begin

        ; open the hdf5 file
        fid = h5f_open(fn, /WRITE)
    
    endif else begin
    
        ; create a new hdf5 file
        fid = h5f_create(fn)
        
    endelse
    
    ;-------------------------------------------------------------------------
    ;   Create the first table
    ;-------------------------------------------------------------------------
    
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
    
    nrecords = 50
    
    databuffer = replicate(record_definition,nrecords)
    
    databuffer.testfield_int[*] = 1
    databuffer.testfield_float[*] = 1.0
    databuffer.testfield_string[*] = 'First'
    databuffer.testfield_long[*] = 1L
    databuffer.testfield_double[*] = 1.0D

    loc_id = fid
    dset_name = 'first_table'
    table_title = 'first table'
    chunk_size = nrecords / 2
    compress = 0
    
    wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
                         
    ;-------------------------------------------------------------------------
    ;   Create the second table (with the same record structure)
    ;-------------------------------------------------------------------------
    
    databuffer.testfield_int[*] = 2
    databuffer.testfield_float[*] = 2.0
    databuffer.testfield_string[*] = 'Second'
    databuffer.testfield_long[*] = 2L
    databuffer.testfield_double[*] = 2.0D
    
    dset_name = 'second_table'
    table_title = 'second table'
    
    wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
                         
    ;-------------------------------------------------------------------------
    ;   Add records from the second table to the first table
    ;-------------------------------------------------------------------------
    
    ; add 10 records from the end of the second table to the beginning 
    ; of the first table
    
    dset_read = 'second_table'
    dset_write = 'first_table'
    read_nrecords = 10
    read_start = 40
    write_start = 0
    
    wmb_h5tb_add_records_from, loc_id, $
                               dset_read, $
                               read_start, $
                               read_nrecords, $
                               dset_write, $
                               write_start
    
    ; close the file

    h5f_close, fid
    
end


;-----------------------------------------------------------------------------
;
;   Combine two tables
;
;-----------------------------------------------------------------------------

pro wmb_h5tb_combine_tables_example

    compile_opt idl2, strictarrsubs

    ; choose an output file name

    fn = dialog_pickfile(FILTER='*.h5', /WRITE)
    
    fn_info = file_info(fn)
    
    if fn_info.exists then begin
        chk_file = file_test(fn, /WRITE)
        if ~chk_file then message, 'Error opening file'
    endif
    
    if fn_info.exists then begin

        ; open the hdf5 file
        fid = h5f_open(fn, /WRITE)
    
    endif else begin
    
        ; create a new hdf5 file
        fid = h5f_create(fn)
        
    endelse
    
    ;-------------------------------------------------------------------------
    ;   Create the first table
    ;-------------------------------------------------------------------------
    
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
    
    nrecords = 50
    
    databuffer = replicate(record_definition,nrecords)
    
    databuffer.testfield_int[*] = 1
    databuffer.testfield_float[*] = 1.0
    databuffer.testfield_string[*] = 'First'
    databuffer.testfield_long[*] = 1L
    databuffer.testfield_double[*] = 1.0D

    loc_id = fid
    dset_name = 'A_table'
    table_title = 'A table'
    chunk_size = nrecords / 2
    compress = 0
    
    wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
                         
    ;-------------------------------------------------------------------------
    ;   Create the second table (with the same record structure)
    ;-------------------------------------------------------------------------
    
    databuffer.testfield_int[*] = 2
    databuffer.testfield_float[*] = 2.0
    databuffer.testfield_string[*] = 'Second'
    databuffer.testfield_long[*] = 2L
    databuffer.testfield_double[*] = 2.0D
    
    dset_name = 'B_table'
    table_title = 'B table'
    
    wmb_h5tb_make_table, table_title, $
                         loc_id, $
                         dset_name, $
                         nrecords, $
                         record_definition, $
                         chunk_size, $
                         compress, $
                         databuffer = databuffer
                         
    ;-------------------------------------------------------------------------
    ;   Combine the two tables
    ;-------------------------------------------------------------------------
    
    ; add 10 records from the end of the second table to the beginning 
    ; of the first table
    
    loc_id1 = loc_id
    loc_id2 = loc_id
    dset_name1 = 'A_table'
    dset_name2 = 'B_table'
    dset_name3 = 'AB_table'
    new_table_title = 'A+B table'
    chunk_size = 30
    compress = 0
    
    wmb_h5tb_combine_tables, loc_id1, $
                             dset_name1, $
                             loc_id2, $
                             dset_name2, $
                             dset_name3, $
                             new_table_title, $
                             chunk_size, $
                             compress
    
    ; close the file

    h5f_close, fid
    
end


pro wmb_h5tb_examples

    wmb_h5tb_make_table_example
;    
;    wmb_h5tb_read_table_example
;
;    wmb_h5tb_read_records_example
;
;    wmb_h5tb_read_records_range_example    
;
;    wmb_h5tb_read_records_index_example
;
;    wmb_h5tb_read_fields_by_name_example
;
;    wmb_h5tb_read_fields_by_index_example
;
;    wmb_h5tb_append_records_example
;
;    wmb_h5tb_insert_records_example
;
;    wmb_h5tb_overwrite_records_example
;
;    wmb_h5tb_overwrite_records_range_example
;    
;    wmb_h5tb_overwrite_records_index_example        
;
;    wmb_h5tb_overwrite_fields_by_name_example
;
;    wmb_h5tb_overwrite_fields_by_index_example
;
;    wmb_h5tb_add_records_from_example
;
;    wmb_h5tb_combine_tables_example
    
end
