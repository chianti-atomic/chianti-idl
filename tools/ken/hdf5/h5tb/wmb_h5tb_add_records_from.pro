
;
; wmb_h5tb_add_records_from
; 
; Purpose: Add records from first table to second table. 
;
; Description: wmb_h5tb_add_records_from adds records from a dataset named 
;              dset_name1 to a dataset named dset_name2. Both tables are 
;              attached to the object specified by the identifier loc_id.  
;              
; Parameters:
; 
; loc_id
;     IN: Identifier of the file or group in which the table is located. 
; dset_name1
;     IN: The name of the dataset to read the records. 
; start1
;     IN: The position to read the records from the first table 
; read_nrecords
;     IN: The number of records to read from the first table. 
; dset_name2
;     IN: The name of the dataset to write the records. 
; start2
;     IN: The position to write the records on the second table 
;

pro wmb_h5tb_add_records_from, loc_id, $
                               dset_name1, $
                               start1, $
                               read_nrecords, $
                               dset_name2, $
                               start2

    compile_opt idl2, strictarrsubs

    dims = ulon64arr(1)
    mem_size = ulon64arr(1)
    count = ulon64arr(1)
    offset = ulon64arr(1)


    ;-------------------------------------------------------------------------
    ; get information about the two tables
    ;-------------------------------------------------------------------------
    
    wmb_h5tb_get_table_info, loc_id, dset_name1, nfields1, nrecords1
    wmb_h5tb_get_table_info, loc_id, dset_name2, nfields2, nrecords2
    
    wmb_h5tb_get_field_info, loc_id, dset_name1, recdef1
    wmb_h5tb_get_field_info, loc_id, dset_name2, recdef2
    
    
    ;-------------------------------------------------------------------------
    ; compare the record definitions to ensure that they match
    ;-------------------------------------------------------------------------

    recdef_matched = wmb_compare_struct(recdef1, $
                                        recdef2, $
                                        /compare_tag_names)
                                             
    if recdef_matched eq 0 then message, 'Unmatched record definitions'

    
    ;-------------------------------------------------------------------------
    ; open the first dataset
    ;-------------------------------------------------------------------------

    ; open the dataset
    
    did_1 = h5d_open(loc_id, dset_name1)
    
    
    ; get the file data space handle
    
    sid_1 = h5d_get_space(did_1)
    
    

    ;-------------------------------------------------------------------------
    ; read data from 1st table
    ;-------------------------------------------------------------------------
    
    ; define a hyperslab in the dataset of the size of the records
    
    offset[0] = start1
    count[0]  = read_nrecords
    
    h5s_select_hyperslab, sid_1, offset, count, /RESET
    
    ; create a memory dataspace handle
    
    mem_size[0] = read_nrecords
    m_sid = h5s_create_simple(mem_size)
    
    tmp_buf = h5d_read(did_1, MEMORY_SPACE = m_sid, FILE_SPACE = sid_1)
    
    
    ;-------------------------------------------------------------------------
    ; save data from 1st table into the 2nd table
    ;-------------------------------------------------------------------------

    wmb_h5tb_insert_records, loc_id, dset_name2, start2, read_nrecords, tmp_buf
    
    
    ;-------------------------------------------------------------------------
    ; release resources from 1st table
    ;-------------------------------------------------------------------------
    
    h5s_close, m_sid
    h5s_close, sid_1
    h5d_close, did_1
    
end



