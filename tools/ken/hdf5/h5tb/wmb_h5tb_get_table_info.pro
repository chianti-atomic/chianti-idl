

; wmb_h5tb_get_table_info
; 
; Purpose: Gets the number of records and fields of a table.  Optionally, 
;          this procedure will also return the values of the TITLE and CLASS
;          attributes of the dataset.
;


pro wmb_h5tb_get_table_info, loc_id, $
                             dset_name, $
                             nfields, $
                             nrecords, $
                             title=title, $
                             class=class

    compile_opt idl2, strictarrsubs

    ; open the dataset
                      
    did = h5d_open(loc_id, dset_name)
    
    ; get the datatype
    
    tid = h5d_get_type(did)
    
    ; get the number of members
    
    num_members = h5t_get_nmembers(tid)
    
    ; get the number of fields
    
    nfields = ulong64(num_members)
    
    ; get the dataspace handle
    
    sid = h5d_get_space(did)
    
    ; get dimensions
    
    dims = h5s_get_simple_extent_dims(sid)
    
    ; terminate access to the dataspace
    
    h5s_close, sid
    
    ; get the number of records
    
    nrecords = ulong64(dims[0])
    
    ; get the table title attribute
    
    wmb_h5tba_get_title, did, table_title
    
    ; get the table class attribute
    
    wmb_h5lt_get_attribute_disk, did, 'CLASS', did_class
    
    title = table_title
    class = did_class
    
    ; close
    
    h5t_close, tid
    h5d_close, did
    
end
