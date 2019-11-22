

; wmb_h5tb_get_field_info
; 
; Purpose: Get information about fields.  This function retrieves the 
;          record_definition structure, which defines the data type of
;          each record in the table.
;


pro wmb_h5tb_get_field_info, loc_id, $
                             dset_name, $
                             record_definition

    compile_opt idl2, strictarrsubs

    ; open the dataset
                      
    did = h5d_open(loc_id, dset_name)
    
    ; get the datatype
    
    tid = h5d_get_type(did)
    
    ; get the structure corresponding to the record data type
    
    record_definition = wmb_h5tb_datatype_to_idl_structure(tid)
    
    ; close
    
    h5t_close, tid
    h5d_close, did
    
end

