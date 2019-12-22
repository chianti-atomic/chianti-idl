

; wmb_h5tb_create_type
; 
; Purpose: Private function that creates a memory type ID.
;


function wmb_h5tb_create_type, loc_id, $
                               dset_name
    
    compile_opt idl2, strictarrsubs
                              
    wmb_h5tb_get_field_info, loc_id, dset_name, record_definition
                             
    mem_type_id = h5t_idl_create(record_definition)
    
    return, mem_type_id
    
end
