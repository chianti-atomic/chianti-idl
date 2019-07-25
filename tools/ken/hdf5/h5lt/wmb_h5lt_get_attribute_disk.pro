
;
; wmb_h5lt_get_attribute_disk
; 
; Purpose: Reads an attribute named attr_name with the datatype stored on disk.
; 



pro wmb_h5lt_get_attribute_disk, loc_id, attr_name, attr_out

    compile_opt idl2, strictarrsubs

    ; Open the attribute

    attr_id = h5a_open_name(loc_id, attr_name)

    ; Get the type

    attr_type = h5a_get_type(attr_id)
    
    ; Read the attribute data
    
    attr_out = h5a_read(attr_id, attr_type)
    
    ; Close
    
    h5t_close, attr_type
    h5a_close, attr_id
    
end

