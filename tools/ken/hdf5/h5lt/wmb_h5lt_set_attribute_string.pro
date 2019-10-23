

; wmb_h5lt_set_attribute_string
; 
; Purpose: Creates and writes a string attribute named attr_name and attaches
;          it to the object specified by the name obj_name.
; 
; Comments: If the attribute already exists, it is overwritten


pro wmb_h5lt_set_attribute_string, loc_id, obj_name, attr_name, attr_data

    compile_opt idl2, strictarrsubs

    ; Open the object

    obj_id = wmb_h5o_open(loc_id, obj_name)
    
    if obj_id lt 0 then message, 'Could not open the specified object'

    ; Create the attribute

    attr_type = h5t_idl_create(attr_data)
    
    attr_space_id = h5s_create_scalar()
    
    ; Check if the attribute already exists
    
    has_attr = wmb_h5lt_find_attribute(obj_id, attr_name)
    
    ; Delete the attribute if it already exists
    
    if has_attr then h5a_delete, obj_id, attr_name
    
    ; Create and write the attribute
    
    attr_id = h5a_create(obj_id, attr_name, attr_type, attr_space_id)
    
    h5a_write, attr_id, attr_data
    
    ; Close
    
    h5a_close, attr_id
    h5s_close, attr_space_id
    h5t_close, attr_type
    
    ; Close the object
    
    wmb_h5o_close, obj_id
    
end
