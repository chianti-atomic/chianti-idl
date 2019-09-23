
; wmb_h5lt_find_attribute
; 
; Purpose: Inquires if an attribute named attr_name exists attached to 
;          the object loc_id.
; 
; Return:  Returns one if the attribute was found, and zero otherwise.
; 

function wmb_h5lt_find_attribute, loc_id, attr_name, match_index = match_index

    compile_opt idl2, strictarrsubs

    num_attr = h5a_get_num_attrs(loc_id)
    
    matched = 0
    match_index = 0
    
    for i = 0, num_attr-1 do begin
    
        aid = h5a_open_idx(loc_id, i)
        
        aname = h5a_get_name(aid)
        
        if aname eq attr_name then begin
        
            matched = 1
            match_index = i
            break
            
        endif
    
    endfor
    
    return, matched
    
end