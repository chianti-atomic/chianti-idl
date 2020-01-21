
;
; wmb_h5tba_get_title
; 
; Purpose: Read the table title.
;


pro wmb_h5tba_get_title, loc_id, table_title

    compile_opt idl2, strictarrsubs

    ; get the "TITLE" attribute
                      
    wmb_h5lt_get_attribute_disk, loc_id, 'TITLE', table_title
    
end

