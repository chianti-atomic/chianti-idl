
; wmb_h5o_close
; 
; Purpose: Closes the group, dataset, or named datatype specified by object_id.
;          This function is the companion to wmb_h5o_open, and has the same 
;          effect as calling H5Gclose, H5Dclose, or H5Tclose.  
; 

pro wmb_h5o_close, obj_id

    compile_opt idl2, strictarrsubs

    ; What type of object is specified?
    
    obj_type = h5i_get_type(obj_id)
     
    case obj_type of
    
        'GROUP': h5g_close, obj_id
        
        'DATASET': h5d_close, obj_id
        
        'DATATYPE': h5t_close, obj_id
        
        else: 
        
    endcase

end