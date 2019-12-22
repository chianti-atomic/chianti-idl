
; wmb_h5o_open
; 
; Purpose: Opens a group, dataset, or named datatype specified by a location, 
;          loc_id, and a path name, name, in an HDF5 file. 
; 
; Return:  Returns an object identifier for the opened object if successful; 
;          otherwise returns a negative value. 
; 

function wmb_h5o_open, loc_id, name

    compile_opt idl2, strictarrsubs

    ; What type of object is specified?
    
    tmp_info = h5g_get_objinfo(loc_id, name)

    obj_type = tmp_info.TYPE
     
    case obj_type of
    
        'GROUP': begin
        
            gid = h5g_open(loc_id, name)
            return, gid
        
        end
        
        'DATASET': begin
        
            did = h5d_open(loc_id, name)
            return, did
        
        end
        
        'TYPE': begin
        
            datatype_id = h5t_open(loc_id, name)
            return, datatype_id
        
        end
        
        else: 
        
    endcase

    ; if the function has not returned yet, then either the structure tag
    ; was not found, or the object type is not a group, dataset, or named
    ; datatype
    
    return, -1
    
end