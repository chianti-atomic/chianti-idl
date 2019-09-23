

; wmb_h5tb_datatype_to_idl_structure
; 
; Purpose: Private function to convert the record datatype from an HDF5 
;          table (an HDF5 compound datatype) to an IDL structure variable.
;          
; Return value: IDL structure variable corresponding to the datatype.
; 
; Notes: String values will be represented in the structure as a string 
;        type field, filled with underscore "_" characters to represent the
;        size in bytes of the field.
;


function wmb_h5tb_datatype_to_idl_structure, datatype_id
     
    compile_opt idl2, strictarrsubs
     
    ; get the structure corresponding to the record data type
                                       
    n_tid = h5t_idltype(datatype_id, structure=recdef)
    
    if n_tid ne 8 then message, 'Incompatible data type'       
                                             
    nmembers = h5t_get_nmembers(datatype_id)
    
    for i = 0, nmembers-1 do begin
    
        m_class = h5t_get_member_class(datatype_id, i)
        
        if m_class eq 'H5T_STRING' then begin
        
            fieldtp = h5t_get_member_type(datatype_id, i)
            fieldsz = h5t_get_size(fieldtp)
            tmpstr = string(replicate(95b,fieldsz-1)) 
            recdef.(i) = tmpstr
            
            h5t_close, fieldtp
        
        endif
    
        if m_class eq 'H5T_COMPOUND' then begin
        
            fieldtp = h5t_get_member_type(datatype_id, i)
            tmpstruct = wmb_h5tb_datatype_to_idl_structure(fieldtp)
            recdef.(i) = tmpstruct
        
            h5t_close, fieldtp
        
        endif
    
    endfor            
   
    return, recdef                              
                                             
end