

; wmb_h5tb_data_to_record_definition
; 
; Purpose: Private function to convert an array of structures to a 
;          record definition appropriate for the wmb_h5tb_* functions.
;          
; Return value: IDL structure variable corresponding to the data.
; 
; Notes: String values will be represented in the structure as a string 
;        type field, with a length equal to the longest element in the
;        table.
;        
; WARNING: This routine makes internal copies of large segments of the
;          data, which could result in poor performance.
;


function wmb_h5tb_data_to_record_definition, data
     
    compile_opt idl2, strictarrsubs
     
    tmp_dtype = size(data, /type)
    tmp_ndims = size(data, /n_dimensions)
    tmp_nelts = size(data, /n_elements)
     
    ; this function accepts only a 1D array of structures
    if tmp_dtype ne 8 or tmp_ndims ne 1 then message, 'Incompatible data type'       
                                             
    nmembers = n_tags(data)
    tmp_tags = tag_names(data)
    
    recdef = data[0]
    
    for i = 0, nmembers-1 do begin
    
        member_type = size(data.(i), /type)
        
        case member_type of 
        
            7: begin
        
                ; grab the entire string array from the data
                
                string_data = data.(i)
                tmpa = max(strlen(string_data),tmpind) > 1
                tmpstr = string(replicate(95b,tmpa)) 
                
                recdef.(i) = tmpstr
        
            end
        
            8: begin
            
                ; this is a structure
                
                ; what are the dimensions of the structure?
                
                tmpstruct = recdef.(i)
                
                tmpstruct_nelts = n_elements(tmpstruct)
                tmpstruct_ndims = size(tmpstruct, /n_dimensions)
                
                ; is this a simple structure, or a 1D array of structures?
                
                if tmpstruct_ndims eq 1 then begin
                
                    if tmpstruct_nelts eq 1 then begin
                
                        ; this is a simple structure
                
                        struct_data = data.(i)
                        
                        tmpstruct = $
                            wmb_h5tb_data_to_record_definition(struct_data)
                            
                        recdef.(i) = tmpstruct
                        
                    endif else begin
                    
                        ; we will not allow an array of structures to be 
                        ; nested within an individual record of the table
                        
                        message, 'Incompatible data type'
                    
                    endelse
                    
                endif else begin
                
                    message, 'Incompatible data type'
                
                endelse
        
            end
    
            else: 
            
        endcase
    
    endfor            
   
    return, recdef                              
                                             
end