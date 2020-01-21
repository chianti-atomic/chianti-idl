

; wmb_h5_form_dataset_path
;
; Purpose: Generate the full dataset name including the group name.
;
; Returns the full dataset path including the group name.


function wmb_h5_form_dataset_path, full_group_name, dataset_name

    compile_opt idl2, strictarrsubs


    ; strip off the leading '/' in dataset_name if it is present

    if strmid(dataset_name,0,1) eq '/' then begin

        dataset_name = strmid(dataset_name,1)

    endif


    ; strip off the leading '/' in full_group_name if it is present

    if strmid(full_group_name,0,1) eq '/' then begin
        
        full_group_name = strmid(full_group_name,1)

    endif
    
    
    ; strip off the trailing '/' in full_group_name if it is present

    if strmid(full_group_name,0,1,/REVERSE_OFFSET) eq '/' then begin
        
        full_group_name = strmid(full_group_name,0,strlen(full_group_name)-1)

    endif


    if full_group_name eq '' then begin
        
        full_dataset_path = '/' + dataset_name
        
    endif else begin
        
        full_dataset_path = '/' + full_group_name + '/' + dataset_name

    endelse
    

    return, full_dataset_path
    
end