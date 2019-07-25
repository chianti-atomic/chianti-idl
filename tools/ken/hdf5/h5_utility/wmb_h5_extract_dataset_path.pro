

; wmb_h5_extract_dataset_path
;
; Purpose: Extract the dataset name and full group name from an HDF5 path.
;
; The full group name and dataset name are output to the variables 
; full_group_name and dataset_name.


pro wmb_h5_extract_dataset_path, full_path, full_group_name, dataset_name

    compile_opt idl2, strictarrsubs


    ; strip off the leading '/' in full_path if it is present

    if strmid(full_path,0,1) eq '/' then begin

        full_path = strmid(full_path,1)

    endif


    slash_pos = strpos(full_path, '/', /reverse_search)

    if slash_pos eq -1 then full_group_name = '/' $
                 else full_group_name = strmid(full_path, 0, slash_pos)

    if slash_pos eq -1 then dataset_name = full_path $
                       else dataset_name = strmid(full_path, slash_pos+1)


    full_dataset_name = '/' + full_group_name + '/' + dataset_name

    
end