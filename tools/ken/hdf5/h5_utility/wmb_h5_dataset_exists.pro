

; wmb_h5_dataset_exists
;
; Purpose: Verify that a given dataset exists within an HDF5 file.
; 
; Note that full_group_name and dataset_name are case sensitive.
;
; Returns 1 if the dataset exists.


function wmb_h5_dataset_exists, filename, full_group_name, dataset_name

    compile_opt idl2, strictarrsubs

    if filename eq '' then message, 'Error: Invalid file name'
    if dataset_name eq '' then message, 'Error: Invalid dataset name'
    if full_group_name eq '' then full_group_name = '/'


    ; check if filename is a valid hdf5 file

    if ~wmb_h5_file_test(filename) then begin
        
        message, 'Error: invalid HDF5 file'
        return, 0
        
    endif
    

    full_dataset_name = wmb_h5_form_dataset_path(full_group_name, dataset_name)


    ; open the HDF5 file and parse the contents

    wmb_h5_list, filename, filter=full_dataset_name, output=h5listoutput

    output_ndims = size(h5listoutput,/n_dimensions)
    outputdims = size(h5listoutput,/dimensions)
    
    ; if nothing was found matching the filter, then a 1D array will be
    ; returned
    
    if output_ndims ne 2 then return, 0
    
    ; if matches to the filter are found, a 2D string array is returned
    
    n_found = outputdims[1] - 1
    if n_found eq 0 then return, 0
    
    dataset_found = 0
    
    for i = 0, n_found-1 do begin
        
        tmptype = h5listoutput[0,i+1]
        tmpname = h5listoutput[1,i+1]
        
        if tmptype eq 'dataset' and tmpname eq full_dataset_name then begin
            
            dataset_found = 1
            
        endif
        
    endfor
    
    return, dataset_found
    
end