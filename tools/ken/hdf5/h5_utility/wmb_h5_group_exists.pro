

; wmb_h5_group_exists
;
; Purpose: Verify that a given group exists within an HDF5 file.
; 
; Note that groupname is case sensitive.
;
; Returns 1 if the group exists.


function wmb_h5_group_exists, filename, groupname

    compile_opt idl2, strictarrsubs

    if filename eq '' then message, 'Error: Invalid file name'
    if groupname eq '' then message, 'Error: Invalid dataset name'


    ; check if filename is a valid hdf5 file

    if ~wmb_h5_file_test(filename) then begin
        
        message, 'Error: invalid HDF5 file'
        return, 0
        
    endif


    ; strip off the leading '/' in groupname if it is present

    if strmid(groupname,0,1) eq '/' then groupname = strmid(groupname,1)


    ; open the HDF5 file and parse the contents

    wmb_h5_list, filename, filter=groupname, output=h5listoutput
   
    output_ndims = size(h5listoutput,/n_dimensions)
    outputdims = size(h5listoutput,/dimensions)
    
    ; if nothing was found matching the filter, then a 1D array will be
    ; returned
    
    if output_ndims ne 2 then return, 0
    
    ; if matches to the filter are found, a 2D string array is returned
    
    n_found = outputdims[1] - 1
    if n_found eq 0 then return, 0
    
    groupname_found = 0
    
    for i = 0, n_found-1 do begin
        
        if h5listoutput[0,i+1] eq 'group' then begin
            
            tmpstr = h5listoutput[1,i+1]
            
            if strmid(tmpstr,0,1) eq '/' then tmpstr = strmid(tmpstr,1)
            
            if tmpstr eq groupname then groupname_found = 1
            
        endif

    endfor
    
    return, groupname_found
    
end