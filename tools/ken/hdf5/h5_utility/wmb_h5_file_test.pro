

; wmb_h5_file_test
;
; Purpose: Test that filename is a valid HDF5 file.
;
; Returns 1 if the file is valid.


function wmb_h5_file_test, filename, $
                           write=chk_write

    compile_opt idl2, strictarrsubs

    if N_elements(chk_write) eq 0 then chk_write = 0

    ; check if filename is a valid hdf5 file

    f_info = file_info(filename)
    fn_exists = f_info.exists
    fn_writeable = f_info.write

    if fn_exists then fn_is_hdf5 = h5f_is_hdf5(filename) $
                 else fn_is_hdf5 = 0

    if chk_write eq 0 then begin
        
        return, fn_is_hdf5

    endif else begin
        
        return, (fn_is_hdf5 && fn_writeable)

    endelse

end
