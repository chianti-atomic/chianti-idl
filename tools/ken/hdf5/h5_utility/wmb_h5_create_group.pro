

; wmb_h5_create_group
;
; Purpose: Create a new group within an HDF5 file.
; 
; Note that groupname is case sensitive.
;
; Returns 1 if successful.


function wmb_h5_create_group, filename, groupname

    compile_opt idl2, strictarrsubs

    if filename eq '' then message, 'Error: Invalid file name'
    if groupname eq '' then message, 'Error: Invalid dataset name'


    ; check if filename exists, is writable, and is a valid hdf5 file

    if ~wmb_h5_file_test(filename, /WRITE) then begin
        
        message, 'Error: invalid HDF5 file'
        return, 0
        
    endif
    
    
    ; strip off the leading '/' in groupname if it is present
    
    if strmid(groupname,0,1) eq '/' then groupname = strmid(groupname,1)
    
    
    ; test if the group exists already
    
    if wmb_h5_group_exists(filename, groupname) then return, 1
    
    
    ; the group does not exist - create new groups
    
    tokens = strsplit(groupname, '/', /extract, /preserve_null, count=ntokens)

    n_groupnames = ntokens

    ; open the file

    file_id = h5f_open(filename, /WRITE)

    loc_id = file_id

    if n_groupnames gt 0 then begin

        group_names = tokens

        group_ids = lonarr(n_groupnames)

        for i = 0, n_groupnames-1 do begin

            if (wmb_h5_object_exists(loc_id, group_names[i])) then begin

                ; open an existing group
                loc_id = h5g_open(loc_id, group_names[i])

            endif else begin

                ; create a new group
                loc_id = h5g_create(loc_id, group_names[i])

            endelse

            group_ids[i] = loc_id

        endfor

    endif


    ; close the groups

    if n_groupnames gt 0 then begin

        for i = n_groupnames-1, 0, -1 do h5g_close, group_ids[i]

    endif


    ; close the file
    
    h5f_close, file_id
    
    return, wmb_h5_group_exists(filename, groupname)
    
end