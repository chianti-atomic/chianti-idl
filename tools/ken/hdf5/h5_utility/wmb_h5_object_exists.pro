;+
;   Determines if an object with a given name exists at a given location.
;   
;   :Returns:
;       1 if exists, 0 if it doesn't
;
;   :Params:
;       loc_id : in, required, type=long
;           file or group identifier
;       name : in, required, type=string
;           name of object to check

function wmb_h5_object_exists, loc_id, name

    compile_opt idl2, strictarrsubs

    nobjs = h5g_get_num_objs(loc_id)

    for i = 0, nobjs-1 do begin

        if h5g_get_obj_name_by_idx(loc_id, i) eq name then return, 1

    endfor

    return, 0

end