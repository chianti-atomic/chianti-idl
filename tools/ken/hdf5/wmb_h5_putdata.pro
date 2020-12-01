; docformat = 'rst'

;+
;   Write a variable or attribute to an HDF 5 file with a simple notation.
;   
;   This code is based on Michael Galloy's mg_h5_putdata function, which 
;   is part of mg_lib hosted at http://github.com/mgalloy.
;
;   :Categories:
;       file i/o, hdf5, sdf
;
;   :Author:
;       Mark Bates
;
;-


;+
;   Create a reference given a string to the reference object.
;
;   :Private:
;
;   :Returns:
;       reference identifier
;
;   :Params:
;       file_id : in, required, type=ulong64
;           file identifier
;       obj_to_reference : in, required, type=string
;           string specifying full path of object to reference
;       refgroup : out, required, type=ulong64
;           identifier for group containing reference object, must be 
;           closed by caller
;-

function wmb_h5_putdata_getreference, file_id, $
                                      obj_to_reference, $
                                      refgroup

    compile_opt idl2, strictarrsubs

    slash_pos = strpos(obj_to_reference, '/', /reverse_search)

    if slash_pos eq -1 then groupname = '/' $
                       else groupname = strmid(obj_to_reference, 0, slash_pos)

    refgroup = h5g_open(file_id, groupname)

    if slash_pos eq -1 then objname = obj_to_reference $
                       else objname = strmid(obj_to_reference, slash_pos+1)
    
    ref_id = h5r_create(refgroup, objname)

    return, ref_id

end


;+
;   Determines if an object with a given name exists at a given location.
;
;   :Private:
;
;   :Returns:
;       1 if exists, 0 if it doesn't
;
;   :Params:
;       loc_id : in, required, type=long
;           file or group identifier
;       name : in, required, type=string
;           name of object to check
;-

function wmb_h5_putdata_varexists, loc_id, name

    compile_opt idl2, strictarrsubs

    nobjs = h5g_get_num_objs(loc_id)
    
    for i = 0, nobjs-1 do begin
    
        if h5g_get_obj_name_by_idx(loc_id, i) eq name then return, 1
        
    endfor

    return, 0
    
end


; +
;   Write a variable to a file.
;
;   :Private:
;
;   :Params:
;       fildename : in, required, type=long
;           HDF5 filename to write variable into
;       name : in, required, type=string
;           name of variable in HDF5 file
;       data : in, optional, type=any
;           IDL variable to write
;
;   :Keywords:
;       reference : in, optional, type=boolean
;           set to indicate that `data` is a reference to an attribute/variable 
;           in the file instead of actual data
;-

pro wmb_h5_putdata_putvariable, file_id, name, data, reference=reference

    compile_opt idl2, strictarrsubs

    ; determine into which group we are writing the data
    
    tokens = strsplit(name, '/', /extract, /preserve_null, count=ntokens)

    n_groupnames = ntokens-1
 
    ; tokens contains an array of strings - the last string in the array
    ; is the dataset name

    slash_pos = strpos(name, '/', /reverse_search)

    if slash_pos eq -1 then fullgroupname = '/' $
                       else fullgroupname = strmid(name, 0, slash_pos)

    if slash_pos eq -1 then objname = name $
                       else objname = strmid(name, slash_pos+1)

    if objname eq '' then message, 'Invalid object name'

    ; determine the group ids and the location id that we are writing to

    loc_id = file_id

    if n_groupnames gt 0 then begin

        group_names = tokens[0:n_groupnames-1]
        
        group_ids = lonarr(n_groupnames)   
    
        for i = 0, n_groupnames-1 do begin
        
            if (wmb_h5_putdata_varexists(loc_id, group_names[i])) then begin
            
                ; open an existing group
                loc_id = h5g_open(loc_id, group_names[i])

            endif else begin

                ; create a new group
                loc_id = h5g_create(loc_id, group_names[i])

            endelse
            
            group_ids[i] = loc_id
            
        endfor

    endif

    ; loc_id now points to where the data will be written

    
    ; check the IDL data type and dimension
    
    tmp_dtype = size(data, /type)
    tmp_ndims = size(data, /n_dimensions)

    if tmp_dtype eq 8 and tmp_ndims eq 1 then begin
    
        ; we are writing a table
    
        table_title = objname
        dset_name = objname
        nrecords = size(data, /n_elements)
        record_definition = wmb_h5tb_data_to_record_definition(data)
        chunk_size = round(nrecords/100) > 1
        compress = 0
    
        wmb_h5tb_make_table, table_title, $
                             loc_id, $
                             dset_name, $
                             nrecords, $
                             record_definition, $
                             chunk_size, $
                             compress, $
                             databuffer = data

    endif else begin

        ; check if we are writing a reference - in this case the data variable
        ; is a string containing the full path to the object to reference
        
        if keyword_set(reference) then begin
    
            ref_id = wmb_h5_putdata_getreference(file_id, data, refgroup)
    
            ; When we write the data we will be writing the reference id - 
            ; note that a reference id does not need to be closed.  We will 
            ; need to close the refgroup id however.

            datatypeId = h5t_reference_create()
        
        endif else begin
    
            ; test for the case of a string array, and use the longest string 
            ; element for datatype creation
      
            if tmp_dtype eq 7 and tmp_ndims gt 0 then begin
      
                ; data is a string array - find the longest string in the array
                tmpa = max(strlen(data),tmpind)
                tmp_longeststring =  data[tmpind]
                datatypeId = h5t_idl_create(tmp_longeststring)
      
            endif else begin
    
                datatypeId = h5t_idl_create(data)
    
            endelse
    
        endelse

        ; scalars and arrays are created differently
        
        if (size(data, /n_dimensions) eq 0L) then begin
    
            dataspaceId = h5s_create_scalar()
    
        endif else begin
    
            dataspaceId = h5s_create_simple(size(data, /dimensions))
    
        endelse

        if wmb_h5_putdata_varexists(loc_id, objname) then begin
        
            datasetId = h5d_open(loc_id, objname)
        
        endif else begin
        
            datasetId = h5d_create(loc_id, objname, datatypeId, dataspaceId)
            
        endelse

        if keyword_set(reference) then h5d_write, datasetId, ref_id $
                                  else h5d_write, datasetId, data

        h5t_close, datatypeId
        h5s_close, dataspaceId
        h5d_close, datasetId

    endelse

    if n_groupnames gt 0 then begin
    
        for i = n_groupnames-1, 0, -1 do h5g_close, group_ids[i]

    endif

    if keyword_set(reference) then h5g_close, refgroup
    
end




;+
;   Write an attribute to a file.
;
;   :Private:
;
;   :Params:
;       filename : in, required, type=long
;           HDF5 filename to write variable into
;       objpath : in, required, type=string
;           full path of object in HDF5 file
;       attname : in, required, type=string
;           name of attribute to write
;       attvalue : in, optional, type=any
;           IDL variable to write
;
;   :Keywords:
;       reference : in, optional, type=boolean
;           set to indicate that `data` is a reference to an attribute/variable
;           in the file instead of actual data
;-

pro wmb_h5_putdata_putattribute, file_id, objpath, attname, attvalue, $
                                 reference=reference

    compile_opt idl2, strictarrsubs

    obj_id = wmb_h5o_open(file_id, objpath)


    ; create the datatype and get the reference id if necessary

    if keyword_set(reference) then begin
    
        ref_id = wmb_h5_putdata_getreference(file_id, attvalue, refgroup)

        ; when we write the attribute value we will be writing the 
        ; reference id

        datatypeId = h5t_reference_create()

    endif else begin
    
        datatypeId = h5t_idl_create(attvalue)

    endelse


    ; create the dataspace - scalars and arrays are created differently
    
    if size(attvalue, /n_dimensions) eq 0 then begin
    
        dataspaceId = h5s_create_scalar()
        
    endif else begin

        dataspaceId = h5s_create_simple(size(attvalue, /dimensions))
        
    endelse

    ; if the attribute already exists then delete it
    
    chk_exists = wmb_h5lt_find_attribute(obj_id, attname)
    if chk_exists then h5a_delete, obj_id, attname

    ; create the attribute

    attributeId = h5a_create(obj_id, attname, datatypeId, dataspaceId)

    ; write the attribute

    if keyword_set(reference) then h5a_write, attributeId, ref_id $
                              else h5a_write, attributeId, attvalue
    
    ; close 
    
    h5t_close, datatypeId
    h5s_close, dataspaceId
    h5a_close, attributeId
    if keyword_set(reference) then h5g_close, refgroup                 
    wmb_h5o_close, obj_id

end


;+
;   Write data to a file.
;
;   :Params:
;       file_id : in, required, type=long
;           File id of an HDF5 file to write the data into.
;           Obtained from h5f_open or h5f_create.
;       name : in, required, type=string
;           name of variable in HDF5 file
;       data : in, optional, type=any
;           IDL variable to write
;
;   :Keywords:
;       reference : in, optional, type=boolean
;           set to indicate that `data` is a reference to an 
;           attribute/variable in the file instead of actual data
;-

pro wmb_h5_putdata, file_id, name, data, reference=reference

    compile_opt idl2, strictarrsubs

    ; check the file id
    
    idtype = h5i_get_type(file_id)
    if idtype ne 'FILE' then message, 'Invalid HDF5 file ID'

    ; strip off the leading '/' in name if it is present
    
    if strmid(name,0,1) eq '/' then name = strmid(name,1)
    
    ; check if we are writing an attribute
    
    dotPos = strpos(name, '.', /reverse_search)
  
    if (dotPos eq -1) then begin
      
        ; write a variable
        
        wmb_h5_putdata_putvariable, file_id, name, data, reference=reference
      
    endif else begin
  
        ; write an attribute
        
        objpath = strmid(name, 0, dotPos)
        attname = strmid(name, dotPos + 1L)
        
        wmb_h5_putdata_putattribute, file_id, objpath, attname, data, $
                                     reference=reference
                                    
    endelse
    
end


;+
;   Simple example of wmb_h5_putdata usage.
;-

pro wmb_h5_putdata_example

    compile_opt idl2, strictarrsubs

    fn = dialog_pickfile(FILTER='*.h5', /WRITE)
    
    fn_info = file_info(fn)
    
    if fn_info.exists then begin
        chk_file = file_test(fn, /WRITE)
        if ~chk_file then message, 'Error opening file'
    endif
    
    if fn_info.exists then begin
    
        ; open the hdf5 file
        fid = h5f_open(fn, /WRITE)
    
    endif else begin
    
        ; create a new hdf5 file
        fid = h5f_create(fn)
        
    endelse
    
    wmb_h5_putdata, fid, 'scalar', 4.0
    wmb_h5_putdata, fid, 'array', findgen(10)
    wmb_h5_putdata, fid, 'reference', 'scalar', /reference
    wmb_h5_putdata, fid, 'group/another_scalar', 1.0
    wmb_h5_putdata, fid, 'group/another_array', findgen(10)
    wmb_h5_putdata, fid, 'array.attribute', 'Attribute of an array'
    wmb_h5_putdata, fid, 'ref2', 'group/another_array', /reference
    
    ; get the reference and check the value of the referenced object
    did = h5d_open(fid, 'ref2')
    buf = h5d_read(did)
    
    ; buf should now be the reference id
    
    obj_did = h5r_dereference(fid, buf)
    newbuf = h5d_read(obj_did)
    
    print, newbuf
    
    h5d_close, did
    h5d_close, obj_did
    h5f_close, fid

end
