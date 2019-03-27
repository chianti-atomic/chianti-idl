; docformat = 'rst'

;+
;   Read a variable or attribute from an HDF 5 file with a simple notation.
;   
;   This code is based on Michael Galloy's mg_h5_getdata function, which 
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
;   Compute the H5D_SELECT_HYPERSLAB arguments from the bounds.
;
;   :Private:
;
;   :Params:
;       bounds : in, required, type=lonarr(3,ndims)
;           bounds
;
;   :Keywords:
;       start : out, optional, type=lonarr(ndims)
;           input for start argument to H5S_SELECT_HYPERSLAB
;       count : out, optional, type=lonarr(ndims)
;           input for count argument to H5S_SELECT_HYPERSLAB
;       stride : out, optional, type=lonarr(ndims)
;           input for stride keyword to H5S_SELECT_HYPERSLAB
;-

pro wmb_h5_getdata_computeslab, bounds, $
                                startselect, $
                                dims, $
                                stride
                               
    compile_opt idl2, strictarrsubs


    if size(bounds, /n_dimensions) eq 1 then begin
   
        ndims = 1
    
    endif else begin
    
        ndims = (size(bounds, /dimensions))[1]
        
    endelse

    ; startselect, endselect, stride, and dims are all 1D arrays

    startselect = reform(bounds[0,*])
    
    endselect = reform(bounds[1,*])
    
    stride = reform(bounds[2,*])

    dims = ceil( (endselect-startselect+1) / double(stride), /L64 ) > 1
    
end


;+
;   Reads data in a dataset.
;
;   :Private:
;
;   :Returns:
;       value of data read from dataset
;
;   :Params:
;       fileId : in, required, type=long
;           HDF 5 indentifier of the file
;       variable : in, required, type=string
;           string navigating the path to the dataset
;
;   :Keywords:
;       bounds : in, optional, type=lonarr(3,ndims)
;           gives start value, end value, and stride for each dimension 
;           of the variable
;       empty : out, optional, type=boolean
;           set to a named variable to return whether the dataset is empty
;-

function wmb_h5_getdata_getvariable, file_id, variable, bounds=bounds, $
                                     empty=empty
                                    
    compile_opt idl2, strictarrsubs

    ; open the dataset
    did = h5d_open(file_id, variable)
    
    ; get the dataspace id
    sid = h5d_get_space(did)

    ; get the data rank
    datarank = h5s_get_simple_extent_ndims(sid)

    ; get the number of elements in the dataspace
    npoints = h5s_get_simple_extent_npoints(sid)

    ; test for scalar
    chkscalar = (datarank eq 0)

    ; test if the bounds have been specified
    chkbounds = (N_elements(bounds) ne 0)

    ; read the data
    
    if (npoints eq 0) then begin
    
        empty = 1
        data = !NULL

    endif else if chkscalar or (npoints eq 1) or ~chkbounds then begin
    
        empty = 0
        data = h5d_read(did)
    
    endif else begin
    
        ; this is the case when the data is an array with greater than
        ; one element, and the bounds have been specified
        
        ; bounds is specified as a (3,rank) array containing the start element,
        ; end element, and stride for each dimension of the array

        ; from the bounds, compute the dimensions of the new array and
        ; the inputs to the h5s_compute_hyperslab function
        
        empty = 0
        
        wmb_h5_getdata_computeslab, bounds, start, dims, stride

        m_sid = h5s_create_simple(dims)                  

        h5s_select_hyperslab, sid, start, dims, STRIDE = stride, /RESET
                              
        data = h5d_read(did, file_space=sid, memory_space=m_sid)
                        
        h5s_close, m_sid
    
    endelse

    h5s_close, sid
    h5d_close, did

    return, data
    
end



;+
;   Get the value of an attribute in a file.
;
;   :Private:
;
;   :Returns:
;       attribute value
;
;   :Params:
;       fileId : in, required, type=long
;           HDF 5 file identifier of the file to read
;       variable : in, required, type=string
;           path to attribute using '/' to navigate groups/datasets 
;           and '.' to indicate the attribute name
;
;   :Keywords:
;       error : out, optional, type=long
;           error value
;-

function wmb_h5_getdata_getattribute, file_id, objpath, attname

    compile_opt idl2, strictarrsubs

    obj_id = wmb_h5o_open(file_id, objpath)
    
    att_id = h5a_open_name(obj_id, attname)
    
    result = h5a_read(att_id)
    
    h5a_close, att_id

    wmb_h5o_close, obj_id

    return, result
  
end



;+
;   Pulls out a section of a HDF5 variable.
;
;   :Returns:
;       data array
;
;   :Params:
;       file_id : in, required, type=long
;           file_id of the HDF5 file
;       variable : in, required, type=string
;           variable name (with path if inside a group)
;
;   :Keywords:
;       bounds : in, optional, type=lonarr(3,ndims)
;           gives start value, end value, and stride for each dimension 
;           of the variable
;       error : out, optional, type=long
;           error value
;-

function wmb_h5_getdata, file_id, variable, bounds=bounds, error=error

    compile_opt idl2, strictarrsubs

    ; check the file id
    
    idtype = h5i_get_type(file_id)
    if idtype ne 'FILE' then message, 'Invalid HDF5 file ID'
    
    if N_elements(variable) eq 0 then message, 'No variable requested'

    ; check if we are writing an attribute
    
    dotPos = strpos(variable, '.', /reverse_search)

    if (dotPos eq -1) then begin
    
        ; read a variable

        result = wmb_h5_getdata_getvariable(file_id, $
                                            variable, $
                                            bounds=bounds, $
                                            empty=empty)

    endif else begin
    
        ; read an attribute
        
        objpath = strmid(variable, 0, dotPos)
        attname = strmid(variable, dotPos + 1L)
        
        result = wmb_h5_getdata_getattribute(file_id, objpath, attname)
        
    endelse

    return, result

end


;+
;   Simple example of wmb_h5_putdata usage.
;-

pro wmb_h5_getdata_example

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
    
    ; write the test data
    testdata = lonarr(10,50,100)
    testdata[5,*,*] = 5
    
    variablename = '/testgroup/3D_long_array'
    classname = variablename + '.CLASS'
    testclass = 'a big array'
    
    wmb_h5_putdata, fid, variablename, testdata
    wmb_h5_putdata, fid, classname, testclass
    
    ; full result is lonarr(10, 50, 100)
    fullResult = wmb_h5_getdata(fid, variablename)
    
    ; pull out a slice of the full result
    bounds = [[3, 3, 1], [5, 49, 2], [0, 49, 3]]
    res1 = wmb_h5_getdata(fid, variablename, bounds=bounds)
    help, res1
    
    ; compare indexing into fullResult versus slice pulled out
    same = array_equal(fullResult[3, 5:*:2, 0:49:3], res1)
    print, same ? 'equal' : 'error'
    
    ; grab an attribute
    print, wmb_h5_getdata(fid, classname)
    
    ; close
    h5f_close, fid

end

