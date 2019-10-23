
; wmb_h5_list
;
; A modified version of the IDL H5_LIST function.
; 



pro _wmb_h5_list_parse, struc, out
  compile_opt idl2, hidden
  
  case struc._type of
    'GROUP' : begin
      path = struc._path+'/'+struc._name
      if (STRMID(path, 0, 2) eq '//') then begin
        path = STRMID(path, 1)
      endif
      out = [[out],['group', path, '']]
      tags = TAG_NAMES(struc)
      for i=0,N_ELEMENTS(tags)-1 do begin
        if (SIZE(struc.(i), /TNAME) eq 'STRUCT') then begin
          _wmb_h5_list_parse, struc.(i), out
        endif
      endfor
    end
    'DATASET' : begin
      path = struc._path+'/'+struc._name
      if (STRMID(path, 0, 2) eq '//') then begin
        path = STRMID(path, 1)
      endif
      dims = '['
      for i=0,struc._ndimensions-1 do begin
        dims += STRTRIM(struc._dimensions[i],2)
        if (i ne struc._ndimensions-1) then begin
          dims += ', '
        endif
      endfor
      dims += ']'
      out = [[out],['dataset', path, struc._datatype+' '+dims]]
    end
    else :
  endcase
  
end

pro wmb_h5_list, file, FILTER=filterIn, OUTPUT=out
  compile_opt idl2, hidden
  on_error, 2

  if ((N_ELEMENTS(file) ne 1) || (SIZE(file, /TNAME) ne 'STRING')) then begin
    message, 'FILE must be a scalar string'
    return
  endif
  if (N_ELEMENTS(filterIn) gt 1) then begin
    message, 'FILTER must be a scalar string'
    return
  endif
  if (~(FILE_INFO(file)).exists) then begin
    message, 'File not found: '+file
    return
  endif
  if (~H5F_IS_HDF5(file)) then begin
    message, 'FILE is not an HDF5 file: '+file
    return
  endif
  
  catch, err
  if (err ne 0) then begin
    catch, /cancel
    message, 'Unable to read file: '+file
  endif
  
  struc = H5_PARSE(file)
  hasOut = ARG_PRESENT(out)
  
  out = !NULL
  ; Parse the structure and build up the out string array
  _wmb_h5_list_parse, struc, out
  if (N_ELEMENTS(out) eq 0) then return
  
  ; WMB: Added the following to handle empty HDF5 files
  if (size(out, /N_DIMENSIONS) ne 2) then return
  
  ; Change the first line to show the name of the file and save it off
  fileline = out[*,0]
  fileline[0] = 'file'
  fileline[1] = STRMID(fileline[1], 1)
  out = out[*,1:*]
  
  ; Filter list
  if (N_ELEMENTS(filterIn) eq 1) then begin
    ; Break into parts based on '*'
    filter = STRSPLIT(filterIn, '*', /EXTRACT)
    keep = !NULL
    for i=0,N_ELEMENTS(out[0,*])-1 do begin
      !NULL = where(STRUPCASE(out[0,i]) eq STRUPCASE(filter), cnt)
      if (cnt ne 0) then begin
        keep = [[keep],[out[*,i]]]
        continue
      endif
      keepFlag = 1b
      curpos = 0
      for j=0,N_ELEMENTS(filter)-1 do begin
        curpos = STRPOS(out[1,i], filter[j], curpos)
        if (curpos eq -1) then begin
          keepFlag = 0b
          break
        endif
      endfor
      if (keepFlag) then begin
        keep = [[keep],[out[*,i]]]
        continue
      endif
    endfor
    out = keep
  endif
  
  ; Add file line back to output
  out = [[fileline], [out]]

  ; Pretty print
  if (~hasOut && (N_ELEMENTS(out) ne 0)) then begin
    space = 5
    mx = N_ELEMENTS(out) eq 3 ? STRLEN(out) : MAX(STRLEN(out), DIMENSION=2)
    fmt = "(A-"+STRTRIM(mx[0]+space,2)+",A-"+STRTRIM(mx[1]+space,2)+",A0)"
    print, out, format=fmt
  endif
  
end