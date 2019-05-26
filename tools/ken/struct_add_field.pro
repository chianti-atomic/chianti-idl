pro struct_add_field, struct, tag, data, before=before, after=after $
                    , itag=itag
;Add a new field to an existing structure. If the field already exists,
; the data will be copied into the existing field. Keyword may be used
; to control placement of a new field in the structure.
;
;Inputs:
; tag (string) Case insensitive name for the new structure field. Leading
;  and trailing spaces will be ignored. If the requested field already
;  exists, the specified data are copied into the existing field.
; data (any) Data used to fill the new or existing field specified in tag.
;  The dimensionality of data must agree with the dimensionality of struct,
;  except that scalar data will be replicated into a vector structure, if
;  necessary.
; [before=] (string) The new structure field will be created immediately
;  before the first occurence of an existing structure field with the tag
;  name specified by this keyword. An error occurs if the specified tag
;  name does not already exist.
; [after=] (string) The new structure field will be created immediately
;  after the first occurence of an existing structure field with the tag
;  name specified by this keyword. An error occurs if the specified tag
;  name does not already exist.
; [itag=] (scalar) The new structure field will be created in position
;  itag, where itag=0 indicates that the new field should be created at
;  the beginning of the structure, itag=1 indicates that the new field
;  should be the second field in the structure, etc.
;
;Input/Output:
; struct (structure) structure to be modified.
;
;Examples:
;
; Add wavelength vector to structure:
;
;   IDL> struct_add_field, sme, 'wave', wave
;
; Add integer mask vector to structure, setting each entry to 1.
;
;   IDL> struct_add_field, sme, 'mask', 1
;
;History:
; 2003-Aug-11 Valenti  Adapted from struct_delete_field.pro.
; 2003-Oct-23 Valenti  Allow new field to be a structure with no nesting.

if n_params() lt 3 then begin
  print, 'syntax: struct_add_field, struct, tag, data [begin=, after=, itag=]'
  return
endif

;Check that input is a structure.
  if size(struct, /tname) ne 'STRUCT' then begin
    message, 'first argument is not a structure'
  endif

;Check that no more than one of before=, after=, and itag= are set.
  if keyword_set(before) + keyword_set(after) $
                         + n_elements(itag) gt 1 then begin
    message, 'specify no more than one of before=, after=, and itag='
  endif

;Check dimensionality of data.
  nstruct = n_elements(struct)
  ndata = n_elements(data)
  if ndata ne 1 and ndata ne nstruct then begin
    message, 'dimensionality of struct and data are incompatible'
  endif

;Get list of structure tags.
  tags = tag_names(struct)
  ntags = n_elements(tags)

;Check whether the requested field exists in input structure.
  ctag = strupcase(strtrim(tag, 2))		;canoncial form of tag
  imatch = where(tags eq ctag, nmatch)
  if nmatch gt 0 then begin
    struct.(imatch[0]) = data			;overwrite data
    return
  endif

;Figure out where to place the new tag.
  if n_elements(itag) eq 0 then itag = ntags	;end of structure is default
  if keyword_set(before) then begin
    iwhr = where(tags eq before, nwhr)
    if nwhr eq 0 then message, 'before=' + before + ' tag not in structure'
    itag = iwhr[0]
  endif
  if keyword_set(after) then begin
    iwhr = where(tags eq after, nwhr)
    if nwhr eq 0 then message, 'after=' + after + ' tag not in structure'
    itag = iwhr[0] + 1
  endif

;Use first record as a template for new structure fields.
  rec = struct[0]

;Copy any fields that precede target for new field. Append new field.
  if itag gt 0 then begin			;target field occurs first
    new = create_struct(tags[0], rec.(0))	;initialize structure
    for i=1, itag-1 do begin			;insert leading unchange
      new = create_struct(new, tags[i], rec.(i))
    endfor
    new = create_struct(new, tag, data[0])	;add new field
  endif else begin
    new = create_struct(tag, data[0])		;start with new field
  endelse

;Replicate remainder of structure after desired tag.
  for i=itag, ntags-1 do begin
    new = create_struct(new, tags[i], rec.(i))
  endfor

;Handle trivial case when structure array contains only one element.
  if nstruct eq 1 then begin
    struct = new
    return
  endif

;Create vector structure to match original dimensionality.
  new = replicate(new, nstruct)

;Copy data into vector structure.
  for i=0, itag-1 do new.(i) = struct.(i)
  if size(data, /tname) eq 'STRUCT' then begin
    for j=0, n_tags(data)-1 do new.(i).(j) = data.(j)
  endif else begin
    new.(i) = data
  endelse
  for i=itag, ntags-1 do new.(i+1) = struct.(i)
  struct = new

end
