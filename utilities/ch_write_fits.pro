;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;       CH_WRITE_FITS
;
; PURPOSE:
;       Write  standard FITS binary table data from CHIANTI input  structure.
;
; CALLING SEQUENCE:
;       CH_WRITE_FITS, Input, Filename
;
; INPUTS:
;
;       Input = Structure to be written to FITS file.
;
;
; OUTPUTS:
;
;
;       Filename = String containing the name of the file to be written.
;                CH_WRITE_FITS creates two binary table extension in a single
;                FITS file. The second one is appended as a new extension.
;
; OPTIONAL INPUTS:  Header COMMENTS.
;
; KEYWORDS: 
;           head1, head2
;          Additional  COMMENTS to be added at the bottom of the two binary tables.
;
;
; NOTES:
; 
;       Any existing FITS file can be over-written or not.
;       Use CH_READ_FITS  to convert the FITS file back into a structure. 
;
; CALLS:
;
;	FXPAR(), FXADDPAR, IS_IEEE_BIG(), HOST_TO_IEEE, DIALOG_MESSAGE
;
; COMMON BLOCKS: none.
;
; RESTRICTIONS:
;       (1)     Limited to 127 columns in tables by IDL structure limits.
;       (2)     String columns with all columns of zero length crash the
;               program
;       (3)     The input structure has to be of the type TRANSITIONS.
;
; PREV. HIST. :
;
;       The subroutines in this procedure are extracted without modifications from
;       the MWRFITS.PRO routine, written by T. McGlynn Version 0.95 2000-11-06
;       and present in the ASTRON library (in SolarSoft under /gen/idl_libs/astron/).
;
; EXAMPLE:
;
;         ch_write_fits, transitions , 'test.fits'
;
; WRITTEN     : 
;       Ver.1, 22-May-02 Giulio Del Zanna (GDZ)
;
;
; MODIfICATION HISTORY:
;
; VERSION     : 1, 22-May-02, GDZ
;
;-
pro chk_and_upd, header, key, value, comment

; Add a keyword as non-destructively as possible to a FITS header

xcomm = ""
if n_elements(comment) gt 0 then xcomm = comment
if n_elements(header) eq 0 then begin
   
   fxaddpar, header, key, value, xcomm
   
endif else begin
   
   oldvalue = fxpar(header, key, count=count, comment=oldcomment)
   
   if (count eq 1) then begin

      qchange = 0               ; Set to 1 if either the type of variable or its
                                ; value changes.
      size1 = size(oldvalue) & size2 = size(value)
      if size1[size1[0]+1] NE size2[size2[0]+1] then qchange = 1 $
      else if (oldvalue ne value) then qchange = 1

      if (qchange) then begin

         if n_elements(oldcomment) gt 0 then xcomm = oldcomment[0]
         fxaddpar, header, key, value, xcomm
         
      endif
      
   endif else begin
      
      fxaddpar, header, key, value, xcomm
   endelse
   
endelse
end


pro mwr_tablehdr, lun, input, header,    $
       no_types=no_types,                $
       logical_cols = logical_cols,	  $
       bit_cols = bit_cols,		  $
       nbit_cols= nbit_cols,             $
       no_comment=no_comment,            $
       silent=silent

;  Create and write the header for a binary table.

if not keyword_set(no_types) then no_types = 0
nfld = n_tags(input[0])
if nfld le 0 then begin
   print, 'CH_WRITE_FITS Error: Input contains no structure fields.'
   return
endif

tags = tag_names(input)

; Get the number of rows in the table.

nrow = n_elements(input)

dims = lonarr(nfld)
tdims = strarr(nfld)
types = strarr(nfld)

;
; Get the type and length of each column.  We do this
; by examining the contents of the first row of the structure.
;

nbyte = 0

for i=0, nfld-1 do begin

   a = input[0].(i)

   sz = size(a)
   
   nelem = sz[sz[0]+2]
   type_ele = sz[sz[0]+1]
   if type_ele eq 7 then begin
      maxstr = max(strlen(input.(i)))
   endif
   
   dims[i] = nelem
   
   if (sz[0] lt 1) or (sz[0] eq 1 and type_ele ne 7) then begin
      tdims[i] = ''
   endif else begin
      tdims[i] = '('
      
      if type_ele eq 7 then begin
         tdims[i] = tdims[i] + strcompress(string(maxstr), /remo) + ','
      endif
      
      for j=1, sz[0] do begin
         tdims[i] = tdims[i] + strcompress(sz[j])
         if j ne sz[0] then tdims[i] = tdims[i] + ','
      endfor
      tdims[i] = tdims[i] + ')'
   endelse
   
   
   case type_ele of
      1: 	begin
         types[i] = 'B'
         nbyte = nbyte + nelem
      end
      2:	begin
         types[i] = 'I'
         nbyte = nbyte + 2*nelem
      end
      3:	begin
         types[i] = 'J'
         nbyte = nbyte + 4*nelem
      end
      4:	begin
         types[i] = 'E'
         nbyte = nbyte + 4*nelem
      end
      5:	begin
         types[i] = 'D'
         nbyte = nbyte + 8*nelem
      end
      6:	begin
         types[i] = 'C'
         nbyte = nbyte + 8*nelem
      end
      7:	begin
         types[i] = 'A'
         nbyte = nbyte + maxstr*nelem
         dims[i] = maxstr*nelem
      end
      9:   begin
         types[i] = 'M'
         nbyte = nbyte + 16*nelem
      end
      0:   begin
         print,'CH_WRITE_FITS Error: Undefined structure element??'
         return
      end
      8:   begin
         print, 'CH_WRITE_FITS Error: Nested structures'
         return
      end
      else:begin
         print, 'CH_WRITE_FITS Error: Cannot parse structure'
         return
      end
   endcase
endfor

; Put in the required FITS keywords.
chk_and_upd, header, 'XTENSION', 'BINTABLE', 'Binary table written by CH_WRITE_FITS'
chk_and_upd, header, 'BITPIX', 8, 'Required value'
chk_and_upd, header, 'NAXIS', 2, 'Required value'
chk_and_upd, header, 'NAXIS1', nbyte, 'Number of bytes per row'
chk_and_upd, header, 'NAXIS2', n_elements(input), 'Number of rows'
chk_and_upd, header, 'PCOUNT', 0, 'Normally 0 (no varying arrays)'
chk_and_upd, header, 'GCOUNT', 1, 'Required value'
chk_and_upd, header, 'TFIELDS', nfld, 'Number of columns in table'

if (not keyword_set(no_comment)) then begin
   fxaddpar, header, 'COMMENT', ' ', after='TFIELDS'
   fxaddpar, header, 'COMMENT', ' *** End of required fields ***', after='TFIELDS'
   fxaddpar, header, 'COMMENT', ' ', after='TFIELDS'
endif
;
; Handle the special cases.
;
if keyword_set(logical_cols) then begin
   nl = n_elements(logical_cols)
   for i = 0, nl-1 do begin
      icol = logical_cols[i]
      if types[icol-1] ne 'A'  then begin
         print,'WARNING: Invalid attempt to create Logical column:',icol
         goto, next_logical
      endif
      types[icol-1] = 'L'
next_logical:
   endfor
endif

if keyword_set(bit_cols) then begin
   nb = n_elements(bit_cols)
   if nb ne n_elements(nbit_cols) then begin
      print,'WARNING: Bit_cols and Nbit_cols not same size'
      print,'         No bit columns generated.'
      goto, after_bits
   endif
   for i = 0, nb-1 do begin
      nbyte = (nbit_cols[i]+7)/8
      icol = bit_cols[i]
      if types[icol-1] ne 'B'  or (dims[icol-1] ne nbyte) then begin
         print,'WARNING: Invalid attempt to create bit column:',icol
         goto, next_bit
      endif
      types[icol-1] = 'X'
      tdims[icol-1] = ''
      dims[icol-1] = nbit_cols[i]
next_bit:
   endfor
after_bits:
endif

; First add in the TTYPE keywords if desired.
;
if not no_types then begin
   for i=0, nfld - 1 do begin
      key = 'TTYPE'+strcompress(string(i+1),/remove)
      if not keyword_set(use_colnums) then begin
         value= tags[i]+' '
      endif else begin
         value = 'C'+strmid(key,5,2)
      endelse
      chk_and_upd, header, key, value
   endfor
   
   if (not keyword_set(no_comment)) then begin
      fxaddpar, header, 'COMMENT', ' ', before='TTYPE1'
      fxaddpar, header, 'COMMENT', ' *** Column Names *** ',before='TTYPE1'
      fxaddpar, header, 'COMMENT', ' ',before='TTYPE1'
   endif
endif
; Now add in the TFORM keywords
for i=0, nfld-1 do begin
   if dims[i] eq 1 then begin
      form = types[i]
   endif else begin
      form=strcompress(string(dims[i]),/remove) + types[i]
   endelse
   
   tfld = 'TFORM'+strcompress(string(i+1),/remove)
   
; Check to see if there is an existing value for this keyword.
; If it has the proper value we will not modify it.
; This can matter if there is optional information coded
; beyond required TFORM information.
   
   oval = fxpar(header, tfld)
   oval = strcompress(string(oval),/remove_all)
   if (oval eq '0')  or  (strmid(oval, 0, strlen(form)) ne form) then begin
      chk_and_upd, header, tfld, form
   endif
endfor

if (not keyword_set(no_comment)) then begin
   fxaddpar, header, 'COMMENT', ' ', before='TFORM1'
   fxaddpar, header, 'COMMENT', ' *** Column formats ***', before='TFORM1'
   fxaddpar, header, 'COMMENT', ' ', before='TFORM1'
endif

; Now write TDIM info as needed.
firsttdim = -1
for i=0, nfld-1 do begin
   if tdims[i] ne '' then begin
      chk_and_upd, header, 'TDIM'+strcompress(string(i+1),/remo), tdims[i]
   endif
   if firsttdim eq -1 then firsttdim = i
endfor

w=where(tdims ne '')
if w[0] ne -1 and not keyword_set(no_comment) then begin
   fxaddpar, header, 'COMMENT', ' ',   $
     before='TDIM'+strcompress(string(firsttdim+1),/remo)
   fxaddpar, header, 'COMMENT', ' *** Column dimensions (2 D or greater) ***',  $
     before='TDIM'+strcompress(string(firsttdim+1),/remo)
   fxaddpar, header, 'COMMENT', ' ', $
     before='TDIM'+strcompress(string(firsttdim+1),/remo)
endif
; Write to the output device.
mwr_header, lun, header

end

pro mwr_tabledat, lun, input, header
;
;  Write binary table data to a FITS file.
;
; file		-- unit to which data is to be written.
; Input		-- IDL structure
; Header	-- Filled header

nfld = n_tags(input)

; Pad out strings to constant length.
for i=0, nfld-1 do begin
   
   sz = size(input.(i))
   nsz = n_elements(sz)
   typ = sz[nsz-2]
   if (typ eq 7) then begin

      siz = max(strlen(input.(i)))
      
      blanks = string(bytarr(siz) + 32b)
      input.(i) = strmid(input.(i)+blanks, 0, siz)

   endif
endfor

; Use Astron library routine to convert to IEEE (since byteorder
; may be buggy).
if not is_ieee_big() then host_to_ieee, input

nbyte = long(fxpar(header, 'NAXIS1'))
nrow = long(fxpar(header, 'NAXIS2'))

siz = nbyte*nrow

padding = 2880 - (siz mod 2880)
if padding eq 2880 then padding = 0

;
; Write the data segment.
;
writeu, lun, input

; If necessary write the padding.
;
if padding gt 0 then begin
   pad = bytarr(padding)        ; Should be null-filled by default.
   writeu, lun, pad
endif

end



pro mwr_header, lun, header
;
; Write a header
;

; Fill strings to at least 80 characters and then truncate.

space = string(replicate(32b, 80))
header = strmid(header+space, 0, 80)

w = where(strmid(header,0,8) eq "END     ")

if w[0] eq -1 then begin

   header = [header, strmid("END"+space,0,80)]
   
endif else begin
   if (n_elements(w) gt 1) then begin 
                                ; Get rid of extra end keywords;
      print,"CH_WRITE_FITS Warning: multiple END keywords found."
      for irec=0L, n_elements(w)-2 do begin
         header[w[irec]] = strmid('COMMENT INVALID END REPLACED'+  $
                                  space, 0, 80)
      endfor
   endif

                                ; Truncate header array at END keyword.
   header = header[0:w[n_elements(w)-1]]
endelse

nrec = n_elements(header)
if nrec mod 36 ne 0 then header = [header, replicate(space,36 - nrec mod 36)]

writeu, lun, byte(header)
end


pro mwr_dummy, lun

; Write a dummy primary header-data unit.

fxaddpar, header, 'SIMPLE', 'T','Dummy primary header created by CH_WRITE_FITS on ' 
fxaddpar, header, 'BITPIX', 8, ' on '+systime() 
fxaddpar, header, 'NAXIS', 0, 'No data is associated with this header'
fxaddpar, header, 'EXTEND', 'T', 'Two extensions are present'
fxaddpar, header, 'COMMENT', ' '
fxaddpar, header, 'COMMENT', ' This FITS file contains CHIANTI atomic data '
fxaddpar, header, 'COMMENT', ' '
fxaddpar, header, 'COMMENT', ' CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of'
fxaddpar, header, 'COMMENT', ' Astrophysical Plasmas. It is a collaborative project involving the Naval'
fxaddpar, header, 'COMMENT', ' Research Laboratory (USA), the University of Florence (Italy), the'
fxaddpar, header, 'COMMENT', ' University of Cambridge and the Rutherford Appleton Laboratory (UK)'
fxaddpar, header, 'COMMENT', ' '


mwr_header, lun, header
end



PRO  ch_write_fits, xinput, file,  head1=head1, head2=head2

on_ioerror, open_error

; Check required keywords.


if n_params() NE 2  THEN  BEGIN 
   print, 'CH_WRITE_FITS: Usage:'
   print, '    CH_WRITE_FITS, struct_name, file, [head=head]'
   return 
endif


; Save the data into an array/structure that we can modify.

;check that it is a 'TRANSITIONS' structure at least by having the minimal tags:


test = ch_check_str (xinput, error=error)

IF NOT  test THEN  message, 'Wrong input structure given!!!!'

input1 = rem_tag(xinput,'lines')
input2 = xinput.lines


found= findfile(file, count=count)

;check if it exists:

IF count EQ 1 THEN BEGIN 

   Result = DIALOG_MESSAGE( 'The file '+file+' already exists.'+$
                            ' Do you want to overwrite this file ? ', /question)

   IF result EQ 'No' THEN return

ENDIF 


; Open the input file.

if !version.os eq 'vms' then openw, lun, file, 2880, /block, /none, /get_lun $
else openw, lun, file, /get_lun


mwr_dummy, lun


; Create a binary table.


header = ["COMMENT  ", $
          "COMMENT This FITS file contains CHIANTI atomic data ", $
          "COMMENT   ", $
          "COMMENT  For each spectral line, the following information is included:", $
          "COMMENT   ", $
          "COMMENT iz           The atomic number of the elements (e.g., 26=Fe) ",$
          "COMMENT ion          The ionisation stage (e.g., 13=XIII) ",$
          "COMMENT snote        The identification of the ion (e.g., 'Fe XXIV d') ",$
          "COMMENT ident        The identification of the transition, ",$
          "COMMENT                configuration and terms in text form. ",$
          "COMMENT ident_latex  The identification of the transition, ",$
          "COMMENT                configuration and terms in latex form. ",$
          "COMMENT lvl1   The lower level of the transition (see .elvlc file for ion) ",$
          "COMMENT lvl2   The upper level for transition. ",$
          "COMMENT tmax   The temperature of maximum emission of the line,",$
          "COMMENT               zero in the isothermal approximation",$
          "COMMENT wvl    Wavelength of the transition, in Angstroms. ",$
          "COMMENT flag   A flag, (-1) if the line has only theoretical energy levels.",$
          "COMMENT               Otherwise flag=0. ",$
          "COMMENT int    Intensity of line (erg/cm2/s/sr or phot/cm2/s/sr), ",$
          "COMMENT               without the  element abundance factor.  ",$
          "COMMENT goft   The contribution function G(T) of the line (optional). ",$
          "COMMENT   "]

IF n_elements(head1) NE 0 THEN header = [header, head1]

header = [header, "COMMENT   ", $
          "COMMENT  Created by CH_WRITE_FITS on  "+systime(), $
          "COMMENT   "]


mwr_tablehdr, lun, input2, header,    $
  no_types=no_types,                $
  logical_cols = logical_cols,	  $
  bit_cols = bit_cols,		  $
  nbit_cols= nbit_cols,             $
  no_comment=no_comment 

mwr_tabledat, lun, input2, header


free_lun, lun

;APPEND THE SECOND PART 

if !version.os eq 'vms' then openu, lun, file, 2880, /block, /none, /get_lun, /append $
else openu, lun, file, /get_lun, /append

header = ["COMMENT  ", $
          "COMMENT This FITS file contains CHIANTI atomic data  ", $
          "COMMENT   ", $
          "COMMENT  The following information is included: ", $
          "COMMENT   ", $
          "COMMENT  ioneq_name     The ion balance file used (full path). ",$
;          "COMMENT  ioneq          The ion balance values. ",$
          "COMMENT  ioneq_logt        The Log10 T values associated. ",$
          "COMMENT  ioneq_ref      The references. ",$
          "COMMENT   --------------------------------------------------------------- ", $
          "COMMENT     IF AN ISOTHERMAL APPROXIMATION IS USED:", $
          "COMMENT   ", $
          "COMMENT  logt_isothermal     The Log10(T) values used.  ",$
          "COMMENT  logem_isothermal    The Log10(EM) values used.  ",$
          "COMMENT   ", $
          "COMMENT     OTHERWISE: ", $
          "COMMENT  dem_name      The differential emission measure file used (full path). ",$
          "COMMENT  dem           The Log10 DEM values  ",$
          "COMMENT  dem_logt         The Log10 T values associated. ",$
          "COMMENT  dem_ref       The references. ",$
          "COMMENT   --------------------------------------------------------------- ", $
          "COMMENT  model_file    A string indicating the (Te,Ne) model file used  ",$
          "COMMENT  model_name    A string indicating the  model used  ",$
          "COMMENT              (e.g.  constant density or constant pressure) ",$
          "COMMENT  model_ne   The Ne value used. ",$
          "COMMENT  model_pe   The Pe value used. ",$
          "COMMENT  model_te   The Te value used. ",$
          "COMMENT  wvl_units     The wavelength units. ",$
          "COMMENT  wvl_limits       The wavelength limits specified by the user. ",$
          "COMMENT  int_units     The intensity  (or G(T) units)  ",$
          "COMMENT  date         The date and time when the structure was created. ",$
          "COMMENT  version      The version number of the CHIANTI database used. ",$
          "COMMENT  add_protons   A flag indicating if proton data were used (1) or not (0) ",$
          "COMMENT  photoexcitation  A flag  to indicate if photoexcitation was ",$
          "COMMENT                     included (1) or not (0). ",$
          "COMMENT  radtemp      The blackbody radiation field temperature used ",$
          "COMMENT                  (if photoexcitation was included). ",$
          "COMMENT  rphot        Distance from the centre of the star in stellar ",$
          "COMMENT                  radius units (if photoexcitation was included). ",$
          "COMMENT   "]

IF n_elements(head2) NE 0 THEN header = [header, head2]

header = [header, "COMMENT   ", $
          "COMMENT  Created by CH_WRITE_FITS on  "+systime(), $
          "COMMENT   "]


; Create a binary table.

mwr_tablehdr, lun, input1, header

;,    $
;  no_types=no_types,                $
;  logical_cols = logical_cols,	  $
;  bit_cols = bit_cols,		  $
;  nbit_cols= nbit_cols,             $
;  no_comment=no_comment 

mwr_tabledat, lun, input1, header


free_lun, lun

return  


; Handle error in opening file.
open_error:
on_ioerror, null
print, ' Error: Cannot open output: ', file
if n_elements(lun) gt 0 then free_lun, lun

return
END 

