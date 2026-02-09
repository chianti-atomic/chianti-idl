

PRO read_scups, splfile, splstr, add_empty=add_empty, verbose=verbose, $
                fits=fits, hdf5=hdf5, ascii=ascii, no_check=no_check

;+
; NAME:
;     READ_SCUPS
;
; PURPOSE:
;     Reads a file in the CHIANTI SCUPS format.
;
; CATEGORY:
;     CHIANTI; read; upsilons.
;
; CALLING SEQUENCE:
;     READ_SCUPS, SplFile, SplStr
;
; INPUTS:
;     SplFile:  The name of the SCUPS file to read.
;
; KEYWORDS:
;
;     ADD_EMPTY: If set, then an additional line of data is added to
;                SPLSTR.DATA. This extra line is empty, and is
;                intended for use by WRITE_SCUPS_GUI.
;     VERBOSE:   If set, then the routine prints out the time taken to
;                read the file.
;     FITS:      read the SCUPS data in gzipped FITS format
;     HDF5:      read the SCUPS data in HDF5 format
;     ASCII:     read the ASCII format of the SCUPS file.
;     NO_CHECK:  By default, the routine checks if the number of
;                temperatures given for the transition matches the actual
;                number of temperatures. Set this keyword to skip this
;                check.
;
; OUTPUTS:
;
;     SPLSTR   A structure containing the spline data. There are two
;              tags called 'info' and 'data' that are each structures.
;
;              UPSSTR.INFO has the following tags:
;
;              ION_NAME        STRING    'ca_2'
;              ION_Z           INT             20
;              ION_N           INT              2
;              ION_ROMAN       STRING    'Ca II'
;              ION_LATEX       STRING    '\ion{Ca}{ii}'
;              ION_LATEX_ALT   STRING    '\ion{Ca}{2}'
;              COMMENTS        STRING    Array[12]
;              CHIANTI_VER     STRING    '7.1'
;              TIME_STAMP      STRING    'Fri May 24 16:54:46 2013'
;              FILENAME        STRING    'ca_2.scups_auto'
;              MISSING         FLOAT          -1.00000
;              NTRANS          LONG               766
;
;              UPSSTR.DATA has the following tags:
;
;              LVL1            INT              1
;              LVL2            INT              2
;              T_TYPE          INT              2
;              DE              FLOAT          0.124400
;              GF              FLOAT           0.00000
;              LIM             FLOAT          -1.00000
;              C_UPS           FLOAT          0.560000
;              NSPL            INT             13
;              STEMP           FLOAT     Array[18]
;              SPL             FLOAT     Array[18]
;
; EXAMPLE:
;     IDL> zion2filename,26,13,fname
;     IDL> read_scups,fname+'.scups',scupstr
;     IDL> read_scups,fname+'.scups',scupstr,/ascii,/verbose
;     IDL> read_scups,fname+'.scups',scupstr,/hdf5,/verbose
;
; MODIFICATION HISTORY:
;     Ver. 1, 24-May-2013, Peter Young
;     Ver. 2, 28-May-2013, Peter Young
;         added /ADD_EMPTY keyword.
;     Ver. 3, 31-May-2013, Peter Young
;         the arrays i1,i2,i3 are now long integers.
;     Ver. 4, 28 Apr 2014, Giulio Del Zanna (GDZ)
;         rewritten to speed up by almost a factor of two.
;     Ver.5, 24-Jun-2020, Peter Young
;         fixed problem with formatting of the comment lines; updated
;         header.
;     Ver.6, 12 Oct 2020, GDZ
;         Added the option to read the FITS or HDF5 SCUPS files.
;         by default read the scups data in FITS format, gzipped.
;     Ver.7, 28-Oct-2020, Peter Young
;         Added /ascii keyword; modified how file options are
;         checked; if fits or hdf5 files don't exist, then
;         check for the ascii file.
;     Ver.8, 28-Oct-2020, Peter Young
;         Uses file_modtime to make sure that the modification time of
;         the hdf5 and fits files is later than the ascii file.
;     Ver.9, 23-Nov-2020, Peter Young
;         Replaced call to file_modtime with file_info (since
;         _modtime only introduced in IDL 8.5.1). Added some
;         additional information messages.
;     Ver.10, 04-Oct-2024, Peter Young
;         Added a check to make sure the number of temperatures/upsilons
;         for each transition matches the actual number. This adds about
;         2% to the read time.
;     Ver.11, 02-Dec-2025, Peter Young
;         Added /no_check to avoid doing the check introduced in v.10
;         (this is useful for checking problem files).
;-



  IF n_params() LT 1 THEN BEGIN
     print,'Use:  IDL>  read_scups, splfile, splstr [, /add_empty, /verbose, /ascii, /hdf5, /fits ]'
     return 
  ENDIF 

  t1=systime(/sec)

  missing_val=-1

  chck=file_info(splfile)
  IF chck.exists EQ 0 THEN BEGIN
     print,'% READ_SCUPS: the file was not found. Returning...'
     return
  ENDIF 
  mtime_ascii=chck.mtime
  
  splfile_fits=splfile+'.fits.gz'
  splfile_h5=splfile+'.h5'

  CASE 1 OF
     keyword_set(ascii): spl_opt=0
     keyword_set(hdf5): BEGIN
        spl_opt=2
        chck=file_info(splfile_h5)
        IF chck.exists EQ 0 THEN BEGIN
           print,'% READ_SCUPS: hdf5 file not found; reading the ascii file ('+file_basename(splfile)+').'
           spl_opt=0
        ENDIF ELSE BEGIN 
           mtime_h5=chck.mtime
           IF mtime_h5 LT mtime_ascii THEN BEGIN
              print,'% READ_SCUPS: hdf5 file needs updating so reading ascii file instead ('+file_basename(splfile)+').'
              spl_opt=0
           ENDIF
        ENDELSE 
     END 
     ELSE: BEGIN    ; default is FITS
        spl_opt=1
        chck=file_info(splfile_fits)
        IF chck.exists EQ 0 THEN BEGIN
           print,'% READ_SCUPS: fits file not found; reading the ascii file ('+file_basename(splfile)+').'
           spl_opt=0
        ENDIF ELSE BEGIN 
           mtime_fits=chck.mtime
           IF mtime_fits LT mtime_ascii THEN BEGIN
              print,'% READ_SCUPS: fits file needs updating so reading ascii file instead ('+file_basename(splfile)+').'
              spl_opt=0
           ENDIF
        ENDELSE 
     END 
  ENDCASE 


 ;
 ; This case statement reads the FITS and HDF5 files directly into the
 ; output structure and exits.
 ;
  CASE spl_opt OF
     1: BEGIN
        info=mrdfits(splfile_fits,1,/silent)
        data=mrdfits(splfile_fits,2,/silent)
       ;
        splstr={info:info, data:data}
       ;
        if keyword_set(verbose) then begin 
           t2=systime(/sec)
           print,'% READ_SCUPS: the scups FITS file was read.'
           print,'% READ_SCUPS: '+splfile+', seconds:', string(format='(f8.3)',t2-t1)
        endif 
        return
     END
     2: BEGIN
        H5_OPEN
        st = H5_PARSE(splfile_h5, /READ_DATA)
    ; the following is necessary to clear the memory:
        H5_CLOSE
        splstr={info:st.scups._data.info, data:st.scups._data.data}
        delvarx,st
        if keyword_set(verbose) then begin 
           print,'% READ_SCUPS: the HDF5 scups file was read.'
           t2=systime(/sec)
           print,'% READ_SCUPS: '+splfile+', seconds:', string(format='(f8.3)',t2-t1)
        endif 
        return
     END
     ELSE:
  ENDCASE 


;
; This routine returns the number of lines in the ascii elvlc
; file. This is useful for speeding up the routine (see later). 
;

     result=query_ascii(splfile,info)
     nlines=info.lines

;
; Read the entire file into a string array
;
     file_string=strarr(nlines)
     openr,lin,splfile,/get_lun
     readf,lin,file_string
     free_lun,lin


; By working out where the -1 is, the number of data lines can be found.
;

; now read the comments at the end. Start from the end going
; backwards. Assume that there are only two lines with '-1'.
     lines_comment=0
     comment=0
     ref=''

     WHILE comment lt 2 DO  BEGIN
        if trim(file_string[nlines-1-lines_comment]) EQ '-1' then $
           comment=comment+1 else ref=[ref,file_string[nlines-1-lines_comment]]
        lines_comment=lines_comment+1
     ENDWHILE 

;
; PRY: 24-Jun-2020
;
     ref=ref[1:*]
     ref=reverse(ref)

     ndata=nlines-lines_comment



     ntrans=ndata/3
     IF 3*ntrans NE ndata THEN BEGIN
        print,'%READ_SCUPS: Each transition should have 3 lines of data!!  Returning...'
        return
     ENDIF 

;
; Each transition has 3 lines of data, so put these into separate
; arrays (data1, data2, data3).
;
     i1=lindgen(ntrans)*3
     i2=lindgen(ntrans)*3+1
     i3=lindgen(ntrans)*3+2
;

     bigarr=fltarr(8, ntrans)
     reads, file_string[i1], bigarr

; number of temperatures for each transition:  
     nspl=reform(bigarr[5,*]) 


; nt_max (maximum number of temperatures for a transition) sets the
; size of the temperature and upsilon arrays 
;
     nt_max=max(nspl)

;
; Create the output data structure.
;
     str2={lvl1: 0, lvl2: 0, t_type: 0, de: 0., gf: 0., lim: 0., c_ups: 0., nspl: 0, $
           stemp: fltarr(nt_max)-1, spl: fltarr(nt_max)-1}

     IF keyword_set(add_empty) THEN BEGIN
        datastr=replicate(str2,ntrans+1)
     ENDIF ELSE BEGIN
        datastr=replicate(str2,ntrans)
     ENDELSE 
;
     datastr[0:ntrans-1].lvl1=reform(bigarr[0,*])
     datastr[0:ntrans-1].lvl2=reform(bigarr[1,*])
     datastr[0:ntrans-1].de=reform(bigarr[2,*])
     datastr[0:ntrans-1].gf=reform(bigarr[3,*])
     datastr[0:ntrans-1].lim=reform(bigarr[4,*])
     datastr[0:ntrans-1].nspl=reform(bigarr[5,*]) 
     datastr[0:ntrans-1].t_type=reform(bigarr[6,*])
     datastr[0:ntrans-1].c_ups=reform(bigarr[7,*])
;

     if min(nspl) eq   nt_max then num_temp=nt_max else  begin 
; find the unique numbers of temperatures
        
        num_temp=-1
        
        for ii=1,nt_max do begin  
           indt=where(nspl eq ii, nnt)
           if nnt gt 0 then num_temp=[num_temp,ii]
        endfor 
        num_temp=num_temp[1:*]
        
     endelse 

     nvt=n_elements(num_temp)

; 
; Fill in each temperature array at a time
;
; PRY, 4-Oct-2024
;  I've added a check on the string lengths to make sure num_temp is consistent
;  with nspl.
;  Use /no_check to skip this check.
;
     IF NOT keyword_set(no_check) THEN BEGIN 
       for ii=0,nvt-1 do begin   
        
         ind=where(nspl eq  num_temp[ii], nn)
         bigarr=fltarr( num_temp[ii], nn)-1
        
         reads,file_string[i2[ind]], bigarr
         len=strlen(file_string[i2[ind]])
         k=where(len NE num_temp[ii]*12,nk)
         IF nk GT 0 THEN BEGIN
           print,file_basename(splfile),' Temperature array does not match no. of temperatures! Returning...'
           return
         ENDIF 
         datastr[ind].stemp[0:num_temp[ii]-1]=bigarr

         reads,file_string[i3[ind]], bigarr
         len=strlen(file_string[i2[ind]])
         k=where(len NE num_temp[ii]*12,nk)
         IF nk GT 0 THEN BEGIN
           print,file_basename(splfile),' Upsilon array does not match no. of upsilons! Returning...'
           return
         ENDIF 
         datastr[ind].spl[0:num_temp[ii]-1]=bigarr
        
       ENDFOR 
     ENDIF 


; Get the CHIANTI version number
;
     vnum=ch_get_version()

;
; Create the various ion identifiers from the filename.
;
     basename=file_basename(splfile)
     bits=str_sep(basename,'.')
     ion_name=bits[0]
     convertname,ion_name,ion_z,ion_n
     ion_z=fix(ion_z)
     zion2spectroscopic,ion_z,ion_n,ion_roman
     ion_roman=strcompress(ion_roman)
     bits=strsplit(ion_roman,' ',/extract)
     ion_latex='\ion{'+bits[0]+'}{'+strlowcase(bits[1])+'}'
     ion_latex_alt='\ion{'+bits[0]+'}{'+trim(ion_n)+'}'


     IF keyword_set(add_empty) THEN ntrans=ntrans+1

;
; Create the info structure
;
     info={ion_name: ion_name, ion_z: ion_z, ion_n: ion_n, ion_roman: ion_roman, $
           ion_latex: ion_latex, ion_latex_alt: ion_latex_alt, comments: ref, $
           chianti_ver: vnum, time_stamp: systime(), $
           filename: splfile, missing: float(missing_val), $
           ntrans: ntrans}


;
; Combine the info and data structures into the final output
; structure. 
;
     splstr={info: info, data: datastr}


  if keyword_set(verbose) then begin 
     print,'% READ_SCUPS: the ASCII scups file was read.'
     t2=systime(/sec)
     print,'% READ_SCUPS: '+splfile+', seconds:', string(format='(f8.3)',t2-t1)
  endif 

END

