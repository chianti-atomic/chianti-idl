;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), 
;       Cambridge University (United Kingdom), George Mason University (USA), and
;       the University of Michigan (USA).
;
;
; NAME:
;   READ_CIRECLVL_STR
;
; PURPOSE:
;
;   read cireclvl files - these contain data that provide rates for radiative recombination 
;       onto bound levels (.reclvl) or rates for the collisional ionization into bound levels
;       (.cilvl)
;
;
; CATEGORY:
;
;   science
;
; CALLING SEQUENCE:
;
;       READ_CIRECLVL, gname, 
;
;
; INPUTS:
;
;   Gname:  the name of the ion , i.e. 'c_4'
;   Ci:  if set, the gname.cilvl file will be read
;   Rec:  if set, the gname.reclvl file will be read
;
;
; OUTPUTS:
;   cireclvlStr, a structure with the following tags
;
;       filename:  the full file name
;       lvl1:  1D array of indices of the lower level (starting at 1)
;       lvl2:  1D array of indices of the upper level (starting at 1)
;       temp:  2D array of temperatures
;       rate;  2D array of rate coefficients
;       ntemp;  1D array of the number of temperatures for each rate
; and
;   Ref:   1D string array of references to the data in the scientific literature
;
;
;
; EXAMPLE:
;
;             > read_cireclvl_str, cireclvlstr, ref
;             
;
; MODIFICATION HISTORY:
;   Written by: Ken Dere (GMU)
;   November 2010:     Version 1.0
;
;
;-
pro read_cireclvl, gname, cireclvlStr, ci=ci, rec=rec
; 
;
    if n_params(0) lt 2 then begin
        print,''
        print,' type> read_cireclvl, cireclvlstr, ref, /ci, /rec'
        print,''
        return
    endif
    convertname, gname, z,ion
    zion2filename,z,ion,fname
    if keyword_set(ci) then begin
        filename = fname+'.cilvl'
    endif else if keyword_set(rec) then begin
        filename = fname+'.reclvl'
    endif else begin
        status = 0
        cireclvlStr = {status:status}
        print,' either ci or rec must be set'
    endelse
;   
    ntemp = 0
    s = ' '
    if file_test(filename) eq 1 then begin
;         print,' file found = ',filename
        status = 1
        openr,lur,filename,/get_lun
        nlines = 0
;         ndata = 0
        notfound = 1
        while not eof(lur) do begin
            readf,lur,s
;             print, s
            tst = strtrim(strmid(s,0,5))
            if strpos(tst,'-1') gt 0 and notfound then begin
                ndata = nlines
                notfound = 0
            endif else begin
                ntemp1 = n_elements(strsplit(strmid(s,14),' ',/extract))
                ntemp = max([ntemp,ntemp1])
            endelse
            nlines = nlines+1
        endwhile
        ;
        temp = fltarr(ntemp,ndata/2)
        rate = fltarr(ntemp,ndata/2)
        ntemp = intarr(ndata/2)
        lvl1 = intarr(ndata/2)
        lvl2 = intarr(ndata/2)
        
        point_lun,lur,0
        
        irate=0
        for i=0,ndata-1,2 do begin
            readf,lur,s
            info = fix(strsplit(strmid(s,0,13),' ',/extract))
            lvl1[irate] = info[2]
            lvl2[irate] = info[3]
            temp[0,irate] = float(strsplit(strmid(s,14),' ',/extract))
            grate = where(temp[*,irate] gt 0.,cnt)
            ntemp[irate] = cnt
            readf,lur,s
            rate[0,irate] = float(strsplit(strmid(s,14),' ',/extract))
            irate = irate + 1
        endfor
        ;
        free_lun,lur
        cireclvlStr = {filename:filename, temp:temp, rate:rate, ntemp:ntemp, lvl1:lvl1, lvl2:lvl2, status:status}
    endif else begin
        status = 0
        cireclvlStr = {status:status}
    endelse
end
    