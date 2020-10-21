;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas.  See www.chiantidatabase.org 
;
;     		          
; Purpose     : converts the CHIANTI ascii SCUPS files to and HDF5 files
;
; Explanation : each CHIANTI ascii SCUPS file is read into an IDL
;               structure, which is then written the the HDF5 file.
;               The routine goes through the masterlist of the ions.
;
; Inputs      : none
;
; Optional Input: the top-level CHIANTI database directory
;
; Outputs     : none
;
; Calls:
;             read_masterlist, ch_scups2hdf5
;
; 
; Written     : 
;       Version 1, Giulio Del Zanna (GDZ)  13 Oct 2020
;
; Modified    :
;
; VERSION     :    V.1, 13 Oct 2020
;
;-        
;---------------------------

pro ch_all_scups2hdf5, dir=dir

  if n_elements(dir) eq 1 then begin

     dir1=!xuvtop 

     !xuvtop=dir
  end
  

  read_masterlist,'',list

  nlist=n_elements(list)

  FOR ilist=0L,nlist-1 do BEGIN

   err = 0

   gname=list[ilist]
   
   ion2filename,gname,fname
   
   file_scups =fname+'.scups'
   file_hdf5=fname+'.scups.h5'

   print,'converting '+file_scups

   ch_scups2hdf5, file_scups, file_hdf5


ENDFOR  ; main loop for the ions

 if n_elements(dir) eq 1 then !xuvtop=dir1
  
end
