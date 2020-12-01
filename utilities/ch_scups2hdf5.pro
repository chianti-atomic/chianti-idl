;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas.  See www.chiantidatabase.org 
;
;     		          
; Purpose     : converts a CHIANTI ascii SCUPS file to and HDF5 file
;
; Explanation : the  CHIANTI ascii SCUPS file is read into an IDL
;               structure, which is then written the the HDF5 file.
;    
; Inputs      : the file names
;
; Outputs     : none
;
; Calls:
;             read_scups and the standard IDL routines to write the
;             HDF5 file.
;
; Note:       to read the structure:
;             s = H5_PARSE('file.h5', /READ_DATA)

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


pro ch_scups2hdf5, file_scups, file_hdf

  ; this the structure with the data:
  read_scups,  file_scups,  st, /verb
  
  H5_OPEN
  
fileID = H5F_CREATE(file_hdf)  

datatypeID = H5T_IDL_CREATE(st)

dataspaceID = H5S_CREATE_SIMPLE(1)

datasetID = H5D_CREATE(fileID, 'SCUPS', datatypeID, dataspaceID)

H5D_WRITE, datasetID, st

H5F_CLOSE, fileID

; the following is necessary to clear the memory:
   H5_CLOSE


end
