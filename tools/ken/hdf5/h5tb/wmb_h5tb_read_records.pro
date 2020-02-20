
;
; wmb_h5tb_read_records
; 
; Purpose: Read records. 
;
; Description: wmb_h5tb_read_records reads some records identified from a 
;              dataset named dset_name attached to the object specified by 
;              the identifier loc_id. 
;              
; Parameters:
; 
; loc_id
;     IN: Identifier of the file or group in which the table is located. 
; dset_name
;     IN: The name of the table. 
; start
;     IN: The position to start reading. 
; nrecords
;     IN: The number of records to read. 
; databuffer
;     OUT: Output data buffer. 
;

pro wmb_h5tb_read_records, loc_id, $
                           dset_name, $
                           start, $
                           nrecords, $
                           databuffer

    compile_opt idl2, strictarrsubs

    ; get the number of records and fields

    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, nrecords_orig
    
    ; open the dataset
    
    did = h5d_open(loc_id, dset_name)
    
    ; read the records
    
    wmb_h5tb_common_read_records, did, $
                                  start, $
                                  nrecords, $
                                  nrecords_orig, $
                                  databuffer

    ; close
    
    h5d_close, did
    
end

