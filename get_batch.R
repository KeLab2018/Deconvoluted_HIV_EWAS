
x=do.call(rbind, strsplit(phe$supplementary_file, "\\/"))[,9]; batch.info = do.call(rbind, strsplit(x, "_"))[,2]

