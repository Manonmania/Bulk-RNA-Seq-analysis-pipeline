# Download ribo-seq data from van Heesch data (TDF files)
 
spack load igvtools@2.3.98
 
igvtools tdftobedgraph hs_lv_ribo_pooled_fwd.tdf hs_lv_ribo_pooled_fwd.bedgraph
 
spack load ucsc-utils
 
sort -k1,1 -k2,2n hs_lv_ribo_pooled_fwd.bedgraph > hs_lv_ribo_pooled_fwd_sorted.bedGraph
 
bedGraphToBigWig hs_lv_ribo_pooled_fwd_sorted.bedGraph chrom.sizes hs_lv_ribo_pooled_fwd_sorted.bw

