from SCellBOW import SCellBOW

scb = SCellBOW("C:/Users/juhip/Desktop/", "C:/Users/juhip/Desktop/")

adata_source = scb.preprocessing("C:/Users/juhip/Desktop/Bh.h5ad", 200, 20, 1e4, 500, 10)
adata_target = scb.preprocessing("C:/Users/juhip/Desktop/smartseq2.h5ad", 200, 3, 1e4, 500, 10)
scb.set_source_data(adata_source)
scb.set_target_data(adata_target)

scb.SCellBOW_source(epochs=40, vec_size=300)

# scb.SCellBOW_test('arpack', 15, 40, 1.0)

