import numpy as np
import pandas as  pd
import scanpy as sc 
import time
import datetime
import multiprocessing
import umap
import nltk
nltk.download('punkt')
from nltk.tokenize import word_tokenize
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from tqdm import tqdm
import random
import pickle

class SCellBOW:

	def __init__(self, read_path, write_path):
		self.adata_source = None
		self.adata_target = None
		self.corpus_not_shuffled_trn = None
		self.corpus_shuffled_trn = None
		self.read_path = ""
		self.write_path = ""

	def set_source_data(self, adata_source):
		self.adata_source = adata_source

	def set_target_data(self, adata_target):
		self.adata_target = adata_target

	def preprocessing(self, path, min_genes, min_cells, target_sum, n_top_genes, max_value):
		adata = sc.read(path)

		print("\n____Create unique index____ ")
		adata.var_names_make_unique()
		print("____Filtering the data____")
		print("pre filtering:",adata.shape)
		sc.pp.filter_cells(adata, min_genes=200)
		sc.pp.filter_genes(adata, min_cells=20)
		print("post filtering:",adata.shape)


		print("____Log normalizing____")
		sc.pp.normalize_total(adata, target_sum=1e4)
		sc.pp.log1p(adata)

		print("____Selecting highly variable genes____")
		print("pre:",adata.shape)
		sc.pp.highly_variable_genes(adata, n_top_genes = 500)
		adata.raw = adata
		adata = adata[:, adata.var.highly_variable]

		print("____Scaling the data____")
		sc.pp.scale(adata, max_value=10)
		adata.to_df().head()
	
		return adata


	## HELPER FUNCTIONS 
	# Create corpus
	def wordfreq_docs_list(self, df):
		corpus = []
		row = []
		s = ''
		names = {e:name for e,name in enumerate(df.columns)}
		for i in tqdm(iterable = df.itertuples()):
			for e,j in enumerate(i[1:]):
				temp = names[e]
				row += [temp] * int(j)
			corpus.append(row)
			s = ''
			row = []
		print('corpus created with size: ',len(corpus))
		return corpus


	## HELPER FUNCTIONS 
	#Shuffle corpus
	def shuf(self, cnsl):
		corpus_shuffled = []
		random.seed(0)
		for l in tqdm(range(len(cnsl))):
			random.shuffle(cnsl[l])
			s = ' '.join(cnsl[l])
			corpus_shuffled.append(s)
		return corpus_shuffled


	## HELPER FUNCTIONS 
	#Train Doc2vec
	def doc2vec(self, corpus=[],method=1,umap_verbose=False,want_ari=False,
							epochs = 20,original_labels = [],model_name='model',
							vec_size = 300):
		# tagging docs
		tagged_data = [TaggedDocument(words=word_tokenize(_d), tags=[str(i)]) for i, _d in enumerate(corpus)]
		print('all docs tagged with len',len(tagged_data) )
			
		max_epochs = epochs 
			
		alpha = 0.025   #The initial learning rate
		model = Doc2Vec(vector_size=vec_size,
						alpha=alpha, 
						min_alpha=0.00025,
						min_count=1,
						window=2,
						workers=8,
						seed = 0,
						dm = method)
			
		# Building a vocabulary
		model.build_vocab(tagged_data,update=False)
		print('vocab built')
		model.train(tagged_data,
					total_examples=model.corpus_count,
					epochs=epochs)
			
			# Save Model
		model.save(self.write_path + model_name)
		print('\nmodel trained')
		return model


	## ACTUAL FUNCTION
	def SCellBOW_source(self,
										epochs,
										vec_size
										):

		#rescale the data
		srcdata = self.adata_source.to_df()
		srclabel = self.adata_source.obs['celltype'].to_frame()
		print(srcdata.shape,srclabel.shape)

	
		scaler = MinMaxScaler(feature_range=(1,10))
		print(scaler.fit(srcdata))
		trainData = scaler.transform(srcdata)
		trainData = pd.DataFrame(trainData,columns=srcdata.columns)
		print(trainData.shape)
		print(trainData.head())


		# Shuffle the corpus
		self.corpus_not_shuffled_trn = self.wordfreq_docs_list(trainData)
		self.corpus_shuffled_trn = self.shuf(self.corpus_not_shuffled_trn)
		start_time = time.time() ## to find embedding

		# Train the model
		model_train = self.doc2vec(self.corpus_shuffled_trn,method = 1,epochs = epochs, 
		                           model_name='train_model',vec_size = vec_size) # model is written to folder here
		train_time = time.time()
		temp = (train_time - start_time)/60
		print("Time : {} mins".format(temp))
	




	def SCellBOW_test(self,
					svd_solver,
					n_neighbors,
					n_pcs,
					resolution,
					):
	
		# Load the source model
		big_model = Doc2Vec.load(self.write_path + 'train_model')

		#rescale the data
		dstdata = self.adata_target.to_df()
		dstlabel = self.adata_target.obs['celltype'].to_frame()
		print(dstdata.shape,dstlabel.shape)
	
		scaler = MinMaxScaler(feature_range=(1,10))
		print(scaler.fit(dstdata))
		trainData = scaler.transform(dstdata)
		trainData = pd.DataFrame(trainData,columns=dstdata.columns)
		print(trainData.shape)
		print(trainData.head())

		corpus_not_shuffled = self.wordfreq_docs_list(trainData)
		corpus_shuffled = self.shuf(corpus_not_shuffled)
	

		def save(data,path):
		  dbfile = open(path, 'wb') 
		  pickle.dump(data, dbfile, protocol=4)                      
		  dbfile.close()
		  return None
		
		# tokenizing the corpus
		tagged_data_tl = [TaggedDocument(words=word_tokenize(_d), tags=[str(i)]) for i, _d in enumerate(corpus_shuffled)]
		print('all docs tagged with len',len(tagged_data_tl) )
		save(tagged_data_tl, self.write_path + 'tagged_data')

		# Updating the vocabulary
		big_model.build_vocab(tagged_data_tl, progress_per=1000,update = True)
		print('vocab updated')

		# Retraining the model
		start_time = time.time()
		big_model.train(tagged_data_tl,
										total_examples= big_model.corpus_count,
										epochs=20)
		train_time = time.time()
		temp = train_time - start_time
		print("Time : {} secs".format(temp))

		# infer vectors for new corpus
		re_vecs_W_new= []
		big_model.random.seed(0)
		for c in tqdm(corpus_shuffled):
			c = c.split(' ')
			re_vecs_W_new.append(big_model.infer_vector(c))
		re_vecs_W_new = np.array(re_vecs_W_new)

		# Create anndata object for the embedding

		## vecdata = adata_test.copy()
		## vecdata.embedding = sc.AnnData(re_vecs_W_new)
		## sc.tl.pca(vecdata.embedding, svd_solver=svd_solver)

		# vecdata = self.adata_target.copy()
		# self.adata_target.embedding = re_vecs_W_new
		self.adata_target.obsm["embedding"] = re_vecs_W_new
		# vecdata.embedding = sc.AnnData(re_vecs_W_new)

		# vecdata_W_new  = sc.AnnData(re_vecs_W_new)

		sc.tl.pca(re_vecs_W_new, svd_solver=svd_solver)
		sc.pp.neighbors(self.adata_target, use_rep="embedding", n_neighbors=n_neighbors, 
		                n_pcs=n_pcs, random_state=1)
		sc.tl.umap(self.adata_target)

		sc.tl.leiden(re_vecs_W_new, key_added='clusters', resolution=resolution)
		with plt.rc_context({'figure.figsize': (8, 8)}):
			sc.pl.umap(self.adata_target, 
			      use_raw=False,
			      layer="embedding",
						color='clusters', 
						add_outline=True, 
						legend_fontsize=14, 
						legend_fontoutline=2,
						frameon=False,
						title='UMAP visualisation', 
						size = 50,
						palette=plt.rcParams["axes.prop_cycle"],
						)

