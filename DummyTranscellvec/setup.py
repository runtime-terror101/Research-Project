from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'SCellBOW'
LONG_DESCRIPTION = 'The first Python package with a longer description'

# Setting up
setup(
        name="SCellBOW", 
        version=VERSION,
        author="Juhi Pandey",
        author_email="<youremail@email.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['gensim', 'leidenalg', 'matplotlib', 'numpy', 'pandas', 'scanpy', 'umap', 'nltk', 'datetime', 'sklearn', 'tqdm'],
                          # add any additional packages that 
                          # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'SCellBOW'],
        classifiers= [ # ???
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)