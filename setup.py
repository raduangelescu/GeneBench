import setuptools

setuptools.setup(
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        'GEOparse==2.0.2',
        'matplotlib==3.3.3',
        'numpy==1.17.4',
        'pandas==1.2.1',
        'pymongo==3.11.2',
        'rpy2==3.3.5',
        'scipy==1.6.0',
        'sklearn==0.0',
        'tensorflow_gpu==2.4.1',
        'xgboost==1.3.3'
    ]
    )
