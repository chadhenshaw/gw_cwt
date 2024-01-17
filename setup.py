from setuptools import setup, find_packages

setup(
    name='gw_cwt',
    version='0.1.0',  
    author='Chad Henshaw',
    author_email='cgh3@gatech.edu',
    description='Continuous wavelet transform for gravitational wave data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/chadhenshaw/gw_cwt',  # Replace with the real URL
    packages=find_packages(),  # This automatically finds the packages in your project
    install_requires=[
        'numpy',
        'scipy',
        'pycbc',
        'h5py',
        'matplotlib'       
    ],
    python_requires='>=3.6',  # Specify the minimum Python version required
)
