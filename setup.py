import os
from setuptools import setup

# Get the directory containing the module
module_dir = os.path.dirname(__file__) or os.path.dirname(__path__[0])

setup(
    name='protein2genome',
    version='0.1.0', 
    packages=['protein2genome', 'protein2genome.resources'],
    include_package_data=True,
    install_requires=[
        'argparse',
        'pandas',
        'pathlib',
        'importlib',
        'tqdm'
    ],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'protein2genome=protein2genome.__main__:main'
        ]
    },
    author='Jiawei Gu',
    author_email='jiawku17@gmail.com',
    description='A library for converting protein coordinates to genomic coordinates.',
    long_description=open(os.path.join(module_dir, 'README.md')).read(),
    long_description_content_type='text/markdown',
    url='https://github.com/jiawku/protein2genome',  # Update with your repository URL
    license='MIT', 
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='Convert protein coordinates to genomic coordinates.',
)
