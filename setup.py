
from setuptools import setup, find_packages

setup(
    name='artemis', author="Jasper Linthorst", author_email="ja.linthorst@gmail.com",
    version='0.1',
    install_requires=[
        'pysam',
        'pybedtools',
        'argparse',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'artemis=artemis.__main__:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)

