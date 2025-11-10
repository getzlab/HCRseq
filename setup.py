from setuptools import setup, find_packages

setup(
    name='hcrseq',
    version='',
    url='',
    license='',
    author='',
    author_email='',
    description='',
    packages=find_packages(),
    install_requires=[
        'Click',
    ],
    entry_points = {
        'console_scripts': [
            'hcrseq = hcrseq.cli:hcrseq',
        ]
    }
)
