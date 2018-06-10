import os

from setuptools import setup


version = __import__('crypko').__version__


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()


setup(
    name='crypko.py',
    version=version,
    author='Bottersnike',
    author_email='bottersnik237@gmail.com',
    description=('A basic wrapper around the Crypko platform.'),
    long_description=read('README.rst'),
    license='MIT',
    install_requires=[
        'requests',
        'web3'
    ],
    packages=['crypko', 'crypko.examples'],
    url='https://gitlab.com/Bottersnike/crypko.py',
    project_urls={
        "Bug Tracker": "https://gitlab.com/Bottersnike/crypko.py/issues",
        "Source Code": "https://gitlab.com/Bottersnike/crypko.py",
    },
    package_data={
        'crypko': ['abi.json'],
    },
)

