from setuptools import setup, find_packages

setup(
    name='pysoftk',
    version='0.1.0',
    author='A. Santana-Bonilla',
    author_email='k2031560@kcl.ac.uk',
    packages= find_packages(exclude=["test", "test.*", "test2", "test.*"]),
    #license='MIT',
    description='KCL_PSO: Kings College London Particle Swarm Optimization',
    #long_description=open('README.txt').read(),
    install_requires=[
        'ase',
    ],
)
