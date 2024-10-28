from setuptools import setup

setup(name='gromologist',
      version='0.310',
      description='Library to handle various GROMACS-related stuff',
      author='Milosz Wieczor',
      author_email='milafternoon@gmail.com',
      license='GNU GPLv3',
      packages=['gromologist'],
      install_requires=['numpy>=1.10.0'],
      zip_safe=False)
