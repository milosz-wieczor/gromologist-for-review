from setuptools import setup, find_packages

setup(name='gromologist',
      version='0.320',
      description='Library to handle various GROMACS-related stuff',
      author='Milosz Wieczor',
      author_email='milafternoon@gmail.com',
      license='GNU GPLv3',
      packages=find_packages(where="src"),
      package_dir={"": "src"},
      install_requires=['numpy>=1.10.0'],
      zip_safe=False)
