from distutils.core import setup, Extension

module = Extension("Example", sources = ["Example.c"])

setup(name="PackageName",
      version='0.01',
      description='This is the example module we created for C',
      ext_modules=[module])