import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
	name="VariaNNNts",
	version="1.0.0",
	author="Alexander G. McFarland",
	author_email="alex.925@gmail.com",
	description="Create synthetic reads with specific variants",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/agmcfarland/VariaNNNts",
	project_urls={
		"Bug Tracker": "https://github.com/agmcfarland/VariaNNNts/issues",
	},
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
	],
	package_dir={"": "src"},
	packages=setuptools.find_packages(where="src", exclude=['docs','test']),
	entry_points={
		# 'console_scripts': ['GeneGrouper = GeneGrouper.__main__:main'], },
		'console_scripts': ['VariaNNNts = VariaNNNts.__main__:main'], },
	package_data={'VariaNNNts' : ['libs/*']},
	python_requires=">=3.6")

# https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments
# https://packaging.python.org/tutorials/packaging-projects/
# https://trstringer.com/easy-and-nice-python-cli/
# https://github.com/pypa/sampleproject/blob/main/setup.py
# https://chriswarrick.com/blog/2014/09/15/python-apps-the-right-way-entry_points-and-scripts/