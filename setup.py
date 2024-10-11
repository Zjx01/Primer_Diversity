from setuptools import setup, find_packages

setup(
    name='Primer_Diversity',
    version='1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'Primer_Diversity=Primer_Diversity.calculate_diversity:main',
        ],
    },
    install_requires=[
	'numpy',
        'pandas',
        'biopython',
        'seaborn',
	'logomaker',
	'matplotlib',
    ],
    python_requires='>=3.12.2',
)

