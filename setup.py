"""
ChaRM Python Package for characterizing reflectivity models


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""

from distutils.core import setup
from sys            import argv
from os             import path

# Find the absolute path to the package folder
script=argv[0]
pathname = path.dirname(script)
packpath= path.abspath(pathname)


packages = [
    "Charm",
    "Charm.Core",
    ]

package_dir={'Charm':packpath}
package_data={'Charm':[
        '*.txt',
        '*.sh',
        'LICENSE',
        'SConstruct',
        'Doc/*.txt','Doc/*.doc',
        'pydata/*.pyd',
        'Data/*.pyd',
        'Results/*.pyd',
        ]}


setup(name='Charm',
      description='Python Package for ChaRM',
      long_description='',
      version='0.1.0',
      author = "Mohammad Maysami",
      author_email='mmaysami@eos.ubc.ca',
      maintainer='Mohammad Maysami',
      maintainer_email='mmaysami@eos.ubc.ca',
      download_url='',
      url='https://wave.eos.ubc.ca ' ,
      packages = packages,
      package_dir=package_dir,
      package_data=package_data,
      )



