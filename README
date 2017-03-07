Python Package for ChaRM (Characterization of Reflectors and Modeling)

For more technical details, read
(https://www.slim.eos.ubc.ca/Publications/Public/Conferences/SEG/2008/maysami08SEGlcf/maysami08SEGlcf.pdf)
and
(https://www.slim.eos.ubc.ca/Publications/Public/Thesis/2008/maysami08THlcs.pdf)


AUTHOR 
	Mohammad Maysami
	Seismic Laboratory for Imaging and Modeling (SLIM)
	Department of Earth & Ocean Sciences (EOSC)
	The University of British Columbia (UBC)
	
LICENSE
	Copyright (C) 2007
	All Rights reserved.

	You may use this code only under the conditions and terms 
	of the license contained in the file LICENSE provided with
	this source code. If you do not agree to these terms you may 
	not use this software package.

	This software package is provided by the copyright holders and contributors 
	"AS IS"  and any express or implied warranties, including, but not limited to, 
	the implied warranties of MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
	PURPOSE are 	disclaimed. In no event shall the copyright holder or contributors 
	be liable for any direct, indirect, incidental, special, exemplary, or consequential 
	damages (including, but not limited to, procurement of substitute goods or 
	services; loss of use, data, or profits; or business interruption) however caused 
	and on any theory of liability, whether in contract, strict liability, or tort (including 
	negligence or otherwise) arising in any way out of the use of this software, even 
	if advised of the possibility of such damage.
	

PREREQUISITES
	-> Python-2.4 or newer (http://www.python.org/)
	-> NumPy-1.0 (http://numpy.scipy.org/)
	-> SciPy-0.5.1 (http://www.scipy.org/)
	-> Matplotlib-0.90.1 or newer (http://matplotlib.sourceforge.net/)
	-> MADAGASCAR-2817 (SVN developer tree at http://rsf.sourceforge.net/)
			Note: Compile with API=python,c++,f90,matlab (only python API is required)
	-> [OPTIONAL] SCons-0.96.95 or newer (http://www.scons.org/)
			Note: It is optional and only required to run SConstruct script.



INSTALLATION NOTES
	-> [Preferred] Just put this folder (Charm) into your $PYTHONPATH. 
		For example, if you are using BASH and $Charm is the path to your 
		Charm folder, then type:
			export PYTHONPATH=$PYTHONPATH:$Charm
		For CShell:
			setenv PYTHONPATH $PYTHONPATH:$Charm

	-> Alternatively, you may use setup.py to install the package in your 
		default python path location, enter this command:
			python setup.py install
	-> To install it into a specific folder type this, where /path/to/mypython 
	   is a folder which is in your $PYTHONPATH environment.
			python setup.py install --install-lib=/path/to/mypython
	
	
	-> Set environment variables(paths) if they are different than default values. 
		These variables are CHARM_Data,CHARM_Results,CHARM_pydata,CHARM_Demos. 
		Otherwise it will be set to defaults which are folders with 
		the same name (Data, Results, and etc.) in $Chamrpy as mentioned above.
		
		For example, if you are using BASH and $DataPATH is the path where you have 
		input rsf data, then type:
			export CHARM_Data=$PYTHONPATH:$DataPATH
		For CShell:
			setenv CHARM_Data $PYTHONPATH:$DataPATH
		
	-> Global variable "_df_input" in __init__.py points to default input rsf file 
		for some of the functions of this package. For your convenience, it is 
		better to be set to a sample rsf file in order to skip declaring it in 
		input argument
				
RUNNING THE PACKAGE IN PYTHON
	-> Go into your python interpreter.
	-> Type "import Charm" at your python command prompt.


PROCESSING WITH STDIN & STDOUT:
	Use following commands in terminal to analyze seismic data and show results
	-> ./sfchar.py <input.rsf >output.rsf args=...
	-> ./sfshow.py <results.rsf           args=...
		Note: Arguments are optional and will be set to default if not provided. 
		For more details about arguments check the scripts header for documentation
 