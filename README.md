#  <span style="color:blueviolet">Genotypes.py</span>
Define alleles and infer genotypes for deeply sequenced genetic loci (PCR amplicons).

 
## <span style="color:mediumorchid">Table of contents</span>  
* [Background](#background)
* [Features](#features)  
* [Requirements](#requirements)
* [Synopsis](#synopsis)
* [System setup](#system%20setup)
* [Code launch notes](#codelaunchnotes)
* [Operation notes](#operationnotes)
* [Input notes](#input%20notes)
* [Output notes](#output%20notes)
* [Visual summary of key script operations](#visual%20summary%20of%20key%20script%20operations)
* [Status](#status)
* [Contact](#contact)

## <span style="color:mediumorchid">Background</span> 
This script returns allele definitions (and extrapolated genotypes) for amplicons sequenced by Illumina® short-read ('next-generation') sequencing technologies. Users input demultiplexed fastq files for sequenced amplicons and specify a BLASTN database to be used as an alignment reference, and the script completes allele definitions (wild-type, deletion, insertion, substitution, *etc.*) for the sequenced locus based on relative read abundance in each fastq file. Based on allele definitions, the script also infers genotype (homozygous, heterozygous) across the sequenced locus for each sample.
<br clear="all" />
<img src="Genotypes_img/Genotypes_thumbnail_sketch.png" align="left" width="750">
<br clear="all" />
## <span style="color:mediumorchid">Features</span> 
* Automates allele definitions and infers genotypes at specific loci, for amplicons deeply sequenced on Illumina® platforms
* Key input: demultiplexed fastq files
* Key outputs: sequence alignments for candidate alleles (reads ranked by abundance), with optional DNA sub-sequence(s) mapped on alignments (*e.g.*, Cas9 guide RNA sequence(s), DNA sub-sequence(s) to test for presence/ablation); inferred genotypes

## <span style="color:mediumorchid">Requirements</span>
* Python 3.7 or higher
* BLASTN (NCBI) (available for OS-appropriate download as part of BLAST+ suite @ <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">Download BLAST Software and Databases</a>)

## <span style="color:mediumorchid">Synopsis</span>
**This script returns inferred genotypes for sample-specific amplicons deeply sequenced on Illumina® sequencing platforms.**
>(see 'Output notes' for file output
details).  


**Users are asked for paths to (1) an output directory, (2) fastq directory (sample sequence sources), (3) BLASTN executable, and (4) BLAST reference sequence database (+prefix common to database files), as well as (optional) DNA sub-sequence(s) to query in read alignments (*e.g.*, Cas9 guide RNA sequence(s), DNA test sequence(s)).**

>(see 'Input notes' for details).

For further usage details, please refer to the following manuscript:  
>*Ehmsen, Knuesel, Martinez, Asahina, Aridomi, Yamamoto (2020)*
    
Please cite usage as:  
>Genotypes.py  
>*Ehmsen, Knuesel, Martinez, Asahina, Aridomi, Yamamoto (2020)*


## <span style="color:mediumorchid">System setup</span>

#### <span style="color:slateblue">1. Confirm that Python 3 and Jupyter Notebook are available on your system, or download & install</span>
##### Mac or Windows OS

**Python 3**  
Mac OS generally comes with Python pre-installed, but Windows OS does not.  Check on your system for the availability of Python version 3.7 or higher by following guidelines below:
 
- First open a console in Terminal (Mac OS) or Command Prompt/CMD (Windows OS), to access the command line interface.
 
-  Check to see which version of Python your OS counts as default by issuing the following command (here, `$` refers to your command-line prompt and is not a character to be typed):  
     
	`$ python --version` 
	 
 	- If the output reads `Python 3.7.3` or any version >=3.7, you are good to go and can proceed to **Jupyter Notebook**.
 
 	- If the output reads `Python 2.7.10` or anything below Python 3, this signifies that a Python version <3 is the default version, and you will need to check whether a Python version >=3.7 is available on your system.
 		- To check whether a Python version >=3.7 is available on your system, issue the following command:  
 
 			`$ python3 --version`  
 
 		- If the output finds a Python version >=3.7 (such as `Python 3.7.3`), you are good to go and can proceed to **Jupyter Notebook**.
		- If the output does *not* find a Python version >3.7, use one of the following two options to download and install Python version >=3.7 on your computer: 
  
 			1- *To install Python 3 alone (separately from Jupyter Notebook)*
 
 			[Python](https://www.python.org/downloads/) https://www.python.org/downloads/
 			- "Download the latest version for Mac OS X", and then follow installation guidelines and prompts when you double-click the downloaded package to complete installation.  
 
 			- Once you have downloaded and installed a Python 3 version >=3.7, double-check in your command-line that Python 3 can be found on your system by issuing the following command:  
 	
 				`$ python3 --version` (Mac OS)  
 				`$ python --version` (Windows OS)
  
 				The output should signify the Python version you just installed.  Proceed to **Jupyter Notebook**.
 
 			2- *To install Python 3 and Jupyter Notebook all at once* (**recommended**)
 
			[Anaconda (with Jupyter Notebook) Download & Installation](https://jupyter.readthedocs.io/en/latest/install/notebook-classic.html) https://jupyter.readthedocs.io/en/latest/install/notebook-classic.html

 			- Download Anaconda with Python 3, and then follow installation guidelines and prompts when you double-click the downloaded package to complete installation.  
 
 			- Once you have downloaded and installed Python 3, double-check in your command-line that Python 3 can be found on your system by issuing the following command: 
  	
				`$ python3 --version` (Mac OS)  
				`$ python --version` (Windows OS)  

 			- Also double-check in your command-line that Jupyter Notebook can be found on your system by issuing the following command:  
 	
 				`$ which jupyter notebook`
 
			- If the output indicates that 'jupyter' is available in the path of your Python 3 installation (such as, `/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter`), you are good to go and can proceed to **2. Configure Python Virtual Environment for Jupyter Notebook**.  

**Jupyter Notebook**  
*Note, these steps are not required if you installed Anaconda with Jupyter Notebook as recommended above*.  Jupyter Notebook is *not* generally pre-installed on Mac OS.  To check whether Jupyter Notebook is available with Python 3 on your machine, issue the following command:
 
 	$ which jupyter notebook
 
If the output indicates that 'jupyter' is available in the path of your Python 3 installation (such as, `/Library/Frameworks/Python.framework/Versions/3.7/bin/jupyter`), you are good to go and can proceed to **2. Configure Python Virtual Environment for Jupyter Notebook**.  If instead you see an error message indicating that 'jupyter notebook' is not available, issue the following commands using the Python program installer (pip) to install Jupyter: 
  
 	$ pip3 install --upgrade pip  
	$ pip3 install jupyter


<br> 
#### <span style="color:slateblue">2. Download the Genotypes repository (or program file subset) from GitHub</span>

Genotypes.py can be accessed as a **Jupyter Notebook** or **Python program file** available at [YamamotoLabUCSF GitHub](https://github.com/YamamotoLabUCSF/Genotypes) (https://github.com/YamamotoLabUCSF/Genotypes).  Please note that access to the file through GitHub requires a personal GitHub account.  

* (a) **Create a personal GitHub account** (free)  
	* follow instructions available at [WikiHow: How to Create an Account on GitHub](https://www.wikihow.com/Create-an-Account-on-GitHub) https://www.wikihow.com/Create-an-Account-on-GitHub. 

* (b) **Navigate to the code repository for Genotypes**
	* The code repository contains:
		* the **Jupyter Notebook** file (.ipynb)
		* the **Python program** file (.py)
		* **image files** associated with the Jupyter Notebook
		* a **requirements** file (Genotypes_requirements.txt), used in the creation of a Python virtual environment to run Genotypes (see **3. Configure Python Virtual Environment for Jupyter Notebook**, below)
		
* (c) **Download** or **clone** the repository to your personal computer:  

	**Download:**   
	* first click on the repository name	 to access the repository
	* then download the entire repository directory and its associated subdirectories and files (**green download icon, labeled "Code"**)
	* alternatively, download only the target files you need for your intended purposes 
	  * for example, download the **Jupyter Notebook** file (.ipynb), **image files directory**, and **requirements** file if you plan to use only the Jupyter Notebook

	**Clone:**   
	  
	* first click on the repository name	 to access the repository
	* then click on the right-hand arrow of the **green download icon, labeled "Code"**, to access the drop-down menu with **Clone** options (HTTPS or GitHub CLI).
	* If selecting the HTTPS option, copy the indicated URL and paste it at your command-line as an argument to the command 'git clone':  
	`$ git clone https://github.com/YamamotoLabUCSF/Genotypes.git` 

* (d) **Choose a directory location** on your machine where you would like to store the downloaded or cloned repository and its files.  This can be any folder/directory location you like.  **Move the repository files from the directory into which they were downloaded or cloned, into this directory**, if you have not already downloaded or cloned the repository/files directly into your target directory.  

	*Example using command line code*  
	*(directory can be created and accessed using command line prompts)*   
	* For example, on Mac OS:  
	
		* To create an empty directory named 'GenotypesCode', in the 'Documents' directory:*  
	`$ mkdir /Users/yourusername/Documents/GenotypesCode`
		* To navigate to the directory named 'Genotypes':	`$ cd /Users/yourusername/Documents/GenotypesCode`

	* For example, on Windows OS:  
	 
		* To create an empty directory named 'GenotypesCode', in the 'Documents' directory:  
	`$ mkdir C:\Documents\Genotypes`
		* To navigate to the directory named 'GenotypesCode':  
		`$ cd C:\Documents\GenotypesCode`

* (e) With the Genotypes repository files now present in this directory you've created and named, now navigate to the **Genotypes** repository in the command line, by issuing the following command:  

	`$ cd /Users/yourusername/Documents/GenotypesCode/Genotypes` (Mac OS)  
	or  
	`$ cd C:\Documents\GenotypesCode\Genotypes` (Windows OS)
 
<br>
 
#### <span style="color:slateblue">3. Create a Python Virtual Environment for Jupyter Notebook</span>

You are now ready to install a couple of **additional Python modules** that Genotypes.py requires for operation.  These Python modules can be installed from the Python Package Index repository ([PyPI](https://pypi.org/)) (https://pypi.org/) by individual download from the PyPI website and installation on your machine.  However, to simplify and streamline installation of these modules, we will create a Python **virtual environment** (self-contained 'directory' with all the Python modules needed to run Genotypes.py), using the following approach:

1.  First, install the Python module **virtualenv** ([virtualenv](https://pypi.org/project/virtualenv/)) (https://pypi.org/project/virtualenv/), by issuing the following command at the command line: 
 
	`$ pip3 install virtualenv`  (Mac OS)  
	or  
	`$ pip install virtualenv`  (Windows OS)
	
	pip3 (Mac OS) or pip (Windows OS) is Python 3's installation manager, and as long as there is an internet connection available, pip3 (Mac OS) or pip (Windows OS) will access the specified module from PyPI (here, virtualenv) and install it for access by Python 3.

2. Next, choose a **directory location** on your machine where you would like to install the files associated with a virtual environment.  This can be any folder/directory location you like (for example, you may have a favorite directory where other Python virtual environments are stored).  Alternatively, simply create the Python virtual environment in the Genotypes directory you created above (in section 2d).  At the command line, navigate to the location of this directory.  
	* For example, on Mac OS:  
		* To navigate to the directory named 'GenotypesCode':	`$ cd /Users/yourusername/Documents/GenotypesCode`

	* For example, on Windows OS:  

		* To navigate to the directory named 'GenotypesCode':  
		`$ cd C:\Documents\GenotypesCode`

3. With this directory set as your working location in the command line, now issue the following commands to **create a virtual environment**:  
	
	(a) Create a Python virtual environment named **Genotypes_env**, specifying that the environment will use Python 3:  
	`virtualenv -p python3 Genotypes_env` 
	
	(b) Activate Genotypes_env:   
	`source Genotypes_env/bin/activate`  (Mac OS)  
	or  
	`.\Genotypes_env\Scripts\activate`  (Windows OS) 	
	
	(c) You should now see that your command line prompt has changed to indicate that you are in the virtual environment, displaying something along the lines of:  
	`(Genotypes_env) $`
	
	(d) Now install the Python modules required by Genotypes.py, using the requirements file named **Genotypes_requirements.txt**, located in the Genotypes repository:
	   
	`$ pip3 install -r Genotypes/Genotypes_requirements.txt` (Mac OS)   
	or  
	`$ pip install -r Genotypes/Genotypes_requirements.txt`  (Windows OS)  	
	
	(e) To check that the required Genotypes.py Python modules were installed properly, now issue the follow command:
	
	`$ pip3 list`  (Mac OS)  
	or  
	`$ pip list`  (Windows OS)  	
	
	If the virtual environment was set up successfully, the output will read: 
	  
  **Package    Version**    
  certifi         2020.6.20  
  cycler          0.10.0  
  fpdf            1.7.2  
  kiwisolver      1.2.0  
  matplotlib      3.3.2  
  numpy           1.19.2  
  pandas          1.1.2  
  Pillow          7.2.0  
  pip             20.3.1  
  psutil          5.7.2  
  pyparsing       2.4.7  
  PyPDF2          1.26.0  
  python-dateutil 2.8.1  
  pytz            2020.1  
  scipy           1.5.2  
  setuptools      50.3.2  
  six             1.15.0  
  wheel           0.36.0   


4.  Finally, Jupyter Notebook needs to be made aware of the Python virtual environment you just created.  To accomplish this, issue the following commands:  

    (Mac OS)   
	```$ pip3 install ipykernel```   
	```$ python -m ipykernel install --name=Genotypes_env```  

    (Windows OS)   
	```$ pip install ipykernel```   
	```$ python -m ipykernel install --user --name Genotypes_env``` 


5.  You should now be ready to **access Genotypes.py** in Jupyter Notebook or at the command line!  See **Using Genotypes.py**.

6. Just for completeness, to **exit the Python virtual environment and return to your 'native' (default) environment**, simply issue the following command at the command line:

	``$ deactivate`` 
	
	To **re-enter the virtual environment** at any time in the future, you would use the command in 3b:
	
	`source Genotypes_env/bin/activate`  (Mac OS)  
	or  
	`.\Genotypes_env\Scripts\activate`  (Windows OS) 

7. Also, note that if you'd like **to remove the Genotypes_env at any time**, you can do so by issuing the following command:  

	```$ rm -rf Genotypes_env``` 

	This would delete the virtual environment from your machine.

#### <span style="color:slateblue">4. Download & install external dependencies (BLAST, MEME)</span>

*Additional setup*:  

* Locally install **BLASTN** (see Requirements)  
<img src="Genotypes_img/BLASTN_thumbnail.png" align="left" width="100">
<br clear="all" />

* Download or create a **reference sequence database** (required for BLASTN alignment operations)  
*this can be obtained in one of two ways:*  
  * create custom database from a fasta file containing reference sequence(s) using MAKEBLASTDB (NCBI)  
(details at <https://www.ncbi.nlm.nih.gov/books/NBK279688/>)
  * download pre-formatted NCBI BLAST database (details at <https://www.ncbi.nlm.nih.gov/books/NBK537770/>)  
   
<img src="Genotypes_img/BLASTN_reference_database_thumbnail.png" align="left" width="150">
<br clear="all" />
    * *Note, separate installation of these external dependencies is not necessary if using the virtual machine file Alleles\_and\_altered\_motifs.ovf*



##<span style="color:mediumorchid">Code launch notes</span>  
Code is available as a Jupyter Notebook file (**Genotypes.ipynb**) or as a Python program file (**Genotypes.py**) for direct use, or pre-packaged with all dependencies as an Open Virtualization Format file for virtual machines (**Alleles\_and\_altered\_motifs.ovf**).  

#### <span style="color:slateblue">Jupyter Notebook or Python program file</span>

In **System Set-Up** above, you downloaded and installed Python 3, Jupyter Notebook, & the Genotypes code repository, and created a Python virtual environment (Genotypes\_env) containing the Python modules that Genotypes.py needs in order to run.  To access Genotypes.py (Jupyter Notebook or Python program file) for interactive work, proceed through guidelines indicated below.  

#### <span style="color:slateblue">Jupyter Notebook</span>

*Note*, Jupyter Notebook file requires *Genotypes_img* directory containing five image files to be available in the directory from which the Jupyter Notebook will be opened.  

**Anaconda Navigator**

1.  If you downloaded **Anaconda (with Jupyter Notebook)** as recommended in System Setup above, launch the Anaconda-Navigator application.  Click on **Jupyter Notebook** to activate Jupyter Notebook.  

2.  You will see a new **tab open automatically in your default web browser** (such as Chrome), in which there is a **directory tree** illustrating the current file contents of the **current working** directory.  Navigate to the directory containing **Genotypes.ipynb**.  Click on the **Genotypes.ipynb** file to open it as a Jupyter Notebook.

3. In the Jupyter Notebook menu bar at the top of the newly opened Genotypes.py Jupyter Notebook file, take the following steps to **specify Genotypes\_env as your Python virtual environment** (see System Setup: 3. Create a Python Virtual Environment for Jupyter Notebook):  

	(a) click on the **'Kernel'** tab to open a drop-down menu  
	(b) hover your mouse over **'Change Kernel'** (bottom option of drop-down menu) to identify the kernels available to the Notebook  
	(c) choose **Genotypes\_env** as the kernel  
	 
	 If Genotypes\_env does not appear as a kernel option, troubleshooting is needed to either create the Genotypes\_env virtual environment, or to make Jupyter Notebook aware of the existence of the Genotypes_env virtual environment.

4.  If these steps have been accomplished successfully, you are now ready to use the Genotypes.py Jupyter Notebook.  Be prepared to provide required input variables as detailed below in **Input Notes**.
<br>

**Command Line**  

1.  If you plan to activate Jupyter Notebook from your command line, first navigate to the directory containing **Genotypes.ipynb**.  For example, if you downloaded the code repository as **Genotypes** and are not already in **Genotypes** as your working (current) directory at the command line, **navigate** there at this time using commands similar to the following:

	* For example, on Mac OS:  
	`$ cd /Users/yourusername/Documents/GenotypesCode/Genotypes`

	* For example, on Windows OS:  
	`$ cd C:\Documents\GenotypesCode/Genotypes`
	
	You can **check your current working directory** by issuing the following command:  
	
	`$ pwd`

2.  To further confirm that you are in the repository directory at the command line, you can **check the files** present in the repository directory by issuing the following command:  

	`$ ls` (Mac OS)  
	`$ dir` (Windows OS)	
	
	You should see a list of files that includes the Genotypes Jupyter Notebook and/or Python program file (depending on whether you downloaded the repository or individual program files), named **Genotypes.ipynb** or  **Genotypes.py**
	
3.  Activate Jupyter Notebook by issuing the following command:  

	`jupyter notebook`  (Mac OS, in Terminal)
	
4.  Refer to steps 2-4 under "Anaconda Navigator" approach above.  

#### <span style="color:slateblue">Python program file</span>
**Anaconda Navigator**

1.  If you downloaded **Anaconda (with Jupyter Notebook)** as recommended in System Setup above, launch the Anaconda-Navigator application.  Click on **Jupyter Notebook** to activate Jupyter Notebook.  

2.  You will see a new **tab open automatically in your default web browser** (such as Chrome), in which there is a **directory tree** illustrating the current file contents of the **current workin

**Command Line**  

1.  If you plan to run Genotypes.py from your command line, first navigate to the directory containing **Genotypes.py**.  Prepare access to a Python virtual environment containing appropriate packages required by Genotypes.py, as described in System Setup.  

    (a) Activate Genotypes\_env:   
	`source Genotypes_env/bin/activate`  (Mac OS)  
	or  
	`.\Genotypes_env\Scripts\activate`  (Windows OS) 	
	
	(b) You should now see that your command line prompt has changed to indicate that you are in the virtual environment, displaying something along the lines of:  
	`(Genotypes_env) $`

2.  To run Genotypes.py, type **python3 Genotypes.py** and hit 'Enter':
	`(Genotypes_env) $ python3 Genotypes.py`  
	
3.  If these steps have been accomplished successfully, you will now encounter the first interactive prompts of Genotypes.py.  Be prepared to provide required input variables as detailed below in **Input Notes**.



#### <span style="color:slateblue">Virtual machine (Alleles\_and\_altered\_motifs.ovf)</span>

1. To access the virtual machine that comes pre-installed with external dependencies, **download the file Alleles\_and\_altered\_motifs.ovf** (**Zenodo, DOI 10.5281/zenodo.3406862**). *Note: this requires virtualization software, such as Oracle VM VirtualBox*

2. Instructions for how to set up the virtual machine (including username and password for ovf file) can be found in the Zenodo repository.    



## <span style="color:mediumorchid">Operation notes</span>
*What does this script do?*  

 1. **classify & count reads:** counts unique read types per well (*i.e.*, sample ID); fastq file name provides the sample name  
 
 
 2. **identify top 10 reads** per sample ID (in terms of read abundance); calculates representation among reads within the well at four levels:
 
   * raw frequency (% read type in question, relative to total reads)  
   * percentile (% of other read types that fall below the frequency of the read type in question)  
   * adjusted frequency @ 1% (% read type in question, relative to reads that occur at >1% total frequency)  
   * adjusted frequency @ 10% (% read type in question, relative to reads that occur at >10% total frequency)  
   
   
 3. **align to reference database:** aligns top 10 reads to reference database using BLASTN  
 *(National Center for Biotechnology Information;
    Altschul S.F. et al. (1990) "Basic local alignment search tool")*    
    
    
 4. **return alignments as alleles & extrapolated genotypes;**  
 **(optional) map sub-sequence(s) onto alleles:**  
    -  for mutants, map location of Cas9 cut(s) and indel(s) relative to wt,
       if Cas9 guide sequence(s) supplied by user  
    -  map location of test sub-sequence(s) and whether sub-sequence is altered (ablated),
       if test sub-sequence(s) supplied by user  
       
       
 5. **provide overall population statistics:**  
 
   (a) total sample # for which genotypes were inferred  
   (b) distribution of genotypes among samples (homozygous, heterozygous, *etc.*)  
   (c) estimated wild-type *vs.* mutant allele frequencies  
   (d) summary of samples and reads that either had 'no hit' in reference database provided to BLASTN,
       or multiple hits (>1)  
<br clear="all" />        
**Operations overview:** *See 
Input notes' and 'Output notes'*  
Files labeled *i-iv* below: Key output files containing script interpretations of sample alleles & genotypes.
<img src="Genotypes_img/Genotypes_thumbnail_2.png" align="left" width="650">
<br clear="all" />
 

## <span style="color:mediumorchid">Input notes</span>
You will be prompted for the following user-specific information (up to 7 items):

   **Required** (4 strings specifying directory or executable locations, 1 string specifying sequence database file prefix):  
      <ul>
      <li>where should output files go?</li>
          *path to* **output directory** *for output files*
      <li>where are input files found?</li>
          *path to single directory containing* **demultiplexed fastq files**                                         
      <li>where is BLASTN executable found?</li>
          *path to* **BLASTN** *installation*
      <li>where is the reference sequence database used for alignment?</li>
          *path to directory containing six files that compose the* **reference sequence database** *used
    for BLASTN alignments (.nhr, .nin, .nog, .nsd, .nsi, .nsg)*
      <li>what prefix is common to the six files that compose the reference sequence database?</li>
          *prefix common to database files .nhr, .nin, .nog, .nsd, .nsi, .nsg*
      </ul>
                                                                                                           
   **Optional** (up to 2 lines of comma-separated strings specifying DNA sub-sequence(s):    
 **DNA sub-sequence(s)** to be mapped onto sequence alignments
      <ul>
      <li>**guide RNA sequence** (in 5'-3' DNA representation, excluding PAM sequence)</li>
      <li>**test sequence** (5'-3' sub-sequence motif(s) of interest, to query whether lost or gained in allele(s))</li>
      </ul>  


## <span style="color:mediumorchid">Output notes</span> 
This script produces **8 output files** in the user-specified output directory.  
These include:  

1. **fasta.fa**  
        (collection of fasta entries representing top 10 most abundant sequences assigned to a single sample ID)  
        
2. **blastn_alignments.txt**  
        (output of blastn operation on fasta.fa)
     
3. **allele\_definitions.txt**  
        (output of script operation on blastn\_alignments.txt, samples returned in order of processing)  
        
4. **allele\_evidence.pdf**  
        (output of script operation on blastn\_alignments.txt, plots of calculated read/allele frequencies)  
        
5. **genotypes.txt**  
        (output of script operation on blastn\_alignments.txt, samples returned in ranked order based on  
        genotype inferences)  
        
6. **allele\_definitions.csv**  
        (tabular representation of allele data for all samples)  

7. **population\_summary.txt**  
        (output of script operation on genotypes.txt)  
        
8. **script\_metrics.txt**  
        (summary/analysis of script operation metrics [metadata])  

  Directory structure under an output directory specified as 'Genotypes', for example,  
           would contain the following files after Genotypes.py operations:  

           /Genotypes  
                          `-----allele_definitions.csv  
                          `-----allele_definitions.txt  
                          `-----allele_evidence.pdf  
                          `-----blastn_alignments.txt  
                          `-----fasta.fa  
                          `-----genotypes.txt  
                          `-----population_summary.txt  
                          `-----script_metrics.txt


## <span style="color:mediumorchid">Visual summary of key script operations</span>  
In short, sequencing data in a sample-specific **fastq file** (*e.g.*, below), are converted to user-interpretable  genotype inferences (**key output files**, below), for 100s to 1000s of samples.    
  
<img src="Genotypes_img/fastq_example.png" align="left" width="700">
<br clear="all" />  
*example*  

------
#### Key output files:  
##### allele_definitions.txt 
Samples are reported with sequence alignments to document alleles, along with extrapolated genotypes. 
<img src="Genotypes_img/genotype_example.png" align="left" width="800">
<br clear="all" />
##### allele_evidence.pdf
Samples are reported with frequency plots as evidence.
<img src="Genotypes_img/frequency_plots_example.png" align="left" width="700">  
<br clear="all" />


## <span style="color:mediumorchid">Status</span>
Project is:  _finished_, _open for contributions_


## <span style="color:mediumorchid">Contact</span>
Created by kirk.ehmsen[at]gmail.com - feel free to contact me!    
Keith Yamamoto laboratory, UCSF, San Francisco, CA.
