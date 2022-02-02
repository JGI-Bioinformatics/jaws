#################################
Using Reference Data in Your WDLs
#################################

.. role:: bash(code)
   :language: bash

Reference data are too large to copy everytime a WDL is run.  You can add reference files to your own folder under the jaws refdata folder. It is guarranteed to be synchronized on all the JAWS compute sites.  

How to Use JAWS Refdata
-----------------------
Everyone with a jaws token should have their own directory under a jaws owned folder. For example, log into CORI and look for your username at :bash:`/global/dna/shared/databases/jaws/refdata`. If you can't find your username, contact jaws-support@lbl.gov.

Once you add files to your directory, a WDL task that is running through a docker container can access them through :bash:`/refdata/<your-user-name>`.  If you were not running tasks through a container, you'd have to use the cori path :bash:`/global/dna/shared/databases/jaws/refdata` and your WDL would not be portable.

For example, let's say you want to run blastn. The required nt database can be accessed within your WDL in two ways.  

1) If you are running blast inside a docker container, your command would look like :bash:`blastn --db /refdata/jfroula/NT/20201117/nt/nt --query fasta.fa --out results.out`.  

2) If you are not running inside a container, you would have to point to the CORI filesystem: :bash:`/global/dna/shared/databases/jaws/refdata/jfroula/NT/20201117/nt/nt` (not portable between compute sites).

.. note::

	Syncing the refdata folder from cori to all the other compute sites happens at 1AM and 1PM.


Adding Something to the Repository
----------------------------------

Requirements:

* No symlinks (e.g. latest -> v10.4). Symlinks will not be maintained when the data files are sync'd between sites.



for questions:

* email: jaws-support@lbl.gov 
* slack: #jaws channel 
* JIRA:  https://intranet.lbl.gov/jgi/services/computers-networking/jaws/

