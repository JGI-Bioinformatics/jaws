#################################
Using Reference Data in Your WDLs
#################################

.. role:: bash(code)
   :language: bash

Reference data are any files that are too large to copy everytime a WDL is run.  It can be something like Uniprot or any custom data you want.  Furthermore, this reference data is guarranteed to be the same on all the available JAWS sites.  

Using Existing Reference Data
-----------------------------
First you need to know what's available.  The best way is to log into CORI and look at :bash:`/global/dna/shared/databases/jaws/refdata`.  Furthermore, there should be a README file in there that describes a little about the provenance of each database.  

Now let's say you want to run blastn. The required nt database can be accessed within your WDL in two ways.  

1) If you are running blast inside a docker container, your command would look like :bash:`blastn --db /refdata/20201117/nt/nt --query fasta.fa --out results.out`.  

2) If you are not running inside a container, you would have to point to the CORI filesystem: :bash:`/global/dna/shared/databases/jaws/refdata/20201117/nt/nt` but now your WDL is not portable to run at different sites.

.. note::

	One of the benifits of running inside a docker container is that the WDL is portable since :bash:`/refdata` points to the appropriate data repository for that site.


Adding Something to the Repository
----------------------------------
If your required reference data is not already in the repository, you will need to contact a JAWS administrator by one of the following methods. They will then add the data to the CORI's repository which will be automatically sync'd up with the other sites. You should recieve a reply when everything is complete.

* email: jaws-support@lbl.gov 
* slack: #jaws channel 
* JIRA:  https://intranet.lbl.gov/jgi/services/computers-networking/jaws/

Requirements:

* No symlinks (e.g. latest -> v10.4). Symlinks will not be maintained when the data files are sync'd between sites.
