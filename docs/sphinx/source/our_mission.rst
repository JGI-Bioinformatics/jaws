We believe that sharing — of tools, data, and infrastructure — is the fastest way that we, as a community, will be able to unlock profound insights hidden in our genomes.

we didn’t want to build yet another siloed, proprietary, homebrewed solution. We wanted to build a solution that was principled in open standards, reproducibility, and scale; we want these things to empower the global community of bioinformaticians to more easily author, test, run, and share their work so that we can spend more time doing actual science, instead of building infrastructure.


open technologies — WDL and Cromwell — that we are using to deliver Workflows, our managed cloud bioinformatics application.


The Workflow Description Language, or WDL (pronounced “widdle”), is a technical language developed at the Broad Institute that allows workflow developers to describe pipelines in a hardware agnostic way — separating logic from the metal on which it runs. A WDL workflow is composed of a sequence of tasks, each having its own inputs, outputs, runtime environment, and commands. The runtime block outlines the shape of the virtual machine image in which the list of commands will run, producing outputs from inputs. Note that each WDL task prescribes a custom runtime environment, which lends to very lean utilization of computational resources.


Anatomy of a WDL Task
As you can see, WDL is extremely simple to read and to write, and eliminates virtually all coordination logic (e.g. marshalling inputs and outputs) that typically litters workflow execution code. The syntax of WDL implies this coordination, the implementation of which is delegated to the execution engine that is responsible for running it.


DNAstack’s Workflows app provides end-to-end encapsulation of WDL, Cromwell, and Google Cloud Platform functionalities to deliver a streamlined bioinformatics service — workflow authoring, testing, and deployment — that is fully integrated with DNAstack’s data management and visualization applications, and programmable through our REST API, Command Line Interface, and Java Client Library. Workflows has been used to reliably execute thousands of pipelines for diverse use cases by customers all over the world.


WDL workflow authoring tool in DNAstack’s Workflows
We believe that the future of precision medicine will be driven through sequencing and statistical analysis of genomes shared from millions of individuals. However, the lack of consistency and reproducibility with which bioinformatics is currently done poses a serious threat to the utility of data generated between sites. Platform-agnostic workflow languages, like WDL, and workflow execution engines that can interpret them, like Cromwell, elegantly separate computational logic from hardware, allowing the next generation of genomics pipelines to be more easily written, shared, reproduced, and deployed at scale.
