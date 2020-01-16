Quick Introduction
==================

The goal of this pipeline is to standardize the data cleaning, featuring extraction, analysis, and evaluation of mobile sensing projects. It leverages Cookiecutter, Snakemake, Sphinx, Scypy, R, and Conda to create an end-to-end reproducible environment that can be published along with research papers. 

At the moment, mobile data can be collected using different sensing frameworks (Aware, Beiwe) and hardware (Fitbit). The pipeline is agnostic to these data sources and can unify their analysis. The current implementation only handles data collected with Aware. However, it should be easy to extend it to other providers. 

We recommend reading Snakemake docs, but the main idea behind the pipeline is that every link in the analysis chain is a rule with an input and an output. Input and output (generally) are files, and these files can be manipulated using any programming language (although Snakemake has wrappers for Python, R, and Julia that can make development slightly more comfortable). Snakemake also allows us to spread the execution of rules over multiple cores, which means that a single analysis pipeline can be executed in parallel for all participants in a  study without any code changes.

Available features:

- SMS
- Calls
- Bluetooth
- Google Activity Recognition
- Battery
- Location (Barnett's)
- Screen
- Light
- Accelerometer

We are updating these docs constantly, but if you think something needs clarification, feel free to reach out or submit a pull request in GitHub.
