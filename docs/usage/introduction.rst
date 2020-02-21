Quick Introduction
==================

The goal of this pipeline is to standardize the data cleaning, featuring extraction, analysis, and evaluation of mobile sensing projects. It leverages Cookiecutter_, Snakemake_, Sphinx_, Scypy_, R_, and Conda_ to create an end-to-end reproducible environment that can be published along with research papers. 

At the moment, mobile data can be collected using different sensing frameworks (AWARE_, Beiwe_) and hardware (Fitbit_). The pipeline is agnostic to these data sources and can unify their analysis. The current implementation only handles data collected with AWARE_. However, it should be easy to extend it to other providers. 

We recommend reading Snakemake_ docs, but the main idea behind the pipeline is that every link in the analysis chain is a rule with an input and an output. Input and output (generally) are files, and these files can be manipulated using any programming language (although Snakemake_ has wrappers for Python_, R_, and Julia_ that can make development slightly more comfortable). Snakemake_ also allows us to spread the execution of rules over multiple cores, which means that a single analysis pipeline can be executed in parallel for all participants in a  study without any code changes.

Available features:

- :ref:`sms` 
- :ref:`calls`
- :ref:`bluetooth`
- :ref:`google-activity-recognition`
- :ref:`battery`
- :ref:`location-features`
- :ref:`screen`
- :ref:`light`
- :ref:`accelerometer`
- :ref:`applications_foreground`
- :ref:`fitbit-heart-rate`
- :ref:`fitbit-steps`

Applications_foreground

We are updating these docs constantly, but if you think something needs clarification, feel free to reach out or submit a pull request in GitHub.


.. _Cookiecutter: http://drivendata.github.io/cookiecutter-data-science/
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _Scypy: https://www.scipy.org/index.html
.. _R: https://www.r-project.org/
.. _Conda: https://docs.conda.io/en/latest/
.. _AWARE: https://awareframework.com/what-is-aware/
.. _Beiwe: https://www.beiwe.org/
.. _Fitbit: https://www.fitbit.com/us/home
.. _Python: https://www.python.org/
.. _Julia: https://julialang.org/