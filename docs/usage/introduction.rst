Quick Introduction
==================

The goal of this pipeline is to standardize the data cleaning, feature extraction, analysis, and evaluation of mobile sensing projects. It leverages Conda_, Cookiecutter_, SciPy_, Snakemake_, Sphinx_, and R_ to create an end-to-end reproducible environment that can be published along with research papers. 

At the moment, mobile data can be collected using different sensing frameworks (AWARE_, Beiwe_) and hardware (Fitbit_). The pipeline is agnostic to these data sources and can unify their analysis. The current implementation only handles data collected with AWARE_. However, it can be easily extended to other providers. 

We recommend reading Snakemake_ docs, but the main idea behind the pipeline is that every link in the analysis chain is a rule with an input and an output. Input and output (generally) are files, which can be manipulated using any programming language (although Snakemake_ has wrappers for Julia_, Python_, and R_ that can make development slightly more comfortable). Snakemake_ also allows the pipeline rules to be executed in parallel on multiple cores without any code changes. This can drastically reduce the time needed to complete an analysis.

Available features:

- :ref:`accelerometer-sensor-doc`
- :ref:`applications-foreground-sensor-doc`
- :ref:`battery-sensor-doc`
- :ref:`bluetooth-sensor-doc`
- :ref:`call-sensor-doc`
- :ref:`fitbit-heart-rate-sensor-doc`
- :ref:`fitbit-steps-sensor-doc`
- :ref:`activity-recognition-sensor-doc`
- :ref:`light-doc`
- :ref:`location-sensor-doc`
- :ref:`screen-sensor-doc`
- :ref:`sms-sensor-doc` 

We are updating these docs constantly, but if you think something needs clarification, feel free to reach out or submit a pull request on GitHub.


.. _Conda: https://docs.conda.io/en/latest/
.. _Cookiecutter: http://drivendata.github.io/cookiecutter-data-science/
.. _SciPy: https://www.scipy.org/index.html
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _R: https://www.r-project.org/

.. _AWARE: https://awareframework.com/what-is-aware/
.. _Beiwe: https://www.beiwe.org/
.. _Fitbit: https://www.fitbit.com/us/home
.. _Python: https://www.python.org/
.. _Julia: https://julialang.org/
