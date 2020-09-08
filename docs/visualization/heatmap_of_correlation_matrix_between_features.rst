.. _heatmap-of-correlation-matrix-between-features

Heatmap of correlation matrix between features
==============================================

Columns and rows are features computed in RAPIDS, cells’ color represents the correlation coefficient between all days of data for every pair of feature of all participants. The user can specify a minimum number of observations required to compute the correlation between two features using the ``MIN_ROWS_RATIO`` parameter (0.5 by default). In addition, this plot can be configured to only display correlation coefficients above a threshold using the ``CORR_THRESHOLD`` parameter (0.1 by default).
