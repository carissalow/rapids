.. _heatmap-of-phone-sensed-bins:

Heatmap of phone sensed bins
============================

See `Heatmap of Phone Sensed Bins Config Code`_

**Rule Chain:**

- Rule: ``rules/preprocessing.smk/download_dataset``
- Rule: ``rules/preprocessing.smk/readable_datetime``
- Rule: ``rules/preprocessing.smk/phone_sensed_bins``
- Rule: ``rules/reports.smk/heatmap_sensed_bins``

.. _figure2-parameters:

**Parameters of heatmap_sensed_bins Rule:**

=======================    =======================
Name                       Description
=======================    =======================
plot                       Whether the rule is executed or not. The available options are ``True`` and ``False``.
bin_size                   Every hour is divided into N bins of size ``BIN_SIZE`` (in minutes).
=======================    =======================

**Observations:**

In this heatmap rows are dates, columns are sensed bins for a participant, and cells’ color shows the number of mobile sensors that logged at least one row of data during that bin. This plot shows the periods of time without data for a participant and can be used as a rough indication of whether time-based sensors were following their sensing schedule (e.g. if location was being sensed every 2 minutes). See Figure 2.

.. figure:: figures/Figure2.png
    :scale: 90 %
    :align: center

    Figure 2 Heatmap of phone sensed bins for a single participant


.. _`Heatmap of Phone Sensed Bins Config Code`: https://github.com/carissalow/rapids/blob/master/config.yaml#L233
