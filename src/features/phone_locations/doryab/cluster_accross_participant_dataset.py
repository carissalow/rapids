import pandas as pd
from doryab_clustering import create_clustering_hyperparameters, cluster



location_data = pd.read_csv(snakemake.input["sensor_input"])
dbscan_eps = snakemake.params["dbscan_eps"]
dbscan_minsamples = snakemake.params["dbscan_minsamples"]
threshold_static = snakemake.params["threshold_static"]
clustering_algorithm = snakemake.params["clustering_algorithm"]

hyperparameters = create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples)
location_data_with_labels = cluster(location_data, clustering_algorithm, threshold_static, **hyperparameters)

location_data_with_labels.to_csv(snakemake.output[0], index=False)
